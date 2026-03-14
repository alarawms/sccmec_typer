# SCCmec Typer: SVG Visualization & Confidence Scoring Design

## Goal

Add IWG-SCC-style color-coded SVG visualizations, an interactive HTML report, and a three-level confidence scoring system to the SCCmec typer. All data must also be present in JSON output for integration with the Staphit master analysis pipeline.

## Architecture

The visualization system is built on Jinja2 templates rendered by a Python visualizer module. A separate confidence scoring module computes scores at three levels (per-gene, per-component, per-type) and injects them into the classifier result dict before visualization. This separation ensures confidence data is available in JSON even when SVG generation is skipped.

## Dependencies

- `jinja2` — template rendering for SVG and HTML

---

## 1. Color Scheme (IWG-SCC Faithful)

| Element | Color | Hex |
|---|---|---|
| mecA / mecC / mecB | Red | `#DC143C` |
| mecR1, mecI | Pink | `#E8A0BF` |
| ccr genes (ccrA, ccrB, ccrC) | Blue-Violet | `#7B68EE` |
| IS elements (IS431, IS1272) | Grey | `#95A5A6` |
| blaZ / accessory genes | Purple | `#8E44AD` |
| J regions (background) | Light grey | `#F0F0F0` |
| DR boundaries | Red arrowheads | `#DC143C` |
| Expected but missing gene | Dashed outline, no fill | — |

Confidence overlay:
- Gene arrow opacity scaled to per-gene confidence score (0.4 min opacity for visibility)
- Confidence tier badge colors: High `#27AE60`, Medium `#F39C12`, Low `#E74C3C`

---

## 2. Confidence Scoring Model

### Level 1: Per-Gene Confidence (0.0–1.0)

```
norm_id  = clamp((id_pct - 90) / 10,  0.0, 1.0)
norm_cov = clamp((cov_pct - 80) / 20, 0.0, 1.0)
gene_confidence = 0.6 * norm_id + 0.4 * norm_cov
```

- Both normalized terms MUST be clamped to [0.0, 1.0] before weighting
- Maps 90–100% identity and 80–100% coverage to 0.0–1.0 each
- Identity weighted 0.6, coverage weighted 0.4
- Tiers: High >= 0.75, Medium 0.40–0.75, Low < 0.40

#### Reads Mode Fallback

In reads mode, `coverage.py` does not compute per-gene identity (`id_pct = 0`). When `id_pct` is 0 or unavailable, confidence uses **coverage-only mode**:

```
gene_confidence = clamp((cov_pct - 80) / 20, 0.0, 1.0)
```

The JSON output sets `confidence.mode` to `"coverage_only"` and the visualizer displays a note: "Confidence based on coverage only (reads mode)." Tier thresholds remain the same.

### Level 2: Component Confidence

```
component_confidence = (genes_detected / genes_expected) * mean(gene_confidences)
```

- Applied separately to mec complex and ccr complex
- IS431 orientation uncertainty penalty: 0.8x when C1/C2 cannot be resolved
  - This penalty applies ONLY to mec classes containing "Class C" in the name (Class C, Class C1, Class C2)
  - It does NOT apply to Class E (mecC-based) despite the old "Class C2" naming

#### ccr = "ANY" case (Plasmid mecB, or any rule with ccr: "ANY")

When the matched SCCmec type rule specifies `ccr: "ANY"`, no ccr complex is expected. In this case:
- `ccr_confidence = 1.0` (neutral — ccr does not factor into type confidence)
- `confidence.ccr_component` in JSON is set to `{"type": "n/a", "score": 1.0, "tier": "High", "genes_expected": [], "genes_detected": [], "genes_missing": [], "note": "No ccr expected for this type"}`

### Level 3: Type-Level Confidence

```
type_confidence = mec_confidence * ccr_confidence * assembly_penalty
```

- `ccr_confidence`: 1.0 when ccr is "ANY" (see above), otherwise from Level 2
- `assembly_penalty`: 1.0 (single contig) or 0.85 (split assembly)
- Tiers: High > 0.70, Medium 0.40–0.70, Low < 0.40

Full methods documented in `docs/confidence_scoring_methods.md`.

---

## 3. Visualization Components

### Diagram A: Schematic Cassette Diagram

Standardized horizontal layout showing the canonical SCCmec structure for the assigned type:

```
┌─ DR ─┐                                                      ┌─ DR ─┐
│  ▶   │  J3  │ ccr complex │  J2  │  mec complex  │  J1  │  ▶   │
│      │      │  ccrA2 ccrB2│      │ IS431-mecA-R1-I│      │      │
└──────┘      │  [violet]   │      │ [red] [pink]   │      └──────┘
              │  conf: 0.94 │      │ conf: 0.99     │
```

- Fixed left-to-right: DR → J3 → ccr → J2 → mec → J1 → DR
- Detected genes: solid fill, opacity scaled to gene confidence
- Missing expected genes: dashed outline, no fill
- Contig origin badge on each gene when split assembly (`[c1]`, `[c2]`)
- Split assembly: dashed separator line + warning banner
- Type assignment banner at top with overall confidence tier and score
- Component confidence badges below each complex

**Negative/Partial results:** When status is "Negative", the schematic shows a "No SCCmec detected" banner with an empty grey cassette outline and no gene arrows. For "Partial" results, only detected components are drawn; missing components show dashed placeholders with "Not detected" labels.

### Diagram B: Coordinate Map (assembly mode only)

Genes at actual genomic positions, one horizontal track per contig:

```
contig_1 (85,000 bp)
├────────────────────────────────────────────────────┤
   ▶ mecA     ▶ mecR1    ▶ IS1272    ▶ IS431
   28,100     29,500     31,200      33,400

contig_2 (42,000 bp)
├────────────────────────────────────────────────────┤
   ◀ ccrB1              ◀ ccrA1
   8,200                 10,500

── ── ── ── ── ── ── ── ── ── ── ── ── ── ── ── ──
⚠ Split Assembly: mec on contig_1, ccr on contig_2
```

- Arrow direction indicates strand
- Proportional gaps between genes
- Stacked contig tracks for split assemblies with dashed connector
- **Skipped in reads mode** (no coordinate data available)
- **Skipped for Negative results** (no genes to map)

### HTML Report

Single-file HTML with embedded SVGs, gene details table, and warnings:

```
┌──────────────────────────────────────────────────┐
│  SCCmec Typing Report: {sample_name}             │
├──────────────────────────────────────────────────┤
│  Type: IV (2B)     Confidence: HIGH (0.92)       │
│  mec: Class B      ccr: Type 2 (ccrA2B2)        │
├──────────────────────────────────────────────────┤
│  [Schematic Cassette Diagram SVG]                │
├──────────────────────────────────────────────────┤
│  [Coordinate Map SVG] (assembly mode only)       │
├──────────────────────────────────────────────────┤
│  Gene Details Table                              │
│  ┌────────┬───────┬──────┬────────┬──────────┐   │
│  │ Gene   │ ID%   │ Cov% │ Strand │ Conf     │   │
│  ├────────┼───────┼──────┼────────┼──────────┤   │
│  │ mecA   │ 99.8  │ 100  │  +     │ 1.0 High │   │
│  │ ...    │       │      │        │          │   │
│  └────────┴───────┴──────┴────────┴──────────┘   │
├──────────────────────────────────────────────────┤
│  Warnings (if any)                               │
└──────────────────────────────────────────────────┘
```

- Hover tooltips on gene arrows show full alignment details
- Self-contained single file (CSS embedded, no external assets)
- For Negative results: shows "No SCCmec detected" with empty cassette outline, no gene table

---

## 4. Enhanced JSON Schema (Staphit Integration)

### 4a. Positive Result Schema

```json
{
    "status": "Positive",
    "sccmec_type": "Type IV",
    "designation": "2B",
    "mec_complex": "Class B",
    "ccr_complex": "Type 2",
    "ccr_genes": ["ccrA2", "ccrB2"],
    "genes_detected": ["IS1272", "IS431", "ccrA2", "ccrB2", "mecA", "mecR1"],
    "mecA_present": true,
    "mecC_present": false,
    "mecB_present": false,
    "confidence": {
        "mode": "full",
        "type_level": {
            "score": 0.923,
            "tier": "High"
        },
        "mec_component": {
            "class": "Class B",
            "score": 0.991,
            "tier": "High",
            "genes_expected": ["mecA", "IS1272"],
            "genes_detected": ["mecA", "IS1272"],
            "genes_missing": [],
            "orientation": null,
            "orientation_resolved": null
        },
        "ccr_component": {
            "type": "Type 2",
            "score": 0.931,
            "tier": "High",
            "genes_expected": ["ccrA2", "ccrB2"],
            "genes_detected": ["ccrA2", "ccrB2"],
            "genes_missing": []
        },
        "per_gene": {
            "mecA": {"identity_pct": 100.0, "coverage_pct": 100.0, "score": 1.0, "tier": "High"},
            "ccrA2": {"identity_pct": 99.7, "coverage_pct": 98.0, "score": 0.942, "tier": "High"}
        }
    },
    "assembly": {
        "contigs": ["contig_1"],
        "is_split": false,
        "split_penalty_applied": false,
        "gene_locations": [
            {"gene": "mecA", "contig": "contig_1", "start": 28100, "end": 30100, "strand": "+"}
        ]
    },
    "warnings": [],
    "hits_summary": [...]
}
```

### 4b. Negative Result Schema

```json
{
    "status": "Negative",
    "reason": "No alignments found",
    "sccmec_type": "Negative",
    "designation": null,
    "mec_complex": "Negative",
    "ccr_complex": "Negative",
    "ccr_genes": [],
    "genes_detected": [],
    "mecA_present": false,
    "mecC_present": false,
    "mecB_present": false,
    "confidence": {
        "mode": "n/a",
        "type_level": {"score": 0.0, "tier": "n/a"},
        "mec_component": null,
        "ccr_component": null,
        "per_gene": {}
    },
    "assembly": {
        "contigs": [],
        "is_split": false,
        "split_penalty_applied": false,
        "gene_locations": []
    },
    "warnings": [],
    "hits_summary": []
}
```

### 4c. Partial Result Schema

Same structure as Positive but with:
- `status`: `"Partial (Orphan ccr)"` or `"Partial (Unclassifiable)"`
- `sccmec_type`: `"Unknown"`
- `designation`: `null`
- Confidence scores reflect missing components (e.g., `mec_component: null` for orphan ccr)

### 4d. Reads Mode Assembly Block

In reads mode, the `assembly` block uses sentinel values:

```json
{
    "assembly": {
        "contigs": ["Reads"],
        "is_split": false,
        "split_penalty_applied": false,
        "gene_locations": []
    }
}
```

`gene_locations` is empty because coordinate data is not available from reads-mode alignment. The `confidence.mode` field is set to `"coverage_only"`.

### 4e. Designation Field Source

The `designation` field (e.g., `"2B"` for Type IV) is sourced from the `"designation"` key in each `sccmec_type_rules` entry in `rules.json`. The classifier copies it to the result dict when a type is matched. For unmatched or negative results, `designation` is `null`.

Key fields for Staphit:
- `confidence.type_level.score` and `.tier` — for filtering/QC in master analysis
- `confidence.mode` — indicates whether confidence is full or coverage-only
- `assembly.is_split` — flag for assembly quality assessment
- `designation` — compact IWG-SCC designation for tabular reports
- `ccr_genes` — explicit gene list for cross-referencing
- `confidence.*.genes_missing` — identifies incomplete elements

---

## 5. File Architecture

| File | Action | Responsibility |
|---|---|---|
| `lib/confidence.py` | Create | Three-level confidence scoring. Pure calculation, no rendering. Handles reads-mode fallback and ccr="ANY" cases. |
| `lib/visualizer.py` | Create | SVG generation and HTML report assembly. Loads Jinja2 templates, computes gene positions, calls confidence module. |
| `templates/cassette.svg.j2` | Create | Jinja2 template for schematic cassette diagram |
| `templates/coordinate_map.svg.j2` | Create | Jinja2 template for genomic coordinate map |
| `templates/report.html.j2` | Create | Jinja2 template for HTML report (embeds SVGs + gene table + confidence) |
| `lib/classifier.py` | Modify | Call confidence scoring, add confidence/assembly/designation data to result dict. Ensure all return paths produce consistent schema. |
| `bin/sccmec_typer.py` | Modify | Call visualizer, write `_map.svg` and `_report.html` output files. Add `--no-viz` flag. |
| `environment.yml` | Modify | Add jinja2 dependency |

Data flow:

```
classifier.py result (consistent schema for all statuses)
       │
       ▼
confidence.py ──▶ enriches result with confidence scores
       │           (handles reads-mode, negative, partial, ccr=ANY)
       ▼
visualizer.py
       ├──▶ cassette.svg.j2          ──▶ {prefix}_map.svg
       ├──▶ coordinate_map.svg.j2         (assembly only, also embedded in HTML)
       └──▶ report.html.j2           ──▶ {prefix}_report.html
```

---

## 6. Output Files

Per sample, the tool produces:

| File | Format | Content |
|---|---|---|
| `{prefix}.json` | JSON | Full classification + confidence + assembly data (all statuses) |
| `{prefix}.tsv` | TSV | One-line summary |
| `{prefix}_elements.csv` | CSV | Per-gene alignment details |
| `{prefix}_map.svg` | SVG | Schematic cassette diagram (standalone) |
| `{prefix}_report.html` | HTML | Interactive report with both diagrams + gene table |

- Optional flag `--no-viz` to skip SVG/HTML generation for pipeline-only use.
- The coordinate map SVG is embedded within the HTML report only. If a standalone coordinate map SVG is needed in the future, a `{prefix}_coordinate_map.svg` output can be added.

---

## 7. Contig / Split Assembly Handling

- **Single contig**: Clean layout, no badges, no penalty
- **Multiple contigs**:
  - Schematic diagram: contig origin badge `[c1]`, `[c2]` on each gene, dashed separator
  - Coordinate map: stacked tracks, one per contig, dashed connector lines
  - Confidence: 0.85 assembly penalty at type level
  - JSON: `assembly.is_split = true`, full `gene_locations` with contig attribution
  - Warning in both visual and JSON output
- **Reads mode**: `assembly.contigs = ["Reads"]`, `is_split = false`, `gene_locations = []`, coordinate map skipped
