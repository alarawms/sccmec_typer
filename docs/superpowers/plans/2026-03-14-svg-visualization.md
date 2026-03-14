# SVG Visualization & Confidence Scoring Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add IWG-SCC-style color-coded SVG visualization, interactive HTML reports, and three-level confidence scoring to the SCCmec typer, with all data in JSON for Staphit pipeline integration.

**Architecture:** A `confidence.py` module scores at gene/component/type levels and enriches the classifier result dict. A `visualizer.py` module renders Jinja2 SVG/HTML templates using the enriched data. The classifier is updated to produce a consistent JSON schema across all result statuses (Positive, Negative, Partial). All confidence and assembly data is available in JSON even when visualization is skipped (`--no-viz`).

**Tech Stack:** Python 3.9+, Jinja2, pytest

**Depends on:** The IWG-SCC Classification Update plan must be executed first (rules.json + classifier rewrite).

**Spec:** `docs/superpowers/specs/2026-03-14-svg-visualization-design.md`

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `lib/confidence.py` | **Create** | Three-level confidence scoring (gene, component, type) |
| `lib/visualizer.py` | **Create** | SVG/HTML generation via Jinja2 templates |
| `templates/cassette.svg.j2` | **Create** | Schematic cassette diagram template |
| `templates/coordinate_map.svg.j2` | **Create** | Genomic coordinate map template |
| `templates/report.html.j2` | **Create** | HTML report template (embeds SVGs + table) |
| `lib/classifier.py` | **Modify** | Consistent JSON schema, designation field, assembly block |
| `bin/sccmec_typer.py` | **Modify** | Call confidence + visualizer, `--no-viz` flag |
| `environment.yml` | **Modify** | Add jinja2 dependency |
| `tests/test_confidence.py` | **Create** | Unit tests for confidence scoring |
| `tests/test_visualizer.py` | **Create** | Unit tests for SVG/HTML generation |

---

## Chunk 1: Confidence Scoring Module

### Task 1: Add jinja2 dependency

**Files:**
- Modify: `environment.yml`

- [ ] **Step 1: Add jinja2 to environment.yml**

Add `- jinja2` under the dependencies section.

- [ ] **Step 2: Install dependency**

Run: `pip install jinja2`

- [ ] **Step 3: Commit**

```bash
git add environment.yml
git commit -m "chore: add jinja2 dependency for SVG visualization"
```

### Task 2: Create confidence scoring module with tests

**Files:**
- Create: `lib/confidence.py`
- Create: `tests/test_confidence.py`

- [ ] **Step 1: Write tests/test_confidence.py**

```python
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.confidence import (
    compute_gene_confidence,
    compute_component_confidence,
    compute_type_confidence,
    enrich_result_with_confidence,
)


class TestGeneConfidence:
    """Level 1: Per-gene confidence scoring."""

    def test_perfect_hit(self):
        score, tier = compute_gene_confidence(100.0, 100.0)
        assert score == 1.0
        assert tier == "High"

    def test_minimum_threshold(self):
        score, tier = compute_gene_confidence(90.0, 80.0)
        assert score == 0.0
        assert tier == "Low"

    def test_mid_range(self):
        score, tier = compute_gene_confidence(95.0, 90.0)
        # norm_id = (95-90)/10 = 0.5, norm_cov = (90-80)/20 = 0.5
        # score = 0.6*0.5 + 0.4*0.5 = 0.5
        assert abs(score - 0.5) < 0.001
        assert tier == "Medium"

    def test_high_threshold(self):
        score, tier = compute_gene_confidence(98.0, 96.0)
        # norm_id = 0.8, norm_cov = 0.8
        # score = 0.6*0.8 + 0.4*0.8 = 0.8
        assert abs(score - 0.8) < 0.001
        assert tier == "High"

    def test_clamp_above_100(self):
        """PAF can produce identity > 100% in edge cases."""
        score, tier = compute_gene_confidence(102.0, 105.0)
        assert score == 1.0

    def test_reads_mode_fallback(self):
        """When id_pct is 0 (reads mode), use coverage-only."""
        score, tier = compute_gene_confidence(0.0, 95.0)
        # coverage-only: (95-80)/20 = 0.75
        assert abs(score - 0.75) < 0.001
        assert tier == "High"

    def test_reads_mode_low(self):
        score, tier = compute_gene_confidence(0.0, 82.0)
        # (82-80)/20 = 0.1
        assert abs(score - 0.1) < 0.001
        assert tier == "Low"


class TestComponentConfidence:
    """Level 2: Component-level confidence."""

    def test_all_genes_present_high(self):
        gene_scores = {"mecA": 1.0, "mecR1": 0.98, "mecI": 0.99}
        expected_genes = ["mecA", "mecR1", "mecI"]
        score, tier, missing = compute_component_confidence(
            gene_scores, expected_genes
        )
        assert score > 0.9
        assert tier == "High"
        assert missing == []

    def test_missing_gene(self):
        gene_scores = {"mecA": 1.0, "mecR1": 0.98}
        expected_genes = ["mecA", "mecR1", "mecI"]
        score, tier, missing = compute_component_confidence(
            gene_scores, expected_genes
        )
        # (2/3) * mean(1.0, 0.98) = 0.667 * 0.99 = 0.660
        assert score < 0.7
        assert "mecI" in missing

    def test_orientation_penalty(self):
        gene_scores = {"mecA": 1.0}
        expected_genes = ["mecA"]
        score_no_penalty, _, _ = compute_component_confidence(
            gene_scores, expected_genes
        )
        score_penalty, _, _ = compute_component_confidence(
            gene_scores, expected_genes, orientation_undetermined=True,
            mec_class="Class C"
        )
        assert abs(score_penalty - score_no_penalty * 0.8) < 0.001

    def test_orientation_penalty_not_applied_to_class_E(self):
        gene_scores = {"mecC": 1.0}
        expected_genes = ["mecC"]
        score, _, _ = compute_component_confidence(
            gene_scores, expected_genes, orientation_undetermined=True,
            mec_class="Class E"
        )
        # No penalty for Class E
        assert score == 1.0

    def test_empty_expected(self):
        """ccr=ANY case: no genes expected."""
        score, tier, missing = compute_component_confidence({}, [])
        assert score == 1.0
        assert tier == "High"


class TestTypeConfidence:
    """Level 3: Type-level confidence."""

    def test_high_confidence(self):
        score, tier = compute_type_confidence(0.99, 0.93, is_split=False)
        assert score > 0.7
        assert tier == "High"

    def test_split_penalty(self):
        score_normal, _ = compute_type_confidence(0.99, 0.93, is_split=False)
        score_split, _ = compute_type_confidence(0.99, 0.93, is_split=True)
        assert abs(score_split - score_normal * 0.85) < 0.001

    def test_low_confidence(self):
        score, tier = compute_type_confidence(0.2, 0.3, is_split=True)
        assert score < 0.4
        assert tier == "Low"


class TestEnrichResult:
    """Integration: enrich_result_with_confidence adds all fields."""

    def test_positive_result_has_confidence(self):
        result = {
            "status": "Positive",
            "sccmec_type": "Type II",
            "mec_complex": "Class A",
            "ccr_complex": "Type 2",
            "genes_detected": ["mecA", "mecR1", "mecI", "ccrA2", "ccrB2"],
            "mecA_present": True,
            "warnings": [],
            "hits_summary": [
                {"gene": "mecA", "id_pct": 100, "cov_pct": 100, "contig": "c1",
                 "start": 0, "end": 1000, "strand": "+"},
                {"gene": "mecR1", "id_pct": 99, "cov_pct": 99, "contig": "c1",
                 "start": 1000, "end": 2000, "strand": "+"},
                {"gene": "mecI", "id_pct": 99, "cov_pct": 100, "contig": "c1",
                 "start": 2000, "end": 3000, "strand": "+"},
                {"gene": "ccrA2", "id_pct": 99, "cov_pct": 98, "contig": "c1",
                 "start": 5000, "end": 6000, "strand": "-"},
                {"gene": "ccrB2", "id_pct": 99, "cov_pct": 97, "contig": "c1",
                 "start": 6000, "end": 7000, "strand": "-"},
            ],
        }
        enriched = enrich_result_with_confidence(result)
        assert "confidence" in enriched
        assert "assembly" in enriched
        assert enriched["confidence"]["mode"] == "full"
        assert enriched["confidence"]["type_level"]["tier"] == "High"
        assert "mecA" in enriched["confidence"]["per_gene"]
        assert enriched["assembly"]["is_split"] is False

    def test_negative_result_has_confidence(self):
        result = {
            "status": "Negative",
            "reason": "No alignments found",
        }
        enriched = enrich_result_with_confidence(result)
        assert enriched["confidence"]["mode"] == "n/a"
        assert enriched["confidence"]["type_level"]["score"] == 0.0
        assert enriched["assembly"]["gene_locations"] == []
```

- [ ] **Step 2: Run tests to verify they FAIL**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/test_confidence.py -v 2>&1 | head -40`
Expected: ImportError (module doesn't exist yet)

- [ ] **Step 3: Create lib/confidence.py**

```python
import json
import os


def _load_rules():
    """Load classification rules for expected gene lookups."""
    rules_path = os.path.join(os.path.dirname(__file__), "../db/rules.json")
    with open(rules_path, "r") as f:
        return json.load(f)


def _clamp(val, lo, hi):
    return max(lo, min(hi, val))


def _tier(score):
    if score >= 0.75:
        return "High"
    elif score >= 0.40:
        return "Medium"
    else:
        return "Low"


def _type_tier(score):
    if score > 0.70:
        return "High"
    elif score >= 0.40:
        return "Medium"
    else:
        return "Low"


def compute_gene_confidence(id_pct, cov_pct):
    """Level 1: Per-gene confidence from identity and coverage percentages.

    When id_pct is 0 (reads mode), uses coverage-only scoring.

    Returns (score, tier) where score is 0.0-1.0.
    """
    norm_cov = _clamp((cov_pct - 80) / 20, 0.0, 1.0)

    if id_pct == 0 or id_pct is None:
        # Reads mode fallback: coverage only
        return norm_cov, _tier(norm_cov)

    norm_id = _clamp((id_pct - 90) / 10, 0.0, 1.0)
    score = 0.6 * norm_id + 0.4 * norm_cov
    return score, _tier(score)


def compute_component_confidence(gene_scores, expected_genes,
                                  orientation_undetermined=False,
                                  mec_class=None):
    """Level 2: Component confidence from gene completeness and quality.

    Args:
        gene_scores: dict of {gene_name: confidence_score} for detected genes
        expected_genes: list of gene names expected for this component
        orientation_undetermined: whether IS431 orientation could not be resolved
        mec_class: the mec complex class name (for orientation penalty logic)

    Returns (score, tier, missing_genes).
    """
    if not expected_genes:
        # ccr = "ANY" case: no genes expected, neutral confidence
        return 1.0, "High", []

    detected = [g for g in expected_genes if g in gene_scores]
    missing = [g for g in expected_genes if g not in gene_scores]

    if not detected:
        return 0.0, "Low", missing

    completeness = len(detected) / len(expected_genes)
    mean_score = sum(gene_scores[g] for g in detected) / len(detected)
    score = completeness * mean_score

    # IS431 orientation penalty: only for Class C variants
    if orientation_undetermined and mec_class and mec_class.startswith("Class C"):
        score *= 0.80

    return score, _tier(score), missing


def compute_type_confidence(mec_confidence, ccr_confidence, is_split=False):
    """Level 3: Overall type confidence combining components and assembly quality.

    Returns (score, tier).
    """
    penalty = 0.85 if is_split else 1.0
    score = mec_confidence * ccr_confidence * penalty
    return score, _type_tier(score)


def _get_expected_genes_for_mec(mec_complex, rules):
    """Look up which genes are expected for a given mec complex class."""
    for rule in rules["mec_complex_rules"]:
        if rule["name"] == mec_complex:
            genes = list(rule["required"])
            if "any_of" in rule:
                # Add the any_of genes as expected (at least one should be present)
                genes.extend(rule["any_of"])
            return genes
    return []


def _get_expected_genes_for_ccr(ccr_type_name, rules):
    """Look up which genes are expected for a given ccr complex type."""
    for rule in rules["ccr_complex_rules"]:
        if rule["name"] == ccr_type_name:
            return list(rule["required_genes"])
    return []


def _get_designation(sccmec_type, rules):
    """Look up the IWG-SCC designation for an SCCmec type."""
    for rule in rules["sccmec_type_rules"]:
        if rule["name"] == sccmec_type:
            return rule.get("designation")
    return None


def _is_ccr_any(sccmec_type, rules):
    """Check if an SCCmec type rule has ccr: ANY."""
    for rule in rules["sccmec_type_rules"]:
        if rule["name"] == sccmec_type:
            return rule.get("ccr") == "ANY"
    return False


def enrich_result_with_confidence(result):
    """Add confidence scores, assembly block, and designation to a classifier result.

    This function handles all result statuses (Positive, Negative, Partial, Error)
    and produces a consistent JSON schema.

    Modifies and returns the result dict in-place.
    """
    # Handle Negative/Error results with minimal schema
    if result.get("status") in ("Negative", "Error"):
        result.setdefault("sccmec_type", "Negative")
        result.setdefault("designation", None)
        result.setdefault("mec_complex", "Negative")
        result.setdefault("ccr_complex", "Negative")
        result.setdefault("ccr_genes", [])
        result.setdefault("genes_detected", [])
        result.setdefault("mecA_present", False)
        result.setdefault("mecC_present", False)
        result.setdefault("mecB_present", False)
        result.setdefault("warnings", [])
        result.setdefault("hits_summary", [])
        result["confidence"] = {
            "mode": "n/a",
            "type_level": {"score": 0.0, "tier": "n/a"},
            "mec_component": None,
            "ccr_component": None,
            "per_gene": {},
        }
        result["assembly"] = {
            "contigs": [],
            "is_split": False,
            "split_penalty_applied": False,
            "gene_locations": [],
        }
        return result

    rules = _load_rules()
    hits = result.get("hits_summary", [])

    # Detect reads mode: all hits have id_pct == 0 or contig == "Reads"
    is_reads_mode = all(
        h.get("id_pct", 0) == 0 or h.get("contig") == "Reads" for h in hits
    ) if hits else False

    # --- Per-gene confidence ---
    per_gene = {}
    for h in hits:
        gene = h["gene"]
        id_pct = h.get("id_pct", 0)
        cov_pct = h.get("cov_pct", 0)
        score, tier = compute_gene_confidence(id_pct, cov_pct)
        # Keep best score if gene appears multiple times
        if gene not in per_gene or score > per_gene[gene]["score"]:
            per_gene[gene] = {
                "identity_pct": id_pct,
                "coverage_pct": cov_pct,
                "score": round(score, 3),
                "tier": tier,
            }

    gene_score_map = {g: v["score"] for g, v in per_gene.items()}

    # --- mec component confidence ---
    mec_complex = result.get("mec_complex", "Negative")
    mec_component = None
    mec_conf_score = 0.0

    if mec_complex != "Negative":
        mec_expected = _get_expected_genes_for_mec(mec_complex, rules)
        # For Class C rules, expected genes include any_of; only count present ones
        # Filter expected to just the required genes for scoring
        for rule in rules["mec_complex_rules"]:
            if rule["name"] == mec_complex:
                core_expected = list(rule["required"])
                break
        else:
            core_expected = mec_expected

        orientation_undetermined = mec_complex == "Class C"
        mec_conf_score, mec_tier, mec_missing = compute_component_confidence(
            gene_score_map, core_expected,
            orientation_undetermined=orientation_undetermined,
            mec_class=mec_complex,
        )
        mec_component = {
            "class": mec_complex,
            "score": round(mec_conf_score, 3),
            "tier": mec_tier,
            "genes_expected": core_expected,
            "genes_detected": [g for g in core_expected if g in gene_score_map],
            "genes_missing": mec_missing,
            "orientation": "undetermined" if orientation_undetermined else None,
            "orientation_resolved": not orientation_undetermined if mec_complex.startswith("Class C") else None,
        }

    # --- ccr component confidence ---
    sccmec_type = result.get("sccmec_type", "Unknown")
    ccr_any = _is_ccr_any(sccmec_type, rules)
    ccr_complex = result.get("ccr_complex", "Negative")
    ccr_component = None
    ccr_conf_score = 1.0  # default neutral

    if ccr_any:
        ccr_component = {
            "type": "n/a",
            "score": 1.0,
            "tier": "High",
            "genes_expected": [],
            "genes_detected": [],
            "genes_missing": [],
            "note": "No ccr expected for this type",
        }
    elif ccr_complex != "Negative":
        # Get the primary ccr type (first if multiple)
        primary_ccr = ccr_complex.split(" / ")[0] if "/" in ccr_complex else ccr_complex
        ccr_expected = _get_expected_genes_for_ccr(primary_ccr, rules)
        ccr_conf_score, ccr_tier, ccr_missing = compute_component_confidence(
            gene_score_map, ccr_expected
        )
        ccr_component = {
            "type": primary_ccr,
            "score": round(ccr_conf_score, 3),
            "tier": ccr_tier,
            "genes_expected": ccr_expected,
            "genes_detected": [g for g in ccr_expected if g in gene_score_map],
            "genes_missing": ccr_missing,
        }

    # --- Type-level confidence ---
    contigs = list(set(h.get("contig", "unknown") for h in hits))
    is_split = len(contigs) > 1 and not is_reads_mode

    type_score, type_tier = compute_type_confidence(
        mec_conf_score, ccr_conf_score, is_split=is_split
    )

    # --- Assembly block ---
    gene_locations = []
    if not is_reads_mode:
        for h in hits:
            gene_locations.append({
                "gene": h["gene"],
                "contig": h.get("contig", "unknown"),
                "start": h.get("start", 0),
                "end": h.get("end", 0),
                "strand": h.get("strand", "+"),
            })

    # --- ccr_genes list ---
    ccr_genes = sorted([g for g in result.get("genes_detected", []) if g.startswith("ccr")])

    # --- Enrich result ---
    result["designation"] = _get_designation(sccmec_type, rules)
    result["ccr_genes"] = ccr_genes
    result["mecC_present"] = "mecC" in result.get("genes_detected", [])
    result["mecB_present"] = "mecB" in result.get("genes_detected", [])
    result["confidence"] = {
        "mode": "coverage_only" if is_reads_mode else "full",
        "type_level": {
            "score": round(type_score, 3),
            "tier": type_tier,
        },
        "mec_component": mec_component,
        "ccr_component": ccr_component,
        "per_gene": per_gene,
    }
    result["assembly"] = {
        "contigs": ["Reads"] if is_reads_mode else sorted(contigs),
        "is_split": is_split,
        "split_penalty_applied": is_split,
        "gene_locations": gene_locations,
    }

    return result
```

- [ ] **Step 4: Run tests**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/test_confidence.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add lib/confidence.py tests/test_confidence.py
git commit -m "feat: add three-level confidence scoring module with tests"
```

---

## Chunk 2: Classifier Schema Updates

### Task 3: Update classifier.py for consistent JSON schema

**Files:**
- Modify: `lib/classifier.py`
- Modify: `bin/sccmec_typer.py`

The classifier's Negative early-return currently omits keys like `sccmec_type`, `mec_complex`, `confidence`, etc. The `enrich_result_with_confidence()` function handles most defaults, but the classifier should ensure the base result dict is consistent. Also wire in the confidence enrichment call.

- [ ] **Step 1: Update the Negative early-return in classifier.py**

In `classify_sccmec()`, change the empty-list early return from:
```python
if not hits_list:
    return {"status": "Negative", "reason": "No alignments found"}
```
to:
```python
if not hits_list:
    return {
        "status": "Negative",
        "reason": "No alignments found",
        "sccmec_type": "Negative",
        "mec_complex": "Negative",
        "ccr_complex": "Negative",
        "genes_detected": [],
        "mecA_present": False,
        "warnings": [],
        "hits_summary": [],
    }
```

Also update the second Negative return (no mec or ccr genes found) to include the same full set of keys.

- [ ] **Step 2: Wire confidence enrichment into sccmec_typer.py**

In `bin/sccmec_typer.py`, after the classification call, add the confidence enrichment:

```python
from lib.confidence import enrich_result_with_confidence

# After: result = classify_sccmec(hits)
result = enrich_result_with_confidence(result)
```

- [ ] **Step 3: Run existing classifier tests + confidence tests**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/ -v`
Expected: All PASS

- [ ] **Step 4: Commit**

```bash
git add lib/classifier.py bin/sccmec_typer.py
git commit -m "feat: consistent JSON schema + confidence enrichment in pipeline"
```

---

## Chunk 3: Jinja2 Templates

### Task 4: Create SVG and HTML Jinja2 templates

**Files:**
- Create: `templates/cassette.svg.j2`
- Create: `templates/coordinate_map.svg.j2`
- Create: `templates/report.html.j2`

These templates receive pre-computed data from the visualizer module. The SVG coordinate math is done in Python; the templates handle layout and styling.

- [ ] **Step 1: Create templates directory**

Run: `mkdir -p /home/alarawms/dev/sccmec_typer/templates`

- [ ] **Step 2: Create templates/cassette.svg.j2**

This template renders the schematic cassette diagram. It receives a `cassette` dict from the visualizer with pre-computed x/y positions, colors, and labels for each gene arrow.

```xml
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{{ width }}" height="{{ height }}" viewBox="0 0 {{ width }} {{ height }}">
  <defs>
    <marker id="arrowhead" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">
      <polygon points="0 0, 8 3, 0 6" fill="currentColor"/>
    </marker>
    <marker id="arrowhead-rev" markerWidth="8" markerHeight="6" refX="0" refY="3" orient="auto">
      <polygon points="8 0, 0 3, 8 6" fill="currentColor"/>
    </marker>
  </defs>

  <style>
    .title { font: bold 16px sans-serif; fill: #333; }
    .subtitle { font: 13px sans-serif; fill: #666; }
    .gene-label { font: 11px monospace; fill: #333; }
    .region-label { font: 10px sans-serif; fill: #999; text-anchor: middle; }
    .badge { font: bold 10px sans-serif; }
    .confidence-note { font: italic 10px sans-serif; fill: #999; }
    .warning-text { font: 11px sans-serif; fill: #E74C3C; }
    .missing { stroke-dasharray: 5,3; fill: none; stroke-width: 2; }
  </style>

  <!-- Title banner -->
  <text x="20" y="30" class="title">
    SCCmec {{ sccmec_type }}{% if designation %} ({{ designation }}){% endif %}
  </text>
  <text x="20" y="50" class="subtitle">
    mec: {{ mec_complex }} | ccr: {{ ccr_complex }} | Confidence: {{ type_confidence.tier }} ({{ "%.2f"|format(type_confidence.score) }})
  </text>

  <!-- Confidence badge -->
  <rect x="{{ width - 120 }}" y="15" width="100" height="28" rx="4"
        fill="{{ type_confidence.color }}" opacity="0.9"/>
  <text x="{{ width - 70 }}" y="34" text-anchor="middle" class="badge" fill="white">
    {{ type_confidence.tier }} {{ "%.0f"|format(type_confidence.score * 100) }}%
  </text>

  <!-- J region backgrounds -->
  {% for region in j_regions %}
  <rect x="{{ region.x }}" y="{{ region.y }}" width="{{ region.w }}" height="{{ region.h }}"
        fill="#F0F0F0" rx="2"/>
  <text x="{{ region.x + region.w / 2 }}" y="{{ region.y + region.h + 12 }}" class="region-label">
    {{ region.label }}
  </text>
  {% endfor %}

  <!-- DR boundaries -->
  {% for dr in dr_markers %}
  <polygon points="{{ dr.points }}" fill="#DC143C" opacity="0.8"/>
  {% endfor %}

  <!-- Gene arrows -->
  {% for gene in genes %}
  {% if gene.detected %}
  <rect x="{{ gene.x }}" y="{{ gene.y }}" width="{{ gene.w }}" height="{{ gene.h }}"
        fill="{{ gene.color }}" opacity="{{ gene.opacity }}" rx="2"/>
  {% if gene.strand == "-" %}
  <polygon points="{{ gene.arrow_points }}" fill="{{ gene.color }}" opacity="{{ gene.opacity }}"/>
  {% else %}
  <polygon points="{{ gene.arrow_points }}" fill="{{ gene.color }}" opacity="{{ gene.opacity }}"/>
  {% endif %}
  {% else %}
  <!-- Missing gene: dashed outline -->
  <rect x="{{ gene.x }}" y="{{ gene.y }}" width="{{ gene.w }}" height="{{ gene.h }}"
        class="missing" stroke="{{ gene.color }}" rx="2"/>
  {% endif %}
  <text x="{{ gene.x + gene.w / 2 }}" y="{{ gene.y - 4 }}" class="gene-label" text-anchor="middle">
    {{ gene.name }}{% if gene.contig_badge %} [{{ gene.contig_badge }}]{% endif %}
  </text>
  {% if gene.detected %}
  <text x="{{ gene.x + gene.w / 2 }}" y="{{ gene.y + gene.h + 12 }}" class="gene-label" text-anchor="middle"
        fill="{{ gene.tier_color }}">
    {{ "%.0f"|format(gene.score * 100) }}%
  </text>
  {% endif %}
  {% endfor %}

  <!-- Component confidence badges -->
  {% for badge in component_badges %}
  <text x="{{ badge.x }}" y="{{ badge.y }}" class="badge" fill="{{ badge.color }}">
    {{ badge.label }}: {{ badge.tier }} ({{ "%.2f"|format(badge.score) }})
  </text>
  {% endfor %}

  <!-- Warnings -->
  {% for i, warning in enumerate(warnings) %}
  <text x="20" y="{{ height - 20 + i * 15 }}" class="warning-text">⚠ {{ warning }}</text>
  {% endfor %}

  {% if confidence_note %}
  <text x="20" y="{{ height - 5 }}" class="confidence-note">{{ confidence_note }}</text>
  {% endif %}
</svg>
```

- [ ] **Step 3: Create templates/coordinate_map.svg.j2**

Renders genomic coordinate tracks, one per contig. Receives `tracks` list from visualizer.

```xml
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{{ width }}" height="{{ height }}" viewBox="0 0 {{ width }} {{ height }}">
  <style>
    .track-label { font: bold 12px monospace; fill: #333; }
    .bp-label { font: 9px monospace; fill: #999; }
    .gene-label { font: 10px monospace; fill: #333; text-anchor: middle; }
    .split-warning { font: 11px sans-serif; fill: #E74C3C; }
  </style>

  {% for track in tracks %}
  <!-- Contig track: {{ track.contig }} -->
  <text x="10" y="{{ track.y - 5 }}" class="track-label">
    {{ track.contig }} ({{ "{:,}".format(track.length) }} bp)
  </text>
  <line x1="{{ track.x }}" y1="{{ track.y }}" x2="{{ track.x + track.w }}" y2="{{ track.y }}"
        stroke="#CCC" stroke-width="2"/>

  {% for gene in track.genes %}
  <!-- Gene: {{ gene.name }} at {{ gene.start }}-{{ gene.end }} -->
  {% if gene.strand == "+" %}
  <polygon points="{{ gene.x }},{{ track.y - 10 }} {{ gene.x + gene.w - 6 }},{{ track.y - 10 }} {{ gene.x + gene.w }},{{ track.y }} {{ gene.x + gene.w - 6 }},{{ track.y + 10 }} {{ gene.x }},{{ track.y + 10 }}"
           fill="{{ gene.color }}" opacity="{{ gene.opacity }}"/>
  {% else %}
  <polygon points="{{ gene.x + 6 }},{{ track.y - 10 }} {{ gene.x + gene.w }},{{ track.y - 10 }} {{ gene.x + gene.w }},{{ track.y + 10 }} {{ gene.x + 6 }},{{ track.y + 10 }} {{ gene.x }},{{ track.y }}"
           fill="{{ gene.color }}" opacity="{{ gene.opacity }}"/>
  {% endif %}
  <text x="{{ gene.x + gene.w / 2 }}" y="{{ track.y - 15 }}" class="gene-label">{{ gene.name }}</text>
  <text x="{{ gene.x + gene.w / 2 }}" y="{{ track.y + 25 }}" class="bp-label" text-anchor="middle">
    {{ "{:,}".format(gene.start) }}
  </text>
  {% endfor %}
  {% endfor %}

  {% if is_split %}
  <line x1="10" y1="{{ split_line_y }}" x2="{{ width - 10 }}" y2="{{ split_line_y }}"
        stroke="#E74C3C" stroke-dasharray="8,4" stroke-width="1.5"/>
  <text x="{{ width / 2 }}" y="{{ split_line_y + 18 }}" class="split-warning" text-anchor="middle">
    ⚠ Split Assembly: components on separate contigs
  </text>
  {% endif %}
</svg>
```

- [ ] **Step 4: Create templates/report.html.j2**

Self-contained HTML report embedding both SVGs.

```html
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>SCCmec Report: {{ sample_name }}</title>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
         background: #f8f9fa; color: #333; padding: 20px; max-width: 1000px; margin: 0 auto; }
  .header { background: white; border-radius: 8px; padding: 24px; margin-bottom: 16px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
  .header h1 { font-size: 20px; margin-bottom: 8px; }
  .type-badge { display: inline-block; padding: 6px 16px; border-radius: 20px;
                font-weight: bold; font-size: 14px; color: white; }
  .type-badge.high { background: #27AE60; }
  .type-badge.medium { background: #F39C12; }
  .type-badge.low { background: #E74C3C; }
  .type-badge.na { background: #95A5A6; }
  .summary { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin-top: 12px; }
  .summary dt { font-weight: 600; color: #666; font-size: 12px; text-transform: uppercase; }
  .summary dd { font-size: 14px; margin-bottom: 8px; }
  .card { background: white; border-radius: 8px; padding: 20px; margin-bottom: 16px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
  .card h2 { font-size: 16px; margin-bottom: 12px; color: #555; }
  .card svg { width: 100%; height: auto; }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  th { background: #f1f3f5; padding: 8px 12px; text-align: left; font-weight: 600;
       border-bottom: 2px solid #dee2e6; }
  td { padding: 8px 12px; border-bottom: 1px solid #eee; }
  tr:hover td { background: #f8f9fa; }
  .conf-high { color: #27AE60; font-weight: bold; }
  .conf-medium { color: #F39C12; font-weight: bold; }
  .conf-low { color: #E74C3C; font-weight: bold; }
  .warning { background: #FFF3CD; border: 1px solid #FFEEBA; border-radius: 6px;
             padding: 12px 16px; margin-bottom: 16px; font-size: 13px; }
  .warning::before { content: "⚠ "; }
  .footer { text-align: center; color: #999; font-size: 11px; margin-top: 20px; }
  .confidence-note { font-style: italic; color: #999; font-size: 12px; margin-top: 8px; }
</style>
</head>
<body>
  <div class="header">
    <h1>SCCmec Typing Report: {{ sample_name }}</h1>
    {% if status == "Negative" %}
    <span class="type-badge na">No SCCmec Detected</span>
    {% else %}
    <span class="type-badge {{ confidence_tier_class }}">
      {{ sccmec_type }}{% if designation %} ({{ designation }}){% endif %}
      — {{ type_confidence.tier }} {{ "%.0f"|format(type_confidence.score * 100) }}%
    </span>
    {% endif %}
    <dl class="summary">
      <dt>Status</dt><dd>{{ status }}</dd>
      <dt>mec Complex</dt><dd>{{ mec_complex }}</dd>
      <dt>ccr Complex</dt><dd>{{ ccr_complex }}</dd>
      <dt>mecA</dt><dd>{{ "Present" if mecA_present else "Absent" }}</dd>
    </dl>
    {% if confidence_mode == "coverage_only" %}
    <p class="confidence-note">Confidence based on coverage only (reads mode — identity data unavailable).</p>
    {% endif %}
  </div>

  {% for warning in warnings %}
  <div class="warning">{{ warning }}</div>
  {% endfor %}

  <div class="card">
    <h2>Schematic Cassette Diagram</h2>
    {{ cassette_svg }}
  </div>

  {% if coordinate_map_svg %}
  <div class="card">
    <h2>Genomic Coordinate Map</h2>
    {{ coordinate_map_svg }}
  </div>
  {% endif %}

  {% if gene_table %}
  <div class="card">
    <h2>Gene Details</h2>
    <table>
      <thead>
        <tr>
          <th>Gene</th><th>Identity %</th><th>Coverage %</th>
          <th>Strand</th><th>Contig</th><th>Position</th><th>Confidence</th>
        </tr>
      </thead>
      <tbody>
        {% for gene in gene_table %}
        <tr>
          <td><strong>{{ gene.name }}</strong></td>
          <td>{{ "%.1f"|format(gene.id_pct) if gene.id_pct else "n/a" }}</td>
          <td>{{ "%.1f"|format(gene.cov_pct) }}</td>
          <td>{{ gene.strand }}</td>
          <td>{{ gene.contig }}</td>
          <td>{{ "{:,}".format(gene.start) }}–{{ "{:,}".format(gene.end) }}</td>
          <td class="conf-{{ gene.tier|lower }}">{{ gene.tier }} ({{ "%.0f"|format(gene.score * 100) }}%)</td>
        </tr>
        {% endfor %}
      </tbody>
    </table>
  </div>
  {% endif %}

  <div class="footer">
    Generated by SCCmec Typer | IWG-SCC Classification Types I–XV
  </div>
</body>
</html>
```

- [ ] **Step 5: Commit templates**

```bash
git add templates/
git commit -m "feat: add Jinja2 templates for cassette SVG, coordinate map, and HTML report"
```

---

## Chunk 4: Visualizer Module

### Task 5: Create visualizer module with tests

**Files:**
- Create: `lib/visualizer.py`
- Create: `tests/test_visualizer.py`

The visualizer computes SVG coordinate positions from the enriched result dict, then renders Jinja2 templates.

- [ ] **Step 1: Write tests/test_visualizer.py**

```python
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.visualizer import generate_cassette_svg, generate_report_html


class TestCassetteSVG:
    """Verify SVG generation for schematic diagram."""

    def test_positive_result_produces_svg(self, tmp_path):
        result = _make_positive_result()
        svg = generate_cassette_svg(result)
        assert svg.startswith("<?xml")
        assert "SCCmec Type II" in svg
        assert "mecA" in svg
        assert "ccrA2" in svg

    def test_negative_result_produces_svg(self, tmp_path):
        result = _make_negative_result()
        svg = generate_cassette_svg(result)
        assert "No SCCmec" in svg

    def test_split_assembly_shows_badges(self):
        result = _make_split_result()
        svg = generate_cassette_svg(result)
        assert "[c1]" in svg or "[c2]" in svg


class TestHTMLReport:
    """Verify HTML report generation."""

    def test_produces_html(self, tmp_path):
        result = _make_positive_result()
        html = generate_report_html(result, sample_name="test_sample")
        assert "<!DOCTYPE html>" in html
        assert "test_sample" in html
        assert "Gene Details" in html

    def test_negative_result_report(self):
        result = _make_negative_result()
        html = generate_report_html(result, sample_name="neg_sample")
        assert "No SCCmec Detected" in html


def _make_positive_result():
    """Helper: enriched positive Type II result."""
    return {
        "status": "Positive",
        "sccmec_type": "Type II",
        "designation": "2A",
        "mec_complex": "Class A",
        "ccr_complex": "Type 2",
        "ccr_genes": ["ccrA2", "ccrB2"],
        "genes_detected": ["IS431", "ccrA2", "ccrB2", "mecA", "mecI", "mecR1"],
        "mecA_present": True,
        "mecC_present": False,
        "mecB_present": False,
        "confidence": {
            "mode": "full",
            "type_level": {"score": 0.92, "tier": "High"},
            "mec_component": {
                "class": "Class A", "score": 0.99, "tier": "High",
                "genes_expected": ["mecA", "mecR1", "mecI"],
                "genes_detected": ["mecA", "mecR1", "mecI"],
                "genes_missing": [],
                "orientation": None, "orientation_resolved": None,
            },
            "ccr_component": {
                "type": "Type 2", "score": 0.93, "tier": "High",
                "genes_expected": ["ccrA2", "ccrB2"],
                "genes_detected": ["ccrA2", "ccrB2"],
                "genes_missing": [],
            },
            "per_gene": {
                "mecA": {"identity_pct": 100, "coverage_pct": 100, "score": 1.0, "tier": "High"},
                "mecR1": {"identity_pct": 99, "coverage_pct": 99, "score": 0.97, "tier": "High"},
                "mecI": {"identity_pct": 99, "coverage_pct": 100, "score": 0.99, "tier": "High"},
                "ccrA2": {"identity_pct": 99, "coverage_pct": 98, "score": 0.94, "tier": "High"},
                "ccrB2": {"identity_pct": 99, "coverage_pct": 97, "score": 0.92, "tier": "High"},
                "IS431": {"identity_pct": 98, "coverage_pct": 95, "score": 0.78, "tier": "High"},
            },
        },
        "assembly": {
            "contigs": ["contig_1"], "is_split": False, "split_penalty_applied": False,
            "gene_locations": [
                {"gene": "mecA", "contig": "contig_1", "start": 28100, "end": 30100, "strand": "+"},
                {"gene": "mecR1", "contig": "contig_1", "start": 30200, "end": 32000, "strand": "+"},
                {"gene": "mecI", "contig": "contig_1", "start": 32100, "end": 33000, "strand": "+"},
                {"gene": "IS431", "contig": "contig_1", "start": 27000, "end": 28000, "strand": "+"},
                {"gene": "ccrA2", "contig": "contig_1", "start": 12000, "end": 14000, "strand": "-"},
                {"gene": "ccrB2", "contig": "contig_1", "start": 10000, "end": 12000, "strand": "-"},
            ],
        },
        "warnings": [],
        "hits_summary": [],
    }


def _make_negative_result():
    return {
        "status": "Negative", "reason": "No alignments found",
        "sccmec_type": "Negative", "designation": None,
        "mec_complex": "Negative", "ccr_complex": "Negative",
        "ccr_genes": [], "genes_detected": [],
        "mecA_present": False, "mecC_present": False, "mecB_present": False,
        "confidence": {"mode": "n/a", "type_level": {"score": 0.0, "tier": "n/a"},
                       "mec_component": None, "ccr_component": None, "per_gene": {}},
        "assembly": {"contigs": [], "is_split": False, "split_penalty_applied": False,
                     "gene_locations": []},
        "warnings": [], "hits_summary": [],
    }


def _make_split_result():
    result = _make_positive_result()
    result["assembly"]["is_split"] = True
    result["assembly"]["contigs"] = ["contig_1", "contig_2"]
    result["assembly"]["gene_locations"][4]["contig"] = "contig_2"
    result["assembly"]["gene_locations"][5]["contig"] = "contig_2"
    result["warnings"] = ["Split Assembly: Components found on 2 different contigs"]
    return result
```

- [ ] **Step 2: Run tests to verify they FAIL**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/test_visualizer.py -v 2>&1 | head -20`
Expected: ImportError

- [ ] **Step 3: Create lib/visualizer.py**

This is a learning opportunity — the gene layout algorithm for the schematic cassette diagram involves positioning genes in the canonical SCCmec order (DR-J3-ccr-J2-mec-J1-DR) regardless of actual genomic position. I'll create the scaffold and the coordinate map logic; you can implement the `_layout_cassette_genes()` function which decides the x-position and width of each gene arrow in the schematic.

```python
import os
from jinja2 import Environment, FileSystemLoader

# IWG-SCC color scheme
GENE_COLORS = {
    "mecA": "#DC143C", "mecC": "#DC143C", "mecB": "#DC143C",
    "mecR1": "#E8A0BF", "mecI": "#E8A0BF",
    "blaZ": "#8E44AD",
}
CCR_COLOR = "#7B68EE"
IS_COLOR = "#95A5A6"
TIER_COLORS = {"High": "#27AE60", "Medium": "#F39C12", "Low": "#E74C3C", "n/a": "#95A5A6"}

# Canonical gene order in SCCmec cassette
# Left to right: DR - J3 - [ccr genes] - J2 - [mec genes + IS] - J1 - DR
MEC_GENES = {"mecA", "mecC", "mecB", "mecR1", "mecI", "IS431", "IS431_1", "IS431_2", "IS1272", "blaZ"}
CCR_GENES_PREFIX = "ccr"


def _get_gene_color(gene_name):
    if gene_name in GENE_COLORS:
        return GENE_COLORS[gene_name]
    if gene_name.startswith("ccr"):
        return CCR_COLOR
    if gene_name.startswith("IS"):
        return IS_COLOR
    return "#95A5A6"


def _get_template_env():
    templates_dir = os.path.join(os.path.dirname(__file__), "../templates")
    return Environment(loader=FileSystemLoader(templates_dir), autoescape=False)


def _layout_cassette_genes(result):
    """Compute x/y positions for genes in the schematic cassette diagram.

    Places genes in canonical SCCmec order:
    DR(left) - J3 - [ccr genes] - J2 - [mec complex genes] - J1 - DR(right)

    Returns a list of gene dicts with layout fields added:
    {name, detected, color, opacity, x, y, w, h, strand, score, tier, tier_color,
     contig_badge, arrow_points}
    """
    genes_detected = set(result.get("genes_detected", []))
    per_gene = result.get("confidence", {}).get("per_gene", {})
    is_split = result.get("assembly", {}).get("is_split", False)

    # Build location lookup for contig badges
    gene_contigs = {}
    for loc in result.get("assembly", {}).get("gene_locations", []):
        gene_contigs[loc["gene"]] = loc["contig"]

    # Collect unique contigs for badge abbreviation
    unique_contigs = sorted(set(gene_contigs.values()))
    contig_abbrev = {c: f"c{i+1}" for i, c in enumerate(unique_contigs)}

    # Separate ccr and mec genes
    ccr_genes_list = sorted([g for g in genes_detected if g.startswith("ccr")])
    mec_genes_list = sorted([g for g in genes_detected if g in MEC_GENES])

    # Also include expected but missing genes
    mec_comp = result.get("confidence", {}).get("mec_component")
    ccr_comp = result.get("confidence", {}).get("ccr_component")

    if mec_comp and mec_comp.get("genes_missing"):
        for g in mec_comp["genes_missing"]:
            if g not in mec_genes_list:
                mec_genes_list.append(g)

    if ccr_comp and ccr_comp.get("genes_missing"):
        for g in ccr_comp["genes_missing"]:
            if g not in ccr_genes_list:
                ccr_genes_list.append(g)

    # Layout constants
    svg_width = 900
    gene_h = 24
    gene_y = 90
    dr_w = 30
    j_w = 40
    gene_spacing = 8
    margin = 20

    # Calculate available width for gene blocks
    structural_w = 2 * dr_w + 3 * j_w + 6 * gene_spacing
    available_w = svg_width - 2 * margin - structural_w
    total_genes = len(ccr_genes_list) + len(mec_genes_list)
    gene_w = min(80, available_w / max(total_genes, 1))

    # Position elements left to right
    x = margin
    layout = []

    # DR left
    dr_left = {"points": f"{x},{gene_y - 12} {x + dr_w},{gene_y} {x},{gene_y + 12}"}
    x += dr_w + gene_spacing

    # J3 region
    j3 = {"x": x, "y": gene_y - 15, "w": j_w, "h": 30, "label": "J3"}
    x += j_w + gene_spacing

    # ccr genes
    for gene_name in ccr_genes_list:
        detected = gene_name in genes_detected
        gene_info = per_gene.get(gene_name, {})
        score = gene_info.get("score", 0)
        tier = gene_info.get("tier", "Low")
        opacity = max(0.4, score) if detected else 1.0

        contig_badge = None
        if is_split and gene_name in gene_contigs:
            contig_badge = contig_abbrev.get(gene_contigs[gene_name])

        gx = x
        # Arrow points: right-pointing arrow
        arrow_points = f"{gx},{gene_y - 12} {gx + gene_w - 6},{gene_y - 12} {gx + gene_w},{gene_y} {gx + gene_w - 6},{gene_y + 12} {gx},{gene_y + 12}"

        layout.append({
            "name": gene_name, "detected": detected,
            "color": _get_gene_color(gene_name), "opacity": opacity,
            "x": gx, "y": gene_y - 12, "w": gene_w, "h": gene_h,
            "strand": "+", "score": score, "tier": tier,
            "tier_color": TIER_COLORS.get(tier, "#95A5A6"),
            "contig_badge": contig_badge, "arrow_points": arrow_points,
        })
        x += gene_w + gene_spacing

    # J2 region
    j2 = {"x": x, "y": gene_y - 15, "w": j_w, "h": 30, "label": "J2"}
    x += j_w + gene_spacing

    # mec genes
    for gene_name in mec_genes_list:
        detected = gene_name in genes_detected
        gene_info = per_gene.get(gene_name, {})
        score = gene_info.get("score", 0)
        tier = gene_info.get("tier", "Low")
        opacity = max(0.4, score) if detected else 1.0

        contig_badge = None
        if is_split and gene_name in gene_contigs:
            contig_badge = contig_abbrev.get(gene_contigs[gene_name])

        gx = x
        arrow_points = f"{gx},{gene_y - 12} {gx + gene_w - 6},{gene_y - 12} {gx + gene_w},{gene_y} {gx + gene_w - 6},{gene_y + 12} {gx},{gene_y + 12}"

        layout.append({
            "name": gene_name, "detected": detected,
            "color": _get_gene_color(gene_name), "opacity": opacity,
            "x": gx, "y": gene_y - 12, "w": gene_w, "h": gene_h,
            "strand": "+", "score": score, "tier": tier,
            "tier_color": TIER_COLORS.get(tier, "#95A5A6"),
            "contig_badge": contig_badge, "arrow_points": arrow_points,
        })
        x += gene_w + gene_spacing

    # J1 region
    j1 = {"x": x, "y": gene_y - 15, "w": j_w, "h": 30, "label": "J1"}
    x += j_w + gene_spacing

    # DR right
    dr_right = {"points": f"{x},{gene_y - 12} {x + dr_w},{gene_y} {x},{gene_y + 12}"}

    j_regions = [j3, j2, j1]
    dr_markers = [dr_left, dr_right]

    return layout, j_regions, dr_markers, svg_width


def _layout_coordinate_tracks(result):
    """Compute coordinate map tracks from gene_locations.

    Groups genes by contig, sorts by start position, computes
    proportional x-positions within each track.

    Returns (tracks, total_height, is_split).
    """
    locations = result.get("assembly", {}).get("gene_locations", [])
    if not locations:
        return [], 0, False

    from collections import defaultdict
    contig_genes = defaultdict(list)
    for loc in locations:
        contig_genes[loc["contig"]].append(loc)

    svg_width = 900
    margin = 60
    track_height = 60
    track_spacing = 40
    tracks = []
    y = 30

    for contig in sorted(contig_genes.keys()):
        genes = sorted(contig_genes[contig], key=lambda g: g["start"])
        if not genes:
            continue

        # Find the range of positions on this contig
        min_pos = min(g["start"] for g in genes)
        max_pos = max(g["end"] for g in genes)
        contig_span = max_pos - min_pos if max_pos > min_pos else 1

        track_w = svg_width - 2 * margin
        per_gene = result.get("confidence", {}).get("per_gene", {})

        track_genes = []
        for g in genes:
            gene_info = per_gene.get(g["gene"], {})
            score = gene_info.get("score", 0.5)
            rel_start = (g["start"] - min_pos) / contig_span
            rel_end = (g["end"] - min_pos) / contig_span
            gx = margin + rel_start * track_w
            gw = max(10, (rel_end - rel_start) * track_w)

            track_genes.append({
                "name": g["gene"],
                "start": g["start"],
                "end": g["end"],
                "strand": g.get("strand", "+"),
                "color": _get_gene_color(g["gene"]),
                "opacity": max(0.4, score),
                "x": gx,
                "w": gw,
            })

        tracks.append({
            "contig": contig,
            "length": max_pos,
            "x": margin,
            "y": y + 20,
            "w": track_w,
            "genes": track_genes,
        })
        y += track_height + track_spacing

    is_split = len(tracks) > 1
    total_height = y + (30 if is_split else 10)
    return tracks, total_height, is_split


def generate_cassette_svg(result):
    """Generate the schematic cassette diagram SVG string."""
    env = _get_template_env()

    if result.get("status") == "Negative":
        # Minimal negative SVG
        return env.get_template("cassette.svg.j2").render(
            width=900, height=120,
            sccmec_type="Negative", designation=None,
            mec_complex="Negative", ccr_complex="Negative",
            type_confidence={"score": 0, "tier": "n/a", "color": "#95A5A6"},
            genes=[], j_regions=[], dr_markers=[], component_badges=[],
            warnings=["No SCCmec detected"],
            confidence_note=None,
        )

    genes, j_regions, dr_markers, svg_width = _layout_cassette_genes(result)

    type_conf = result.get("confidence", {}).get("type_level", {"score": 0, "tier": "n/a"})
    type_conf["color"] = TIER_COLORS.get(type_conf["tier"], "#95A5A6")

    # Component badges
    badges = []
    mec_comp = result.get("confidence", {}).get("mec_component")
    ccr_comp = result.get("confidence", {}).get("ccr_component")

    if mec_comp:
        badges.append({
            "x": svg_width / 2 + 100, "y": 140,
            "label": f"mec ({mec_comp['class']})",
            "score": mec_comp["score"], "tier": mec_comp["tier"],
            "color": TIER_COLORS.get(mec_comp["tier"], "#95A5A6"),
        })
    if ccr_comp and ccr_comp.get("type") != "n/a":
        badges.append({
            "x": 100, "y": 140,
            "label": f"ccr ({ccr_comp['type']})",
            "score": ccr_comp["score"], "tier": ccr_comp["tier"],
            "color": TIER_COLORS.get(ccr_comp["tier"], "#95A5A6"),
        })

    conf_mode = result.get("confidence", {}).get("mode", "full")
    conf_note = "Confidence based on coverage only (reads mode)" if conf_mode == "coverage_only" else None

    height = 170 + len(result.get("warnings", [])) * 15
    return env.get_template("cassette.svg.j2").render(
        width=svg_width, height=height,
        sccmec_type=result.get("sccmec_type", "Unknown"),
        designation=result.get("designation"),
        mec_complex=result.get("mec_complex", "Unknown"),
        ccr_complex=result.get("ccr_complex", "Unknown"),
        type_confidence=type_conf,
        genes=genes, j_regions=j_regions, dr_markers=dr_markers,
        component_badges=badges,
        warnings=result.get("warnings", []),
        confidence_note=conf_note,
        enumerate=enumerate,
    )


def generate_coordinate_map_svg(result):
    """Generate the genomic coordinate map SVG string. Returns None for reads mode."""
    if result.get("assembly", {}).get("contigs") == ["Reads"]:
        return None
    if not result.get("assembly", {}).get("gene_locations"):
        return None

    env = _get_template_env()
    tracks, total_height, is_split = _layout_coordinate_tracks(result)

    if not tracks:
        return None

    split_line_y = total_height - 40 if is_split else 0

    return env.get_template("coordinate_map.svg.j2").render(
        width=900, height=total_height,
        tracks=tracks, is_split=is_split, split_line_y=split_line_y,
    )


def generate_report_html(result, sample_name="sample"):
    """Generate the full HTML report with embedded SVGs and gene table."""
    env = _get_template_env()

    cassette_svg = generate_cassette_svg(result)
    coord_svg = generate_coordinate_map_svg(result)

    # Build gene table from per_gene confidence data
    gene_table = []
    per_gene = result.get("confidence", {}).get("per_gene", {})
    locations = {loc["gene"]: loc for loc in result.get("assembly", {}).get("gene_locations", [])}

    for gene_name in sorted(per_gene.keys()):
        info = per_gene[gene_name]
        loc = locations.get(gene_name, {})
        gene_table.append({
            "name": gene_name,
            "id_pct": info.get("identity_pct", 0),
            "cov_pct": info.get("coverage_pct", 0),
            "strand": loc.get("strand", "n/a"),
            "contig": loc.get("contig", "n/a"),
            "start": loc.get("start", 0),
            "end": loc.get("end", 0),
            "score": info.get("score", 0),
            "tier": info.get("tier", "n/a"),
        })

    type_conf = result.get("confidence", {}).get("type_level", {"score": 0, "tier": "n/a"})
    tier_class = type_conf["tier"].lower() if type_conf["tier"] != "n/a" else "na"

    return env.get_template("report.html.j2").render(
        sample_name=sample_name,
        status=result.get("status", "Unknown"),
        sccmec_type=result.get("sccmec_type", "Unknown"),
        designation=result.get("designation"),
        mec_complex=result.get("mec_complex", "Unknown"),
        ccr_complex=result.get("ccr_complex", "Unknown"),
        mecA_present=result.get("mecA_present", False),
        type_confidence=type_conf,
        confidence_tier_class=tier_class,
        confidence_mode=result.get("confidence", {}).get("mode", "full"),
        cassette_svg=cassette_svg,
        coordinate_map_svg=coord_svg,
        gene_table=gene_table if gene_table else None,
        warnings=result.get("warnings", []),
    )


def write_visualization(result, output_prefix, sample_name="sample"):
    """Write SVG and HTML report files.

    Produces:
      {output_prefix}_map.svg    — Schematic cassette diagram
      {output_prefix}_report.html — Full HTML report
    """
    svg = generate_cassette_svg(result)
    svg_path = f"{output_prefix}_map.svg"
    with open(svg_path, "w") as f:
        f.write(svg)

    html = generate_report_html(result, sample_name=sample_name)
    html_path = f"{output_prefix}_report.html"
    with open(html_path, "w") as f:
        f.write(html)

    return svg_path, html_path
```

- [ ] **Step 4: Run tests**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/test_visualizer.py -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add lib/visualizer.py tests/test_visualizer.py
git commit -m "feat: add SVG/HTML visualization module with Jinja2 templates"
```

---

## Chunk 5: Pipeline Integration

### Task 6: Wire visualization into main entry point

**Files:**
- Modify: `bin/sccmec_typer.py`

- [ ] **Step 1: Add --no-viz flag and visualization calls**

After the existing output section in `bin/sccmec_typer.py`, add:

```python
# Add argument to parser
parser.add_argument("--no-viz", action="store_true",
                    help="Skip SVG/HTML visualization generation")

# After JSON/TSV/CSV output section, add:
if not args.no_viz:
    from lib.visualizer import write_visualization
    sample_name = os.path.basename(args.input1)
    svg_path, html_path = write_visualization(result, args.output, sample_name=sample_name)
    print(f"SVG cassette diagram written to {svg_path}")
    print(f"HTML report written to {html_path}")
```

- [ ] **Step 2: Verify full pipeline with a reference strain**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 bin/sccmec_typer.py --1 tests/data/COL.fna -d db/sccmec_targets.fasta -o /tmp/COL_viz_test`
Expected: Produces `/tmp/COL_viz_test_map.svg`, `/tmp/COL_viz_test_report.html`, and JSON with confidence block.

- [ ] **Step 3: Verify --no-viz flag**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 bin/sccmec_typer.py --1 tests/data/COL.fna -d db/sccmec_targets.fasta -o /tmp/COL_noviz_test --no-viz`
Expected: No `_map.svg` or `_report.html` files produced. JSON still has confidence block.

- [ ] **Step 4: Run all tests**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/ -v`
Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add bin/sccmec_typer.py
git commit -m "feat: wire SVG visualization into pipeline with --no-viz flag"
```

---

## Summary of Changes

| Component | Before | After |
|-----------|--------|-------|
| Confidence scoring | None | **Three-level (gene/component/type)** |
| SVG output | None | **Schematic cassette diagram** |
| HTML output | None | **Interactive report with coordinate map** |
| JSON schema | Basic | **Enhanced with confidence + assembly + designation** |
| Reads mode | No confidence | **Coverage-only confidence with mode flag** |
| Negative results | Incomplete schema | **Consistent schema across all statuses** |
| Staphit integration | TSV only | **Full JSON with confidence for pipeline consumption** |
| Dependencies | minimap2 + pandas | **+ jinja2** |
