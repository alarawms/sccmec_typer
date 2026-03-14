# SCCmec Typer: Confidence Scoring Methods

## Overview

The confidence scoring system provides a three-tiered assessment of classification reliability, from individual gene alignments up to the final SCCmec type assignment. This enables users to gauge the trustworthiness of each result and identify which components may require further investigation.

## Level 1: Per-Gene Confidence

Each gene detected by alignment is assigned a confidence score ranging from 0.0 to 1.0, derived from alignment identity and reference coverage:

```
norm_id  = clamp((identity_percent - 90) / 10,  0.0, 1.0)
norm_cov = clamp((coverage_percent - 80) / 20, 0.0, 1.0)
gene_confidence = (W_id * norm_id) + (W_cov * norm_cov)
```

**Parameters:**
- `norm_id = clamp((identity_percent - 90) / 10, 0.0, 1.0)`
  - Maps the identity range 90%–100% to 0.0–1.0
  - 90% is the minimum identity threshold for a valid hit
  - Clamped to [0.0, 1.0] to handle edge cases where PAF identity exceeds 100%
- `norm_cov = clamp((coverage_percent - 80) / 20, 0.0, 1.0)`
  - Maps the coverage range 80%–100% to 0.0–1.0
  - 80% is the minimum coverage threshold for a valid hit
  - Clamped to [0.0, 1.0] to handle split alignments where coverage may exceed 100%
- `W_id = 0.6` (identity weight)
- `W_cov = 0.4` (coverage weight)

**Rationale for weighting:** Sequence identity is weighted more heavily than coverage because a high-identity partial alignment more reliably identifies the correct gene than a full-length but divergent alignment. Divergent identity may indicate a novel allele or paralogous gene, while reduced coverage may simply reflect assembly fragmentation.

**Classification thresholds:**
| Tier | Score Range | Interpretation |
|------|-------------|---------------|
| High | >= 0.75 | Approximately >97.5% identity and >95% coverage |
| Medium | 0.40 – 0.75 | Moderate divergence or incomplete coverage |
| Low | < 0.40 | Borderline hits near filter cutoffs |

### Reads Mode Fallback

In reads mode, `coverage.py` aggregates read alignments into breadth-of-coverage per gene but does not compute per-gene identity (`id_pct = 0`). When `id_pct` is 0 or unavailable, confidence uses **coverage-only mode**:

```
gene_confidence = clamp((coverage_percent - 80) / 20, 0.0, 1.0)
```

The identity term is dropped entirely rather than contributing a zero. The JSON output includes `confidence.mode = "coverage_only"` to indicate reduced-information scoring. Tier thresholds remain the same.

## Level 2: Component Confidence

Each structural component (mec complex class and ccr complex type) receives a confidence score that accounts for both completeness and alignment quality:

```
component_confidence = (genes_detected / genes_expected) * mean(gene_confidences)
```

**Parameters:**
- `genes_expected`: The number of genes defined as required or characteristic for the matched rule in `rules.json`. For example, Class A mec complex expects mecA, mecR1, and mecI (3 genes); ccr Type 2 expects ccrA2 and ccrB2 (2 genes).
- `genes_detected`: The number of those expected genes actually found in the alignment results.
- `mean(gene_confidences)`: The arithmetic mean of per-gene confidence scores (Level 1) for the detected genes.

**Orientation uncertainty penalty:** When the mec complex classification involves IS431 orientation analysis (distinguishing Class C1 from Class C2) and the orientation cannot be resolved — either because only one IS431 copy was detected, or because strand information is unavailable (reads mode) — a multiplicative penalty of 0.8 is applied:

```
if orientation_undetermined and mec_class.startswith("Class C"):
    # Applies to Class C, Class C1, Class C2 only — NOT Class E (mecC-based)
    component_confidence *= 0.80
```

**Rationale:** The penalty reflects that the C1/C2 distinction is based on IS431 copy orientation, which requires either (a) two distinctly named IS431 copies with strand information in assembly mode, or (b) structural analysis not available from short/long read mapping alone. When this information is absent, the mec class assignment carries higher uncertainty.

## Level 3: Type-Level Confidence

The overall confidence in the SCCmec type assignment combines the two component scores with an assembly quality factor:

```
type_confidence = mec_confidence * ccr_confidence * assembly_penalty
```

**Parameters:**
- `mec_confidence`: Component confidence for the mec complex (Level 2)
- `ccr_confidence`: Component confidence for the ccr complex (Level 2). When the matched SCCmec type rule specifies `ccr: "ANY"` (e.g., Plasmid-borne mecB), no ccr complex is expected and `ccr_confidence = 1.0` (neutral factor).
- `assembly_penalty`:
  - `1.0` if all SCCmec components are detected on a single contig
  - `0.85` if components are split across multiple contigs (split assembly)

**Rationale for assembly penalty:** SCCmec is a single mobile genetic element. When its components are detected on different contigs, this may indicate assembly fragmentation, or it may indicate that the detected genes belong to separate genomic loci rather than a single SCCmec cassette. The 0.85 penalty reflects this ambiguity without outright rejecting the classification.

**Classification thresholds:**
| Tier | Score Range | Interpretation |
|------|-------------|---------------|
| High | > 0.70 | Strong evidence: all expected genes present with high-quality alignments on a single contig |
| Medium | 0.40 – 0.70 | Moderate evidence: some genes may be missing, borderline alignments, or split assembly |
| Low | < 0.40 | Weak evidence: significant missing components, poor alignments, or multiple confounding factors |

## Example Calculations

### Example 1: High-Confidence Type II (N315)

| Gene | Identity (%) | Coverage (%) | Gene Confidence |
|------|-------------|-------------|-----------------|
| mecA | 100.0 | 100.0 | 0.6 * 1.0 + 0.4 * 1.0 = **1.00** |
| mecR1 | 99.8 | 99.5 | 0.6 * 0.98 + 0.4 * 0.975 = **0.978** |
| mecI | 99.9 | 100.0 | 0.6 * 0.99 + 0.4 * 1.0 = **0.994** |
| ccrA2 | 99.7 | 98.0 | 0.6 * 0.97 + 0.4 * 0.90 = **0.942** |
| ccrB2 | 99.5 | 97.5 | 0.6 * 0.95 + 0.4 * 0.875 = **0.920** |

- mec_confidence = (3/3) * mean(1.00, 0.978, 0.994) = **0.991**
- ccr_confidence = (2/2) * mean(0.942, 0.920) = **0.931**
- assembly_penalty = 1.0 (single contig)
- **type_confidence = 0.991 * 0.931 * 1.0 = 0.923 (High)**

### Example 2: Low-Confidence Fragmented Assembly

| Gene | Identity (%) | Coverage (%) | Gene Confidence |
|------|-------------|-------------|-----------------|
| mecA | 92.0 | 85.0 | 0.6 * 0.20 + 0.4 * 0.25 = **0.220** |
| mecR1 | 91.5 | 82.0 | 0.6 * 0.15 + 0.4 * 0.10 = **0.130** |
| ccrA2 | 93.0 | 88.0 | 0.6 * 0.30 + 0.4 * 0.40 = **0.340** |
| ccrB2 | — | — | not detected |

- mec_confidence = (2/3) * mean(0.220, 0.130) = 0.667 * 0.175 = **0.117** (mecI missing)
- ccr_confidence = (1/2) * mean(0.340) = 0.500 * 0.340 = **0.170** (ccrB2 missing)
- assembly_penalty = 0.85 (split across 2 contigs)
- **type_confidence = 0.117 * 0.170 * 0.85 = 0.017 (Low)**
