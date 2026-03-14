import json
import os


def _load_rules():
    """Load classification rules for expected gene lookups."""
    rules_path = os.path.join(os.path.dirname(__file__), "../db/rules.json")
    try:
        with open(rules_path, "r") as f:
            return json.load(f)
    except Exception as e:
        raise RuntimeError(f"confidence.py: Failed to load rules.json: {e}") from e


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

        orientation_undetermined = mec_complex.startswith("Class C")
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
