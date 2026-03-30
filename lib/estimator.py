import json
import os


# Gene importance weights (IWG-SCC biological hierarchy)
GENE_WEIGHTS = {
    "mecA": 1.0, "mecC": 1.0, "mecB": 1.0,
    "ccrA1": 0.8, "ccrA2": 0.8, "ccrA3": 0.8, "ccrA4": 0.8,
    "ccrB1": 0.8, "ccrB2": 0.8, "ccrB3": 0.8, "ccrB4": 0.8, "ccrB6": 0.8,
    "ccrC1": 0.8, "ccrC2": 0.8,
    "mecR1": 0.5, "mecI": 0.5,
    "IS431": 0.3, "IS431_1": 0.3, "IS431_2": 0.3, "IS1272": 0.3,
    "blaZ": 0.2,
}
DEFAULT_WEIGHT = 0.2


def _load_rules():
    rules_path = os.path.join(os.path.dirname(__file__), "../db/rules.json")
    with open(rules_path, "r") as f:
        return json.load(f)


def _clamp(val, lo, hi):
    return max(lo, min(hi, val))


def _get_expected_genes(type_rule, rules):
    """Get the expected gene set for an SCCmec type rule."""
    mec_name = type_rule["mec"]
    ccr_name = type_rule["ccr"]
    genes = []

    # mec complex genes (required + any_of)
    for rule in rules["mec_complex_rules"]:
        if rule["name"] == mec_name:
            genes.extend(rule["required"])
            # Include any_of genes — at least one should be present
            if "any_of" in rule:
                genes.extend(rule["any_of"])
            break

    # ccr complex genes (skip if "ANY")
    if ccr_name != "ANY":
        for rule in rules["ccr_complex_rules"]:
            if rule["name"] == ccr_name:
                genes.extend(rule["required_genes"])
                break

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique.append(g)
    return unique


def _compute_gene_contribution(gene, hits_by_gene, soft_hits_by_gene):
    """Compute a gene's contribution to a candidate score.

    Returns (contribution, status, id_pct, cov_pct).
    """
    weight = GENE_WEIGHTS.get(gene, DEFAULT_WEIGHT)

    if gene in hits_by_gene:
        h = hits_by_gene[gene]
        return weight * 1.0, "present", h["id_pct"], h["cov_pct"]

    if gene in soft_hits_by_gene:
        h = soft_hits_by_gene[gene]
        # Normalized quality: how close to passing thresholds
        norm_id = _clamp((h["id_pct"] - 70) / 20, 0.0, 1.0)
        norm_cov = _clamp((h["cov_pct"] - 50) / 30, 0.0, 1.0)
        quality = (norm_id + norm_cov) / 2
        return weight * quality, "sub_threshold", h["id_pct"], h["cov_pct"]

    return 0.0, "absent", 0.0, 0.0


def _generate_ruling(type_rule, gene_evidence, result):
    """Generate a human-readable differential diagnosis ruling."""
    sub_threshold = [g for g in gene_evidence if g["status"] == "sub_threshold"]
    absent = [g for g in gene_evidence if g["status"] == "absent"]
    type_name = type_rule["name"]
    mec_name = type_rule["mec"]
    ccr_name = type_rule["ccr"]

    parts = []

    if sub_threshold:
        for g in sub_threshold:
            metric = "identity" if g["id_pct"] < 90 else "coverage"
            value = g["id_pct"] if metric == "identity" else g["cov_pct"]
            threshold = 90 if metric == "identity" else 80
            parts.append(
                f"{g['gene']} {metric} ({value:.0f}%) is below the "
                f"{threshold}% threshold"
            )

    if absent:
        missing_names = ", ".join(g["gene"] for g in absent)
        parts.append(f"{missing_names} not detected")

    if result.get("mec_complex") == "Negative":
        parts.append(
            f"no mec gene detected; if {mec_name} genes were present "
            f"with {ccr_name}, this would classify as {type_name}"
        )
    elif result.get("ccr_complex") == "Negative" and not sub_threshold:
        parts.append(
            f"no ccr genes detected; {type_name} requires {ccr_name}"
        )

    if parts:
        ruling = "; ".join(parts)
        if sub_threshold and not absent:
            ruling += (
                f"; if included, {mec_name} + {ccr_name} "
                f"would classify as {type_name}"
            )
        return ruling

    return f"detected components do not match {type_name} ({mec_name} + {ccr_name})"


def _determine_candidate_status(gene_evidence):
    """Determine the candidate match status from gene evidence."""
    statuses = set(g["status"] for g in gene_evidence)
    if "absent" not in statuses and "sub_threshold" not in statuses:
        return "full_match"
    if "absent" not in statuses and "sub_threshold" in statuses:
        return "full_match_sub_threshold"
    if "present" in statuses:
        return "partial_match"
    return "no_match"


def estimate_closest_types(result, hits, soft_hits, max_candidates=3):
    """Estimate the closest SCCmec type(s) when definitive classification fails.

    Args:
        result: classifier result dict
        hits: list of above-threshold hit dicts
        soft_hits: list of sub-threshold hit dicts
        max_candidates: maximum number of candidates to return

    Returns:
        estimation dict with candidates, or None if not applicable.
    """
    status = result.get("status", "")
    sccmec_type = result.get("sccmec_type", "")

    # Don't estimate for definitive Positive (with known type) or Negative
    if status == "Negative" or sccmec_type == "Negative":
        return None
    if status == "Positive" and sccmec_type != "Unknown" and "Composite" not in sccmec_type:
        return None

    rules = _load_rules()

    # Build gene lookup dicts (best hit per gene)
    hits_by_gene = {}
    for h in hits:
        gene = h["gene"]
        if gene not in hits_by_gene or h.get("cov_pct", 0) > hits_by_gene[gene].get("cov_pct", 0):
            hits_by_gene[gene] = h

    soft_hits_by_gene = {}
    for h in soft_hits:
        gene = h["gene"]
        if gene not in soft_hits_by_gene or h.get("cov_pct", 0) > soft_hits_by_gene[gene].get("cov_pct", 0):
            soft_hits_by_gene[gene] = h

    # Score each type
    scored = []
    for type_rule in rules["sccmec_type_rules"]:
        if type_rule["name"] == "n/a (Plasmid)":
            continue  # Skip plasmid — not a genomic SCCmec type

        expected_genes = _get_expected_genes(type_rule, rules)
        if not expected_genes:
            continue

        total_weight = sum(GENE_WEIGHTS.get(g, DEFAULT_WEIGHT) for g in expected_genes)
        total_contribution = 0.0
        gene_evidence = []

        for gene in expected_genes:
            contribution, gene_status, id_pct, cov_pct = _compute_gene_contribution(
                gene, hits_by_gene, soft_hits_by_gene
            )
            total_contribution += contribution
            gene_evidence.append({
                "gene": gene,
                "status": gene_status,
                "id_pct": round(id_pct, 1),
                "cov_pct": round(cov_pct, 1),
                "weight": GENE_WEIGHTS.get(gene, DEFAULT_WEIGHT),
                "contribution": round(contribution, 3),
            })

        score = total_contribution / total_weight if total_weight > 0 else 0

        # Complex-match bonus: if the result's already-classified complex
        # matches this type's expected complex, boost the score. This rewards
        # types whose mec/ccr complex aligns with what the classifier found.
        detected_mec = result.get("mec_complex", "")
        detected_ccr = result.get("ccr_complex", "")
        mec_match = (detected_mec != "Negative" and detected_mec == type_rule["mec"])
        ccr_match = (detected_ccr != "Negative" and detected_ccr == type_rule.get("ccr", ""))
        if mec_match:
            score = score * 0.7 + 0.3  # 30% bonus for mec complex match
        if ccr_match:
            score = score * 0.7 + 0.3  # 30% bonus for ccr complex match

        if score <= 0:
            continue

        candidate_status = _determine_candidate_status(gene_evidence)
        ruling = _generate_ruling(type_rule, gene_evidence, result)

        # Component matching summary
        mec_genes = [g for g in gene_evidence if not g["gene"].startswith("ccr")]
        ccr_genes = [g for g in gene_evidence if g["gene"].startswith("ccr")]

        mec_status = "matched" if all(g["status"] == "present" for g in mec_genes) else "partial" if any(g["status"] != "absent" for g in mec_genes) else "missing"
        ccr_status = "matched" if all(g["status"] == "present" for g in ccr_genes) else "partial" if any(g["status"] != "absent" for g in ccr_genes) else "missing"

        mec_detail = ", ".join(f"{g['gene']} {g['status']}" for g in mec_genes)
        ccr_detail = ", ".join(f"{g['gene']} {g['status']}" for g in ccr_genes)

        scored.append({
            "type": type_rule["name"],
            "designation": type_rule.get("designation", ""),
            "score": round(score, 3),
            "status": candidate_status,
            "matching_components": {
                "mec_complex": {
                    "expected": type_rule["mec"],
                    "status": mec_status,
                    "detail": mec_detail,
                },
                "ccr_complex": {
                    "expected": type_rule.get("ccr", "ANY"),
                    "status": ccr_status,
                    "detail": ccr_detail,
                },
            },
            "gene_evidence": gene_evidence,
            "ruling": ruling,
        })

    # Sort by score descending
    scored.sort(key=lambda c: c["score"], reverse=True)
    candidates = scored[:max_candidates]

    if not candidates:
        return None

    # Collect sub-threshold genes used in evidence
    sub_threshold_genes = set()
    for c in candidates:
        for g in c["gene_evidence"]:
            if g["status"] == "sub_threshold":
                sub_threshold_genes.add(g["gene"])

    best = candidates[0]
    return {
        "best_guess": best["type"],
        "best_guess_designation": best["designation"],
        "best_guess_score": best["score"],
        "best_guess_ruling": best["ruling"],
        "n_candidates": len(candidates),
        "method": "weighted_component_matching",
        "sub_threshold_genes_used": sorted(sub_threshold_genes),
        "note": "Estimation based on partial evidence; not a definitive classification",
        "candidates": candidates,
    }
