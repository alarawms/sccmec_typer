import json
import os


def load_rules():
    """Load classification rules from JSON file."""
    rules_path = os.path.join(os.path.dirname(__file__), "../db/rules.json")
    try:
        with open(rules_path, "r") as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading rules: {e}")
        return None


IS431_NAMES = ("IS431", "IS431_1", "IS431_2")


def _flanking_is431_orientation(hits_list):
    """C1/C2 from the IS431 copies that flank the mec gene (IS431-mecA-dmecR1-IS431).

    Uses the nearest IS431 locus on each side of mecA on the same contig, so extra IS431
    copies elsewhere (common in composites) cannot flip the call. Returns "same",
    "opposite", or None when coordinates are missing or one side has no IS431.
    """
    mec = [h for h in hits_list if h["gene"] in MEC_GENES and h.get("contig_len")]
    if not mec:
        return None
    m = max(mec, key=lambda h: h.get("id_pct", 0))
    iss = [h for h in hits_list if h["gene"] in IS431_NAMES and h.get("contig") == m["contig"]
           and h.get("strand") in ("+", "-")]
    left = [h for h in iss if h["end"] <= m["start"] + 100]
    right = [h for h in iss if h["start"] >= m["end"] - 100]
    if not left or not right:
        return None
    # nearest locus on each side; within a locus several database IS431 entries can align
    # in different orientations, so take the longest (full-length) alignment there
    near_l = max(h["end"] for h in left)
    near_r = min(h["start"] for h in right)
    span = lambda h: h.get("aln_len") or (h["end"] - h["start"])
    a = max((h for h in left if h["end"] >= near_l - 300), key=span)
    b = max((h for h in right if h["start"] <= near_r + 300), key=span)
    return "same" if a["strand"] == b["strand"] else "opposite"


def _determine_is431_orientation(hits_list):
    """Determine IS431 orientation from PAF strand information.

    Compares strands of IS431_1 and IS431_2 hits to distinguish
    mec complex Class C1 (same orientation) from C2 (opposite orientation).

    Returns:
        "same"     - both IS431 copies on same strand
        "opposite" - IS431 copies on different strands
        None       - cannot determine (only one copy, or no numbered variants)
    """
    is431_strands = {}
    for h in hits_list:
        gene = h["gene"]
        if gene in ("IS431_1", "IS431_2"):
            is431_strands[gene] = h.get("strand", "+")

    if "IS431_1" in is431_strands and "IS431_2" in is431_strands:
        s1 = is431_strands["IS431_1"]
        s2 = is431_strands["IS431_2"]
        # Guard against reads mode where strand is "N/A"
        if s1 not in ("+", "-") or s2 not in ("+", "-"):
            return None
        if s1 == s2:
            return "same"
        else:
            return "opposite"

    return None


def _classify_mec_complex(genes_found, hits_list, rules):
    """Determine mec complex class using priority-ordered rules.

    The rules are evaluated in order. The first matching rule wins.
    Rules may specify:
      - required: all these genes must be present
      - any_of: at least one of these genes must be present
      - excluded: none of these genes may be present
      - orientation: IS431 orientation must match ("same" or "opposite")
    """
    is431_orientation = (_flanking_is431_orientation(hits_list)
                         or _determine_is431_orientation(hits_list))

    for rule in rules["mec_complex_rules"]:
        # Check required genes
        if not all(g in genes_found for g in rule["required"]):
            continue

        # Check any_of genes
        if "any_of" in rule:
            if not any(g in genes_found for g in rule["any_of"]):
                continue

        # Check excluded genes
        if "excluded" in rule:
            if any(g in genes_found for g in rule["excluded"]):
                continue

        # Check IS431 orientation (for C1/C2 distinction)
        if "orientation" in rule:
            if is431_orientation is None:
                # Cannot determine orientation; skip this specific rule,
                # the generic "Class C" fallback will catch it
                continue
            if is431_orientation != rule["orientation"]:
                continue

        return rule["name"]

    return "Negative"


# A mec gene this close to a contig end means the cassette was likely split by the
# assembler (IS431 repeats flank mecA); class D / missing-ccr calls are then unreliable.
CONTIG_EDGE_BP = 1500

MEC_GENES = ("mecA", "mecC", "mecB")

# mec-complex accessory genes only define the class when they sit in the mec locus.
# IS431 / IS1272 also occur elsewhere in S. aureus genomes; genome-wide presence gave
# false class B / C2 calls on fragmented assemblies.
MEC_ACCESSORY = ("IS431", "IS431_1", "IS431_2", "IS1272", "mecR1", "mecI")
MEC_LOCUS_BP = 15000   # IS431 can sit ~11 kb downstream of mecA (class C2)
UNDETERMINED_CLASS = "Undetermined (fragmented)"


def _mec_locus_hits(hits_list):
    """Hits used for the mec complex: mec genes + accessory genes within MEC_LOCUS_BP
    of a mec gene on the same contig. Read-mode hits (no contig length) are kept as is."""
    mec = [h for h in hits_list if h["gene"] in MEC_GENES]
    if not mec or not all(h.get("contig_len") for h in mec):
        return hits_list
    def near(h):
        return any(h.get("contig") == m.get("contig")
                   and h.get("start", 0) < m.get("end", 0) + MEC_LOCUS_BP
                   and h.get("end", 0) > m.get("start", 0) - MEC_LOCUS_BP for m in mec)
    return [h for h in hits_list if h["gene"] not in MEC_ACCESSORY or near(h)]


def collapse_loci(hits_list):
    """One hit per locus: overlapping hits of the same gene on the same contig (different
    database alleles) collapse to the best one (identity x coverage). Read-mode hits without
    coordinates are returned unchanged."""
    if not hits_list or not all("start" in h and "end" in h and h.get("contig_len") for h in hits_list):
        return list(hits_list)
    best = []
    for h in sorted(hits_list, key=lambda h: (h["gene"], h.get("contig", ""), h["start"])):
        last = best[-1] if best else None
        if last and last["gene"] == h["gene"] and last.get("contig") == h.get("contig") and h["start"] <= last["end"]:
            score = lambda x: x.get("id_pct", 0) * x.get("cov_pct", 0)
            if score(h) > score(last):
                h = dict(h, start=min(h["start"], last["start"]), end=max(h["end"], last["end"]))
                best[-1] = h
            else:
                best[-1] = dict(last, end=max(last["end"], h["end"]))
        else:
            best.append(dict(h))
    return best


def _locus_count(hits_list, gene):
    """Number of distinct (non-overlapping) loci at which `gene` was hit.

    Several database alleles often hit the same locus; those collapse to one.
    Read-mode hits have no coordinates and count as a single locus.
    """
    spans = sorted((h.get("contig", ""), h.get("start", 0), h.get("end", 0))
                   for h in hits_list if h["gene"] == gene)
    loci, last = 0, None
    for contig, start, end in spans:
        if last is None or contig != last[0] or start > last[2]:
            loci += 1
            last = (contig, start, end)
        else:
            last = (contig, last[1], max(last[2], end))
    return loci


def _ccr_copies(hits_list, ccr_type, rules):
    """Copies of a ccr complex = fewest loci among its required genes."""
    rule = next(r for r in rules["ccr_complex_rules"] if r["name"] == ccr_type)
    return min(_locus_count(hits_list, g) for g in rule["required_genes"])


def _mec_at_contig_edge(hits_list):
    """True if any mec gene hit lies within CONTIG_EDGE_BP of its contig end."""
    for h in hits_list:
        if h["gene"] in MEC_GENES and h.get("contig_len"):
            if h["start"] < CONTIG_EDGE_BP or h["contig_len"] - h["end"] < CONTIG_EDGE_BP:
                return True
    return False


def _iwg_designations(mec_complex, ccr_types, hits_list, rules):
    """IWG-SCC style type/designation, including composite elements.

    Notation follows the literature convention for extra ccr complexes, e.g.
    V(5C2&5) = type V with a second ccrC (ccr5) copy, VII(5C1&1) = type VII plus ccr1.
    When several defined types fit (e.g. class B with ccr1 and ccr2 -> I or IV),
    all are listed instead of picking the first rule.

    Returns (iwg_type, iwg_designation, type_candidates).
    """
    if mec_complex in ("Negative", "Plasmid-borne (mecB)") or not ccr_types:
        return "nt", None, []
    num0 = lambda t: t.replace("Type ", "")
    if mec_complex == UNDETERMINED_CLASS:      # ccr known, mec class not (split cassette)
        d = "&".join(num0(t) for t in sorted(ccr_types)) + "?"
        return f"nt({d})", d, []
    mec_ok = [mec_complex] + (["Class C1", "Class C2"] if mec_complex == "Class C" else [])
    num = lambda t: t.replace("Type ", "")
    copies = {t: _ccr_copies(hits_list, t, rules) for t in ccr_types}
    labels, desigs, cands = [], [], []
    for rule in rules["sccmec_type_rules"]:
        if rule["mec"] not in mec_ok or rule["ccr"] not in ccr_types:
            continue
        extras = [num(t) for t in sorted(ccr_types) if t != rule["ccr"]]
        extras += [num(rule["ccr"])] * (copies[rule["ccr"]] - 1)
        d = rule["designation"] + "".join(f"&{e}" for e in sorted(extras))
        roman = rule["name"].replace("Type ", "")
        labels.append(f"{roman}({d})")
        desigs.append(d)
        cands.append(rule["name"])
    if not cands:  # mec + ccr present but no defined type combines them
        mec_letter = mec_complex.replace("Class ", "")
        d = "&".join(num(t) for t in sorted(ccr_types)) + mec_letter
        return f"nt({d})", d, []
    return " / ".join(labels), " / ".join(desigs), cands


def _read_mode_classes(genes_found, rules):
    """mec classes whose required/any_of genes are present, ignoring exclusions and IS431
    orientation. Read-mode hits have no coordinates, so an IS1272/IS431 copy anywhere in
    the genome would otherwise force (or exclude) a class. C1/C2 collapse to 'Class C'."""
    out = []
    for rule in rules["mec_complex_rules"]:
        if not all(g in genes_found for g in rule["required"]):
            continue
        if "any_of" in rule and not any(g in genes_found for g in rule["any_of"]):
            continue
        name = "Class C" if rule["name"] in ("Class C1", "Class C2") else rule["name"]
        if name not in out:
            out.append(name)
    return out


def _has_defined_type(mec_class, ccr_types, rules):
    mecs = [mec_class] + (["Class C1", "Class C2"] if mec_class == "Class C" else [])
    return any(r["mec"] in mecs and r["ccr"] in ccr_types for r in rules["sccmec_type_rules"])


def _classify_ccr_complex(genes_found, rules):
    """Determine ccr complex type(s) using pair-based gene matching.

    Each ccr rule specifies required_genes (e.g., ["ccrA1", "ccrB6"] for Type 7).
    A ccr type matches only when ALL its required genes are present.

    Pair-based rules are checked first (ccrAx + ccrBx), then single-gene rules
    (ccrC1, ccrC2). When a pair matches, the individual genes are "consumed"
    so they don't also match a single-component rule.

    Returns a set of matched ccr type names.
    """
    ccr_genes = {g for g in genes_found if g.startswith("ccr")}
    if not ccr_genes:
        return set()

    ccr_types = set()
    consumed_genes = set()

    # Sort rules: multi-gene pairs first, then single-gene rules
    pair_rules = [r for r in rules["ccr_complex_rules"] if len(r["required_genes"]) > 1]
    single_rules = [r for r in rules["ccr_complex_rules"] if len(r["required_genes"]) == 1]

    # Match pair rules first
    for rule in pair_rules:
        required = set(rule["required_genes"])
        available = ccr_genes - consumed_genes
        if required.issubset(available):
            ccr_types.add(rule["name"])
            consumed_genes.update(required)

    # Match single-gene rules (ccrC1, ccrC2) only if gene not consumed
    for rule in single_rules:
        required = set(rule["required_genes"])
        available = ccr_genes - consumed_genes
        if required.issubset(available):
            ccr_types.add(rule["name"])
            consumed_genes.update(required)

    return ccr_types


def classify_sccmec(hits_list):
    """Classify SCCmec type from parsed alignment hits.

    Takes parsed hits (list of dicts with at minimum a "gene" key) and
    determines the SCCmec type using the IWG-SCC classification scheme
    defined in db/rules.json.

    Algorithm:
    1. Identify mec complex class (A, B, C1, C2, D, E) from gene presence
       and IS431 orientation
    2. Identify ccr complex type(s) (1-9) from ccr gene pairs
    3. Match mec+ccr combination to one of the 15 approved SCCmec types
    4. Handle edge cases: composite, orphan ccr, partial, negative
    """
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

    rules = load_rules()
    if not rules:
        return {"status": "Error", "reason": "Could not load classification rules"}

    # Unique genes detected
    genes_found = set(h["gene"] for h in hits_list)

    # Split assembly check
    contigs_found = set(h.get("contig", "unknown") for h in hits_list)
    warnings = []
    if len(contigs_found) > 1:
        warnings.append(
            f"Split Assembly: Components found on {len(contigs_found)} "
            f"different contigs: {', '.join(sorted(contigs_found))}"
        )

    # mecA/mecC/mecB presence flags
    mecA_present = "mecA" in genes_found
    mecC_present = "mecC" in genes_found
    mecB_present = "mecB" in genes_found

    # Step 1: mec complex — only from genes in the mec locus (see MEC_LOCUS_BP)
    locus_hits = _mec_locus_hits(hits_list)
    locus_genes = set(h["gene"] for h in locus_hits)
    mec_complex = _classify_mec_complex(locus_genes, locus_hits, rules)
    if locus_genes != genes_found & (locus_genes | set(MEC_ACCESSORY)):
        ignored = sorted((genes_found & set(MEC_ACCESSORY)) - locus_genes)
        if ignored:
            warnings.append(f"Ignored for mec class (not within {MEC_LOCUS_BP} bp of mec gene): "
                            f"{', '.join(ignored)}")
    # A mec gene at a contig end with no class-defining neighbours: the class is unknown,
    # not D ("no IS") — the flanking IS431/IS1272 is simply on another contig.
    mec_split = _mec_at_contig_edge(hits_list)
    if mec_split and mec_complex in ("Class D", "Unclassifiable (mecA)"):
        mec_complex = UNDETERMINED_CLASS

    # Step 2: ccr complex
    ccr_types = _classify_ccr_complex(genes_found, rules)
    ccr_complex = " / ".join(sorted(ccr_types)) if ccr_types else "Negative"

    # Read mode: if the first-matching mec class forms no defined type with the ccr found,
    # but another class whose markers are also present does, use that one.
    read_mode = not any(h.get("contig_len") for h in hits_list)
    if read_mode and ccr_types and mec_complex not in ("Negative",) \
            and not _has_defined_type(mec_complex, ccr_types, rules):
        for alt in _read_mode_classes(genes_found, rules):
            if alt != mec_complex and _has_defined_type(alt, ccr_types, rules):
                warnings.append(
                    f"Read mode: {mec_complex} + {ccr_complex} is not a defined type; using {alt} "
                    f"(its markers are also present; reads cannot localise IS431/IS1272 to the mec locus)")
                mec_complex = alt
                break

    # Step 3: SCCmec type assignment
    sccmec_type = "Unknown"
    status = "Positive"

    if mec_complex == "Negative" and not ccr_types:
        return {
            "status": "Negative",
            "reason": "No mec or ccr genes found",
            "sccmec_type": "Negative",
            "mec_complex": "Negative",
            "ccr_complex": "Negative",
            "genes_detected": sorted(genes_found),
            "mecA_present": mecA_present,
            "warnings": [],
            "hits_summary": hits_list,
        }

    if mec_complex == UNDETERMINED_CLASS and ccr_types:
        status = "Partial (Assembly-limited)"
    elif mec_complex == "Negative":
        status = "Partial (Orphan ccr)"
        warnings.append("Found ccr genes but no mecA/mecB/mecC")
    elif not ccr_types:
        # Check if mec complex allows missing ccr (e.g., Plasmid-borne)
        matched_special = False
        for rule in rules["sccmec_type_rules"]:
            if rule["mec"] == mec_complex and rule.get("ccr") == "ANY":
                sccmec_type = rule["name"]
                status = "Positive"
                matched_special = True
                break
        if not matched_special:
            status = "Partial (Unclassifiable)"
            warnings.append("Found mec gene(s) but no ccr genes")

    # Match mec+ccr combination to SCCmec type
    if status == "Positive" and sccmec_type == "Unknown":
        # For mec classes that have C1/C2 fallback to generic C:
        # If mec_complex is "Class C" (orientation undetermined), try matching
        # both C1 and C2 type rules
        mec_candidates = [mec_complex]
        if mec_complex == "Class C":
            mec_candidates = ["Class C1", "Class C2", "Class C"]

        for mec_candidate in mec_candidates:
            for rule in rules["sccmec_type_rules"]:
                if rule["mec"] != mec_candidate:
                    continue
                target_ccr = rule["ccr"]
                if target_ccr == "ANY":
                    sccmec_type = rule["name"]
                    break
                if target_ccr in ccr_types:
                    sccmec_type = rule["name"]
                    # When falling back to generic Class C, report which
                    # specific class the matched type expects
                    if mec_complex == "Class C" and mec_candidate != "Class C":
                        warnings.append(
                            f"IS431 orientation undetermined; type {rule['name']} "
                            f"expects {mec_candidate}"
                        )
                    break
            if sccmec_type != "Unknown":
                break

    iwg_type, iwg_designation, type_candidates = _iwg_designations(
        mec_complex, ccr_types, hits_list, rules)

    # Class C (IS431 orientation unresolved, e.g. read mode): C1 and C2 types both fit;
    # name them all instead of the first one tried.
    if mec_complex == "Class C" and len(ccr_types) == 1 and len(type_candidates) > 1:
        sccmec_type = " / ".join(type_candidates)

    # Composite detection: multiple ccr types. Name every defined type that fits
    # rather than the first rule in db order (class B + ccr1 + ccr2 is I or IV).
    if len(ccr_types) > 1 and type_candidates:
        sccmec_type = f"Composite ({' / '.join(type_candidates)})"
        warnings.append(f"Multiple ccr types detected: {', '.join(sorted(ccr_types))}")
    elif len(ccr_types) > 1:
        sccmec_type = "Composite"
        warnings.append(f"Multiple ccr types detected: {', '.join(sorted(ccr_types))}")

    # Fragmented assemblies: a mec gene at a contig end means the downstream IS431
    # (class C/D distinction) and the ccr genes may simply be on other contigs.
    assembly_limited = mec_split and (not ccr_types or mec_complex in ("Class D", UNDETERMINED_CLASS))
    if not ccr_types:   # mec found, no ccr: distinguish a split cassette from an intact one
        iwg_type = "nt(?)" if assembly_limited else "nt(no ccr)"
    if assembly_limited:
        if not ccr_types and status == "Partial (Unclassifiable)":
            status = "Partial (Assembly-limited)"
        warnings.append(
            f"mec gene within {CONTIG_EDGE_BP} bp of a contig end: cassette likely split "
            f"by the assembly (IS431 repeats); missing ccr / class D may be artefacts. "
            f"Type from reads or a long-read assembly")

    return {
        "status": status,
        "sccmec_type": sccmec_type,
        "iwg_type": iwg_type,
        "iwg_designation": iwg_designation,
        "type_candidates": type_candidates,
        "ccr_copies": {t: _ccr_copies(hits_list, t, rules) for t in sorted(ccr_types)},
        "assembly_limited": assembly_limited,
        "mec_complex": mec_complex,
        "ccr_complex": ccr_complex,
        "genes_detected": sorted(genes_found),
        "mecA_present": mecA_present,
        "warnings": warnings,
        "hits_summary": hits_list,
    }
