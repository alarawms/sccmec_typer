def classify_sccmec(hits_list):
    """
    Takes parsed hits (list of dicts) and determines the SCCmec type.
    """
    if not hits_list:
        return {"status": "Negative", "reason": "No alignments found"}
    
    # Get unique genes found
    genes_found = set(h["gene"] for h in hits_list)
    
    # Check for split assembly
    contigs_found = set(h.get("contig", "unknown") for h in hits_list)
    warnings = []
    if len(contigs_found) > 1:
        warnings.append(f"Split Assembly: Components found on {len(contigs_found)} different contigs: {', '.join(contigs_found)}")

    # 1. Identify mec Complex
    mec_complex = "Negative"
    has_mecA = "mecA" in genes_found or "mecC" in genes_found
    has_mecR1 = "mecR1" in genes_found
    has_mecI = "mecI" in genes_found
    has_IS1272 = "IS1272" in genes_found
    has_IS431 = "IS431" in genes_found or "IS431_1" in genes_found or "IS431_2" in genes_found
    
    if has_mecA:
        if has_mecR1 and has_mecI:
            mec_complex = "Class A"
        elif has_IS1272: # Class B usually has truncated mecR1 and IS1272
            mec_complex = "Class B"
        elif has_IS431: # Class C usually has IS431
            mec_complex = "Class C"
        else:
            mec_complex = "Class Unclear"
    
    # 2. Identify ccr Complex
    # Look for ccr genes
    ccr_genes = [g for g in genes_found if g.startswith("ccr")]
    ccr_types = set()
    
    for ccr in ccr_genes:
        # Extract type from gene name (e.g., ccrA1 -> 1)
        if "1" in ccr: ccr_types.add("Type 1")
        if "2" in ccr: ccr_types.add("Type 2")
        if "3" in ccr: ccr_types.add("Type 3")
        if "4" in ccr: ccr_types.add("Type 4")
        if "5" in ccr: ccr_types.add("Type 5")
        if "8" in ccr: ccr_types.add("Type 8")
    
    ccr_complex = " / ".join(sorted(ccr_types)) if ccr_types else "Negative"
    
    # 3. Determine SCCmec Type & Status
    sccmec_type = "Unknown"
    status = "Positive"

    if mec_complex == "Negative" and ccr_complex == "Negative":
        return {"status": "Negative", "reason": "No mec or ccr genes found", "genes_detected": list(genes_found)}
    
    if mec_complex == "Negative":
        status = "Partial (Orphan ccr)"
        warnings.append("Found ccr genes but no mecA/mecC")
    elif ccr_complex == "Negative":
        status = "Partial (Unclassifiable)"
        warnings.append("Found mecA/mecC but no ccr genes")
    
    # Standard combinations
    if mec_complex == "Class B" and "Type 1" in ccr_types:
        sccmec_type = "Type I"
    elif mec_complex == "Class A" and "Type 2" in ccr_types:
        sccmec_type = "Type II"
    elif mec_complex == "Class A" and "Type 3" in ccr_types:
        sccmec_type = "Type III"
    elif mec_complex == "Class B" and "Type 2" in ccr_types:
        sccmec_type = "Type IV"
    elif "Class C" in mec_complex and "Type 5" in ccr_types:
        sccmec_type = "Type V"
    elif mec_complex == "Class B" and "Type 4" in ccr_types:
        sccmec_type = "Type VI"
        
    # Handling composites or novel types
    if len(ccr_types) > 1:
        sccmec_type = f"Composite ({sccmec_type})" if sccmec_type != "Unknown" else "Composite"
        warnings.append("Multiple ccr types detected")

    results = {
        "status": status,
        "sccmec_type": sccmec_type,
        "mec_complex": mec_complex,
        "ccr_complex": ccr_complex,
        "genes_detected": list(genes_found),
        "warnings": warnings,
        "hits_summary": hits_list
    }
    
    return results
