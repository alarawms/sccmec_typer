import json
import os

def load_rules():
    """Load classification rules from JSON file."""
    rules_path = os.path.join(os.path.dirname(__file__), '../db/rules.json')
    try:
        with open(rules_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading rules: {e}")
        return None

def classify_sccmec(hits_list):
    """
    Takes parsed hits (list of dicts) and determines the SCCmec type using rules.json.
    """
    if not hits_list:
        return {"status": "Negative", "reason": "No alignments found"}
    
    rules = load_rules()
    if not rules:
        return {"status": "Error", "reason": "Could not load classification rules"}

    # Get unique genes found
    genes_found = set(h["gene"] for h in hits_list)
    
    # Check for split assembly
    contigs_found = set(h.get("contig", "unknown") for h in hits_list)
    warnings = []
    if len(contigs_found) > 1:
        warnings.append(f"Split Assembly: Components found on {len(contigs_found)} different contigs: {', '.join(contigs_found)}")

    # 1. Identify mec Complex
    mec_complex = "Negative"
    
    for rule in rules["mec_complex_rules"]:
        # Check required genes
        required_met = all(g in genes_found for g in rule["required"])
        
        # Check 'any_of' genes (if present)
        any_of_met = True
        if "any_of" in rule:
            any_of_met = any(g in genes_found for g in rule["any_of"])
            
        if required_met and any_of_met:
            mec_complex = rule["name"]
            break # Priority is determined by order in JSON
            
    # 2. Identify ccr Complex
    ccr_genes = [g for g in genes_found if g.startswith("ccr")]
    ccr_types = set()
    
    for ccr in ccr_genes:
        for rule in rules["ccr_complex_rules"]:
            if rule["match_substring"] in ccr:
                ccr_types.add(rule["name"])
    
    ccr_complex = " / ".join(sorted(ccr_types)) if ccr_types else "Negative"
    
    # 3. Determine SCCmec Type & Status
    sccmec_type = "Unknown"
    status = "Positive"

    if mec_complex == "Negative" and ccr_complex == "Negative":
        return {"status": "Negative", "reason": "No mec or ccr genes found", "genes_detected": list(genes_found)}
    
    if mec_complex == "Negative":
        status = "Partial (Orphan ccr)"
        warnings.append("Found ccr genes but no mecA/mecB/mecC")
    elif ccr_complex == "Negative":
        # Check if this mec complex allows missing ccr (e.g. Plasmid)
        # We check the sccmec_type_rules for a match where ccr is "ANY" or ignored
        matched_special = False
        for rule in rules["sccmec_type_rules"]:
            if rule["mec"] == mec_complex and rule.get("ccr") == "ANY":
                sccmec_type = rule["name"]
                status = "Positive"
                matched_special = True
                break
        
        if not matched_special:
            status = "Partial (Unclassifiable)"
            warnings.append("Found mecA/mecC but no ccr genes")
    
    # Standard combinations
    if status == "Positive" and sccmec_type == "Unknown":
        for rule in rules["sccmec_type_rules"]:
            if rule["mec"] == mec_complex:
                # Check ccr match
                target_ccr = rule["ccr"]
                if target_ccr == "ANY":
                    sccmec_type = rule["name"]
                    break
                elif target_ccr in ccr_types:
                    sccmec_type = rule["name"]
                    break
        
    # Handling composites or novel types
    if len(ccr_types) > 1:
        # If we already found a type (e.g. Type XI which has multiple ccr), we might not want to prefix Composite?
        # But Type XI rule says ccr="ANY".
        # If it's Type XI, we might want to leave it.
        # But for others, e.g. Type IV + Type 5 ccr -> Composite (Type IV).
        if sccmec_type != "Type XI": # Hardcoded exception for now, or add to rules?
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
