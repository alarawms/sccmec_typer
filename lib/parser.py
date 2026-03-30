import csv

SOFT_HIT_MIN_IDENTITY = 0.70
SOFT_HIT_MIN_COVERAGE = 0.50

def parse_paf(paf_file):
    """
    Parses a PAF file and filters hits.
    Returns a (hits, soft_hits) tuple where:
      - hits: identity >= 90% AND coverage >= 80%
      - soft_hits: identity >= SOFT_HIT_MIN_IDENTITY AND coverage >= SOFT_HIT_MIN_COVERAGE
                   but below the hits thresholds
    """
    hits = []
    soft_hits = []
    
    try:
        with open(paf_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 12: continue
                
                # Extract fields
                query_name = row[0]
                query_len = int(row[1])
                query_start = int(row[2])
                query_end = int(row[3])
                strand = row[4]
                target_name = row[5]
                target_len = int(row[6])
                target_start = int(row[7])
                target_end = int(row[8])
                matches = int(row[9])
                block_len = int(row[10])
                
                # Calculate identity & coverage
                # Identity = matches / block_len
                identity = matches / block_len if block_len > 0 else 0
                
                # Coverage of the GENE (Query)
                coverage = (query_end - query_start) / query_len if query_len > 0 else 0
                
                # Filter: check if it meets hard or soft thresholds
                is_hard = identity >= 0.90 and coverage >= 0.80
                is_soft = (not is_hard
                           and identity >= SOFT_HIT_MIN_IDENTITY
                           and coverage >= SOFT_HIT_MIN_COVERAGE)

                if is_hard or is_soft:
                    # Parse query header (Gene__Accession|Type)
                    # The database headers are now >Gene__Accession|Type to ensure unique IDs for minimap2
                    if "__" in query_name:
                        parts = query_name.split("__", 1)
                        gene = parts[0]
                        meta = parts[1]
                    else:
                        # Fallback for old format or unexpected headers
                        parts = query_name.split(" ", 1)
                        gene = parts[0]
                        meta = parts[1] if len(parts) > 1 else ""

                    if "|" in meta:
                        accession, scc_type = meta.split("|", 1)
                    else:
                        accession, scc_type = meta, "Unknown"

                    hit = {
                        "gene": gene,
                        "accession": accession,
                        "scc_type": scc_type,
                        "identity": identity,
                        "coverage": coverage,
                        "contig": target_name,
                        "start": target_start,
                        "end": target_end,
                        "strand": strand,
                        "len": query_len,
                        "aln_len": block_len,
                        "id_pct": identity * 100,
                        "cov_pct": coverage * 100
                    }

                    if is_hard:
                        hits.append(hit)
                    else:
                        soft_hits.append(hit)

    except FileNotFoundError:
        return [], []

    return hits, soft_hits
