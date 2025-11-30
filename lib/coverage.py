import csv
from collections import defaultdict

def calculate_gene_coverage(paf_file, min_coverage=0.90):
    """
    Parses a PAF file where Target = Genes and Query = Reads.
    Calculates the breadth of coverage for each gene.
    Returns a list of "hits" (genes with > min_coverage).
    """
    gene_intervals = defaultdict(list)
    gene_lengths = {}
    
    try:
        with open(paf_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 12: continue
                
                # In Read Mode:
                # Query = Read
                # Target = Gene
                
                target_name = row[5] # Gene name in DB
                target_len = int(row[6])
                target_start = int(row[7])
                target_end = int(row[8])
                
                # Store gene length (should be consistent)
                gene_lengths[target_name] = target_len
                
                # Store interval [start, end]
                gene_intervals[target_name].append((target_start, target_end))
                
    except FileNotFoundError:
        return []

    hits = []
    
    for gene, intervals in gene_intervals.items():
        # Merge overlapping intervals
        intervals.sort(key=lambda x: x[0])
        
        merged = []
        if intervals:
            curr_start, curr_end = intervals[0]
            for next_start, next_end in intervals[1:]:
                if next_start < curr_end: # Overlap
                    curr_end = max(curr_end, next_end)
                else:
                    merged.append((curr_start, curr_end))
                    curr_start, curr_end = next_start, next_end
            merged.append((curr_start, curr_end))
            
        # Calculate total covered length
        covered_len = sum(end - start for start, end in merged)
        total_len = gene_lengths[gene]
        
        breadth = covered_len / total_len if total_len > 0 else 0
        
        if breadth >= min_coverage:
            # Parse gene name for classification compatibility
            # DB headers: Gene__Accession|Type
            if "__" in gene:
                parts = gene.split("__", 1)
                gene_short = parts[0]
                meta = parts[1]
            else:
                parts = gene.split(" ", 1)
                gene_short = parts[0]
                meta = parts[1] if len(parts) > 1 else ""

            if "|" in meta:
                accession, scc_type = meta.split("|", 1)
            else:
                accession, scc_type = meta, "Unknown"

            hits.append({
                "gene": gene_short,
                "accession": accession,
                "scc_type": scc_type,
                "contig": "Reads", # Placeholder
                "start": 0,
                "end": total_len,
                "strand": "N/A",
                "id_pct": 0, # Not easily calculated from aggregate
                "cov_pct": breadth * 100,
                "len": total_len,
                "aln_len": covered_len,
                "breadth": breadth
            })
            
    return hits
