#!/usr/bin/env python3
import argparse
import sys
import os

# Add lib to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from lib.aligner import run_minimap2
from lib.parser import parse_paf
from lib.classifier import classify_sccmec

def main():
    parser = argparse.ArgumentParser(description="SCCmec Typer: A minimap2-based tool for typing SCCmec elements.")
    parser.add_argument("-i", "--input", required=True, help="Input genome assembly (FASTA) or long reads (FASTQ)")
    parser.add_argument("-d", "--db", required=True, help="Path to reference database (mfa)")
    parser.add_argument("-o", "--output", default="sccmec_results", help="Output prefix")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    
    args = parser.parse_args()
    
    print(f"Starting SCCmec Typer on {args.input}...")
    
    # Determine input type and preset
    # Default to asm5 for assembly
    preset = "asm5"
    is_reads = False
    
    if args.input.lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
        print("Note: Input appears to be FASTQ. Switching to Read Mode.")
        preset = "map-ont" # Default to ONT, could be sr if user specified but for now auto-detect extension? 
        # Ideally we'd have a flag for short reads vs long reads. 
        # Assuming long reads for now as per "sccmec typing" often uses nanopore.
        is_reads = True

    # 1. Alignment
    paf_file = f"{args.output}.paf"
    
    if is_reads:
        # Read Mode: Map Reads (Query) -> Genes (Target)
        run_minimap2(args.input, args.db, paf_file, args.threads, preset=preset, target_is_db=True)
        
        # 2. Parsing (Coverage Calculation)
        from lib.coverage import calculate_gene_coverage
        hits = calculate_gene_coverage(paf_file)
        
    else:
        # Assembly Mode: Map Genes (Query) -> Genome (Target)
        run_minimap2(args.input, args.db, paf_file, args.threads, preset=preset, target_is_db=False)
        
        # 2. Parsing
        hits = parse_paf(paf_file)
    
    # 3. Classification
    result = classify_sccmec(hits)
    
    # Output results
    import json
    import csv

    # JSON Output
    json_output_file = f"{args.output}.json"
    with open(json_output_file, 'w') as f:
        json.dump(result, f, indent=4)
    print(f"JSON results written to {json_output_file}")

    # CSV Output (Detailed Elements)
    csv_output_file = f"{args.output}_elements.csv"
    with open(csv_output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Get headers from first hit if available, else default
        headers = ["gene", "contig", "start", "end", "strand", "id_pct", "cov_pct", "len", "aln_len"]
        writer.writerow(headers)
        
        for hit in result.get("hits_summary", []):
            writer.writerow([
                hit.get("gene", ""),
                hit.get("contig", ""),
                hit.get("start", ""),
                hit.get("end", ""),
                hit.get("strand", ""),
                f"{hit.get('id_pct', 0):.2f}",
                f"{hit.get('cov_pct', 0):.2f}",
                hit.get("len", ""),
                hit.get("aln_len", "")
            ])
    print(f"Detailed elements CSV written to {csv_output_file}")

    # TSV Output (Summary)
    tsv_output_file = f"{args.output}.tsv"
    with open(tsv_output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Header
        writer.writerow(['Sample', 'Status', 'mecA_Present', 'SCCmec_Type', 'Mec_Complex', 'Ccr_Complex', 'Genes_Detected', 'Warnings'])
        # Data
        sample_name = os.path.basename(args.input)
        genes_str = ",".join(result.get('genes_detected', []))
        warnings_str = "; ".join(result.get('warnings', []))
        writer.writerow([
            sample_name,
            result.get('status', 'Unknown'),
            "Yes" if result.get('mecA_present', False) else "No",
            result.get('sccmec_type', 'Unknown'),
            result.get('mec_complex', 'Unknown'),
            result.get('ccr_complex', 'Unknown'),
            genes_str,
            warnings_str
        ])
    print(f"TSV summary written to {tsv_output_file}")

if __name__ == "__main__":
    main()
