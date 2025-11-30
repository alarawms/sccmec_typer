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
    
    # 1. Alignment
    paf_file = f"{args.output}.paf"
    run_minimap2(args.input, args.db, paf_file, args.threads)
    
    # 2. Parsing
    hits = parse_paf(paf_file)
    
    # 3. Classification
    result = classify_sccmec(hits)
    
    # Output results
    import json
    import csv

    # JSON Output
    json_output_file = f"{args.output}.json" if args.output else None
    if json_output_file:
        with open(json_output_file, 'w') as f:
            json.dump(result, f, indent=4)
        print(f"JSON results written to {json_output_file}")
    else:
        print(json.dumps(result, indent=4))

    # TSV Output (Summary)
    if args.output:
        tsv_output_file = f"{args.output}.tsv"
        with open(tsv_output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            # Header
            writer.writerow(['Sample', 'Status', 'SCCmec_Type', 'Mec_Complex', 'Ccr_Complex', 'Genes_Detected', 'Warnings'])
            # Data
            sample_name = os.path.basename(args.input)
            genes_str = ",".join(result.get('genes_detected', []))
            warnings_str = "; ".join(result.get('warnings', []))
            writer.writerow([
                sample_name,
                result.get('status', 'Unknown'),
                result.get('sccmec_type', 'Unknown'),
                result.get('mec_complex', 'Unknown'),
                result.get('ccr_complex', 'Unknown'),
                genes_str,
                warnings_str
            ])
        print(f"TSV summary written to {tsv_output_file}")

if __name__ == "__main__":
    main()
