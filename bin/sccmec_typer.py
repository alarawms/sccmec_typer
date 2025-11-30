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
    
    print("Typing complete.")
    print(result)

if __name__ == "__main__":
    main()
