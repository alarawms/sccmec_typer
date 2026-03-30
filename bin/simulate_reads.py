#!/usr/bin/env python3
import random
import argparse
import sys

def simulate_reads(fasta_file, output_fastq, coverage=20, mean_length=5000, error_rate=0.05):
    """
    Simulates Nanopore-like reads from a FASTA file.
    """
    # Read genome
    genome = ""
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                genome += line.strip()
    
    genome_len = len(genome)
    total_bases = genome_len * coverage
    generated_bases = 0
    
    print(f"Simulating reads for {genome_len}bp genome with {coverage}x coverage...")
    
    with open(output_fastq, 'w') as out:
        read_count = 0
        while generated_bases < total_bases:
            read_count += 1
            # Pick length (approximate normal distribution around mean, but simpler)
            read_len = int(random.gauss(mean_length, mean_length/4))
            read_len = max(500, min(read_len, len(genome)))
            
            # Pick start
            start = random.randint(0, genome_len - read_len)
            seq = list(genome[start : start + read_len])
            
            # Introduce errors
            # Simple error model: 
            # - Mismatch
            # - Insertion
            # - Deletion
            final_seq = []
            qual = []
            
            i = 0
            while i < len(seq):
                r = random.random()
                if r < error_rate:
                    # Error
                    err_type = random.choice(['sub', 'ins', 'del'])
                    if err_type == 'sub':
                        final_seq.append(random.choice('ACGT'))
                        qual.append('!') # Low quality
                        i += 1
                    elif err_type == 'ins':
                        final_seq.append(random.choice('ACGT'))
                        qual.append('!')
                        # Don't increment i, we inserted
                        # But to avoid infinite loops in this simple logic, just insert and move on
                        # Actually, insertion means we add a base and don't consume source
                        # But let's just append and continue
                    elif err_type == 'del':
                        i += 1 # Skip source base
                else:
                    final_seq.append(seq[i])
                    qual.append('I') # High quality (approx Q40)
                    i += 1
            
            # Write FASTQ
            header = f"@SIM_READ_{read_count} start={start} len={len(final_seq)}"
            out.write(f"{header}\n")
            out.write("".join(final_seq) + "\n")
            out.write("+\n")
            out.write("".join(qual) + "\n")
            
            generated_bases += len(final_seq)
            
    print(f"Generated {read_count} reads in {output_fastq}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate Nanopore reads")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ")
    parser.add_argument("--cov", type=int, default=20, help="Coverage")
    args = parser.parse_args()
    
    simulate_reads(args.input, args.output, coverage=args.cov)
