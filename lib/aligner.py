import subprocess
import sys
import os

def run_minimap2(input_file, db_file, output_paf, threads=4):
    """
    Runs minimap2 to align input sequences against the reference database.
    """
    # Determine path to minimap2
    # Check if it's in the PATH
    minimap2_cmd = "minimap2"
    
    # Check if it's in the bin directory relative to this file
    # lib/aligner.py -> ../bin/minimap2
    local_bin = os.path.join(os.path.dirname(os.path.dirname(__file__)), "bin", "minimap2")
    if os.path.exists(local_bin):
        minimap2_cmd = local_bin

    cmd = [
        minimap2_cmd,
        "-cx", "asm5", # Optimized for assembly-to-assembly
        "-t", str(threads),
        input_file, # Target (Genome)
        db_file     # Query (Genes)
    ]
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        with open(output_paf, "w") as out_f:
            subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2: {e.stderr.decode()}")
        sys.exit(1)
