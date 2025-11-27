#!/usr/bin/env python3
"""
Script to remove redundancy from predicted AMPs using CD-HIT.
Generates a non-redundant FASTA file and parses cluster information.
"""

import subprocess
from pathlib import Path

# =======================
# User-configurable paths
# =======================
#input_fasta = "../data/database_amps/DATABASE_AMPs_clean.fasta"
input_fasta = "../data/predicted_amps/predicted_precursors_with_flag.fasta"  # Path to input FASTA with predicted AMPs
output_dir = "../data/non_redundant_datasets"                        # Directory where CD-HIT output will be saved
cdhit_threshold = 0.95                        # Identity threshold for CD-HIT clustering

# Output filenames

cdhit_output_fasta = Path(output_dir) / "predicted_AMPs_nr.fasta"
cdhit_cluster_file = Path(output_dir) / "predicted_AMPs_nr.fasta.clstr"

"""
cdhit_output_fasta = Path(output_dir) / "database_AMPs_nr.fasta"
cdhit_cluster_file = Path(output_dir) / "database_AMPs_nr.fasta.clstr"
"""


# =====================================
# Step 1: Run CD-HIT to cluster sequences
# =====================================
print("Running CD-HIT...")
subprocess.run([
    "cd-hit", 
    "-i", input_fasta,
    "-o", str(cdhit_output_fasta),
    "-c", str(cdhit_threshold),
    "-n", "5"  # Word length recommended for 0.9-1.0 threshold; adjust if needed
], check=True)

print(f"CD-HIT finished. Non-redundant FASTA saved to {cdhit_output_fasta}")
print(f"Cluster information saved to {cdhit_cluster_file}")

# Número de secuencias representativas

num_representative_seqs = 0

with open(cdhit_output_fasta, "r") as fasta:
    for line in fasta:
        if line.startswith(">"):
            num_representative_seqs += 1

print(f"Number of representative sequences: {num_representative_seqs}")
