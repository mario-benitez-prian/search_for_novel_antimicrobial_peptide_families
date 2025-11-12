#!/usr/bin/env python3
"""
run_interproscan.py
Combine database and predicted AMPs, then classify them into families using InterProScan.
"""

import subprocess
from pathlib import Path

# -----------------------------
# CONFIGURATION
# -----------------------------
INTERPROSCAN_DIR = Path.home() / "software/interproscan-5.75-106.0"
DATABASE_FASTA = Path("../data/database_amps/DATABASE_AMPs_clean.fasta")
PREDICTED_FASTA = Path("../data/predicted_amps/predicted_precursors_with_flag.fasta")
OUTPUT_DIR = Path("../results/interpro_results")
OUTPUT_PREFIX = OUTPUT_DIR / "predictions_and_database_families_interpro"

# Create output directory if it doesn't exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# COMBINE FASTA FILES
# -----------------------------
COMBINED_FASTA = OUTPUT_DIR / "combined_database_predicted.fasta"

print("🧬 Combining FASTA files...")
with open(COMBINED_FASTA, "w") as outfile:
    for fasta_file in [DATABASE_FASTA, PREDICTED_FASTA]:
        with open(fasta_file, "r") as infile:
            outfile.write(infile.read().rstrip("\n") + "\n")

print(f"✅ Combined FASTA written to: {COMBINED_FASTA}")

# -----------------------------
# RUN INTERPROSCAN
# -----------------------------
cmd = [
    str(INTERPROSCAN_DIR / "interproscan.sh"),
    "-i", str(COMBINED_FASTA),
    "-f", "tsv",
    "-dp",        # enable domain/family annotation
    "-goterms",   # include GO terms
    "-pathways",
    "-iprlookup", # include InterPro annotations
    "-appl", "AntiFam, CDD, FunFam, Gene3D, Hamap, NCBIfam, PANTHER, PIRSF, PRINTS, Pfam, SMART, SUPERFAMILY",
    "-b", str(OUTPUT_PREFIX)
]

print("⚡ Running InterProScan...")
try:
    subprocess.run(cmd, check=True)
    print(f"✅ InterProScan finished successfully.\nResults saved in {OUTPUT_PREFIX}.tsv")
except subprocess.CalledProcessError as e:
    print("❌ InterProScan failed:")
    print(e)

