#!/usr/bin/env python3
"""
run_interprorx_analysis.py

Full reproducible pipeline:
1. Runs interprorx.py to classify InterProScan results into clusters.
2. Summarizes cluster composition (DATABASE vs PRED).
3. Detects novel AMP clusters (PRED only).
"""

import subprocess
import pandas as pd
from pathlib import Path

# -----------------------------
# CONFIGURATION
# -----------------------------
# Paths to inputs and outputs
INTERPRO_TSV = Path("../results/interpro_results/without_domains_predictions_and_database_families_interpro.tsv")
COMBINED_FASTA = Path("../results/interpro_results/combined_database_predicted.fasta")
OUTDIR = Path("../results/interpro_results")

# interprorx command
INTERPRORX_CMD = [
    "python3", "interprorx.py",
    "-i", str(INTERPRO_TSV),
    "-f", str(COMBINED_FASTA),
    "--min-cluster-size", "1",
    "--add-annotation",
    "--outdir", str(OUTDIR)
]

# Paths to cluster analysis outputs
CLUSTERS_DIR = OUTDIR / "clusters"
SUMMARY_FILE = OUTDIR / "clusters_summary.tsv"
NOVEL_FILE = OUTDIR / "novel_clusters.tsv"

# -----------------------------
# STEP 1: Run interprorx.py
# -----------------------------
print("🚀 Running interprorx.py...")
try:
    subprocess.run(INTERPRORX_CMD, check=True)
    print("✅ interprorx.py completed successfully.")
except subprocess.CalledProcessError as e:
    print("❌ interprorx.py failed:")
    print(e)
    exit(1)

# -----------------------------
# STEP 2: Parse cluster files
# -----------------------------
print("📂 Reading cluster files from:", CLUSTERS_DIR)
clusters = []

for cluster_file in sorted(CLUSTERS_DIR.glob("cluster_*.txt")):
    cluster_id = cluster_file.stem  # e.g. cluster_1
    with open(cluster_file) as f:
        seq_ids = [line.strip() for line in f if line.strip()]

    if not seq_ids:
        continue

    # Determine sources
    sources = set()
    for sid in seq_ids:
        if sid.startswith("DATABASE|"):
            sources.add("DATABASE")
        elif sid.startswith("PRED|"):
            sources.add("PRED")
        else:
            sources.add("UNKNOWN")

    clusters.append({
        "cluster_id": cluster_id,
        "cluster_size": len(seq_ids),
        "sources": ";".join(sorted(sources)),
        "members": ";".join(seq_ids),
        "has_DATABASE": "DATABASE" in sources,
        "has_PRED": "PRED" in sources,
    })

if not clusters:
    print("⚠️ No cluster files found in:", CLUSTERS_DIR)
    exit(1)

# -----------------------------
# STEP 3: Generate summary
# -----------------------------
df = pd.DataFrame(clusters)
df.to_csv(SUMMARY_FILE, sep="\t", index=False)
print(f"✅ Cluster summary written to: {SUMMARY_FILE}")

# -----------------------------
# STEP 4: Identify novel clusters
# -----------------------------
novel_clusters = df[(df["has_PRED"]) & (~df["has_DATABASE"])]
novel_clusters = novel_clusters.sort_values("cluster_size", ascending=False)
novel_clusters.to_csv(NOVEL_FILE, sep="\t", index=False)

print(f"🌱 Novel clusters found: {len(novel_clusters)}")
print(f"✅ Novel cluster list written to: {NOVEL_FILE}")

print("🎉 Done! Full analysis completed successfully.")
