#!/usr/bin/env python3
"""
group_with_mmseqs_shortAMPs_v2.py

Clusters predicted AMP sequences with a unified reference database using MMseqs2,
optimized for short peptide detection (AMPs embedded in longer precursors).
"""

import subprocess
import shutil
from pathlib import Path
from datetime import datetime
import pandas as pd
from Bio import SeqIO

# -----------------------
# CONFIGURATION
# -----------------------
database_fasta = "../data/database_amps/DATABASE_AMPs_clean.fasta"  # unified DB: APD3, DBAASP, CAMP, UniProt...
predicted_fasta = "../data/predicted_amps/predicted_precursors_with_flag.fasta"  # predicted sequences (from your model)
OUTDIR = Path("../results/mmseqs_results")

# MMseqs2 parameters (tuned for short peptides)
MIN_IDENTITY = 0.3       # minimal sequence identity
SENSITIVITY = 7.5        # max sensitivity
COV_MODE = 1             # coverage relative to shorter sequence
COVERAGE = 0.5
THREADS = 8
TMP_DIR = OUTDIR / "tmp"

# -----------------------
# PREPARATION
# -----------------------
OUTDIR.mkdir(parents=True, exist_ok=True)
TMP_DIR.mkdir(parents=True, exist_ok=True)
LOG = OUTDIR / f"log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

combined_fasta = OUTDIR / "combined_all.fasta"
if combined_fasta.exists():
    combined_fasta.unlink()

mapping = {}
seq_count = 0

def write_clean_fasta(in_path, out_path):
    """Append sequences to output FASTA ensuring IDs are clean."""
    global seq_count
    with open(out_path, "a") as out_handle:
        for rec in SeqIO.parse(in_path, "fasta"):
            if len(rec.seq) < 10:
                continue
            clean_id = rec.id.replace(",", "").replace(" ", "_")
            rec.id = clean_id
            rec.name = clean_id
            rec.description = ""
            SeqIO.write(rec, out_handle, "fasta")
            seq_count += 1

print("🧬 Merging FASTA inputs...")
write_clean_fasta(database_fasta, combined_fasta)
write_clean_fasta(predicted_fasta, combined_fasta)
print(f"✅ Combined FASTA: {seq_count} sequences → {combined_fasta}")

# -----------------------
# MMSEQS2 CLUSTERING
# -----------------------
print("🚀 Running MMseqs2 clustering...")

DB = OUTDIR / "db"
CLUST = OUTDIR / "clusters"

cmds = [
    ["mmseqs", "createdb", str(combined_fasta), str(DB)],
    [
        "mmseqs", "cluster", str(DB), str(CLUST), str(TMP_DIR),
        "--min-seq-id", str(MIN_IDENTITY),
        "-s", str(SENSITIVITY),
        "--cov-mode", str(COV_MODE),
        "-c", str(COVERAGE),
        "--alignment-mode", "3",
        "--remove-tmp-files", "1",
        "--threads", str(THREADS)
    ],
    ["mmseqs", "createtsv", str(DB), str(DB), str(CLUST), str(OUTDIR / "clusters.tsv")]
]

for cmd in cmds:
    print("+", " ".join(map(str, cmd)))
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(res.stderr)
        raise RuntimeError(f"MMseqs2 command failed: {' '.join(cmd)}")

print("✅ MMseqs2 clustering completed")

# -----------------------
# PARSE CLUSTERS
# -----------------------
print("📊 Parsing MMseqs2 results...")
df = pd.read_csv(OUTDIR / "clusters.tsv", sep="\t", header=None, names=["rep", "member"])

# Identify source from prefix
df["source_rep"] = df["rep"].apply(lambda x: "DATABASE" if x.startswith("DATABASE|") else "PRED" if x.startswith("PRED|") else "UNKNOWN")
df["source_member"] = df["member"].apply(lambda x: "DATABASE" if x.startswith("DATABASE|") else "PRED" if x.startswith("PRED|") else "UNKNOWN")

# summary per representative
summary = (
    df.groupby("rep")
    .agg(
        cluster_size=("member", "size"),
        sources_rep=("source_rep", lambda s: ";".join(sorted(set(s)))),
        sources_members=("source_member", lambda s: ";".join(sorted(set(s)))),
        members=("member", lambda s: ";".join(s))
    )
    .reset_index()
)

# Detect presence/absence
summary["has_DATABASE"] = summary.apply(lambda row: "DATABASE" in row["sources_rep"].split(";") or "DATABASE" in row["sources_members"].split(";"), axis=1)
summary["has_PRED"] = summary.apply(lambda row: "PRED" in row["sources_rep"].split(";") or "PRED" in row["sources_members"].split(";"), axis=1)

summary.to_csv(OUTDIR / "clusters_summary.tsv", sep="\t", index=False)
print(f"✅ Cluster summary → {OUTDIR / 'clusters_summary.tsv'}")

# -----------------------
# NOVEL CLUSTERS (PRED only)
# -----------------------
novel_clusters = summary[(summary["has_PRED"]) & (~summary["has_DATABASE"])]
novel_clusters.sort_values("cluster_size", ascending=False, inplace=True)
novel_clusters.to_csv(OUTDIR / "novel_clusters.tsv", sep="\t", index=False)
print(f"🌱 Novel clusters found: {len(novel_clusters)}")

print("🎉 Done! Results in:", OUTDIR)
