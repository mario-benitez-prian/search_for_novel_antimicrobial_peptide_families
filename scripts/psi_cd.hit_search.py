#!/usr/bin/env python3
"""
group_with_psicdhit_database_pred.py

Clusters predicted precursors (PRED) and known AMPs (DATABASE)
using psi-cd-hit, identifies novel clusters composed only of
predicted sequences, and outputs results as TSV tables.
"""

import subprocess
import shutil
from pathlib import Path
from datetime import datetime
import pandas as pd
import sys
from Bio import SeqIO
import re

# ==========================================================
# CONFIGURATION (edit these variables)
# ==========================================================
database_fasta = "../data/database_amps/DATABASE_AMPs_clean_without_database.fasta"       # unified DB: APD3, DBAASP, CAMP, UniProt...
predicted_fasta = "../data/predicted_amps/predicted_precursors.fasta" # predicted sequences (from your model)
psi_cdhit_exe = "/home/mario/software/cd-hit-v4.8.1-2019-0228/psi-cd-hit/psi-cd-hit.pl"               # must be in PATH or provide full path
outdir = Path("../results/psicdhit_results")
combined_prefix = "combined_all"
identity = 0.30       # sequence identity threshold
coverage = 0.7        # coverage (-aS)
threads = 8
dry_run = False

# ==========================================================
# Helper functions
# ==========================================================
def find_executable(name):
    return shutil.which(name)

def run_cmd(cmd, dry_run=False):
    print("+", " ".join(cmd))
    if dry_run:
        return None
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print("ERROR running command:", " ".join(cmd), file=sys.stderr)
        print("stdout:", res.stdout, file=sys.stderr)
        print("stderr:", res.stderr, file=sys.stderr)
        raise RuntimeError(f"Command failed (code {res.returncode})")
    return res

def parse_cdhit_clstr(clstr_path):
    """
    Parse cd-hit/psi-cd-hit .clstr file using regex.
    Returns dict: cluster_id -> list of (seq_id, is_rep)
    """
    clusters = {}
    cur = None
    pat = re.compile(r'>(.+?)\.\.\.')
    with open(clstr_path, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>Cluster'):
                cur = line.split()[1]
                clusters[cur] = []
            else:
                m = pat.search(line)
                if m:
                    sid = m.group(1)
                    is_rep = line.strip().endswith('*')
                    clusters[cur].append((sid, is_rep))
    return clusters

# ==========================================================
# MAIN PIPELINE
# ==========================================================
outdir.mkdir(parents=True, exist_ok=True)
logf = outdir / f"run_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

with open(logf, "w") as L:
    L.write("PSI-CD-HIT pipeline run\n")
    L.write(f"Database: {database_fasta}\n")
    L.write(f"Predicted: {predicted_fasta}\n")

psi_cd_hit = find_executable(psi_cdhit_exe)
if not psi_cd_hit:
    raise FileNotFoundError(f"psi-cd-hit.pl binary '{psi_cdhit_exe}' not found in PATH.")

# ----------------------------------------------------------
# Combine FASTA inputs
# ----------------------------------------------------------
combined_fasta = outdir / f"{combined_prefix}.fasta"
if combined_fasta.exists():
    combined_fasta.unlink()

mapping = {}
seq_count = 0

def write_fasta_with_prefix(in_path, prefix, out_path):
    """Write FASTA sequences adding prefix and cleaning IDs (removes commas)."""
    global seq_count
    with open(out_path, "a") as out_handle:
        for rec in SeqIO.parse(in_path, "fasta"):
            clean_id = rec.id.replace(",", "").replace(" ", "_")
            sid = f"{prefix}|{clean_id}"
            rec.id = sid
            rec.name = sid
            rec.description = ""
            mapping[sid] = prefix
            SeqIO.write(rec, out_handle, "fasta")
            seq_count += 1

print("Merging FASTA inputs...")
write_fasta_with_prefix(database_fasta, "DATABASE", combined_fasta)
write_fasta_with_prefix(predicted_fasta, "PRED", combined_fasta)
print(f"Wrote combined FASTA with {seq_count} sequences to {combined_fasta}")

# Quick stats
lengths = [len(rec.seq) for rec in SeqIO.parse(combined_fasta, "fasta")]
short_count = sum(1 for l in lengths if l < 10)
print(f"{short_count} sequences shorter than 10 aa (warning: psi-cd-hit ignores these)")

# ----------------------------------------------------------
# Run PSI-CD-HIT
# ----------------------------------------------------------
out_prefix = outdir / f"psi_c{int(identity*100)}"
cmd = [
    psi_cdhit_exe,
    "-i", str(combined_fasta),
    "-o", str(out_prefix),
    "-c", str(identity),
    "-aS", str(coverage),
    "-g", "1"
]
run_cmd(cmd, dry_run=dry_run)
final_clstr = str(out_prefix) + ".clstr"

# ----------------------------------------------------------
# Parse clusters and generate tables
# ----------------------------------------------------------
print("Parsing final cluster file:", final_clstr)
clusters = parse_cdhit_clstr(final_clstr)

rows = []
for clid, members in clusters.items():
    rep = next((sid for sid, is_rep in members if is_rep), None)
    for sid, is_rep in members:
        src = mapping.get(sid, "UNKNOWN")
        rows.append({
            "cluster_id": clid,
            "seq_id": sid,
            "is_rep": bool(is_rep),
            "rep": rep,
            "source": src
        })

df_map = pd.DataFrame(rows)
map_tsv = outdir / "all_sequences_cluster_map.tsv"
df_map.to_csv(map_tsv, sep="\t", index=False)
print("Wrote cluster map to", map_tsv)

# Summary per cluster
def join_unique(x): return ";".join(sorted(set(x)))
summary = df_map.groupby("cluster_id").agg(
    cluster_size=("seq_id", "size"),
    rep=("rep", lambda s: s.mode().iat[0] if not s.mode().empty else ""),
    sources=("source", join_unique),
    members=("seq_id", lambda s: ";".join(s))
).reset_index()

summary["has_DATABASE"] = summary["sources"].apply(lambda s: "DATABASE" in s.split(";"))
summary["has_PRED"] = summary["sources"].apply(lambda s: "PRED" in s.split(";"))
summary_tsv = outdir / "clusters_summary.tsv"
summary.to_csv(summary_tsv, sep="\t", index=False)
print("Wrote clusters summary to", summary_tsv)

# ----------------------------------------------------------
# Novel clusters = PRED-only
# ----------------------------------------------------------
novel_clusters = summary[(summary["has_PRED"]) & (~summary["has_DATABASE"])].copy()
novel_clusters.sort_values("cluster_size", ascending=False, inplace=True)
novel_clusters_tsv = outdir / "novel_clusters.tsv"
novel_clusters.to_csv(novel_clusters_tsv, sep="\t", index=False)
print(f"Wrote {len(novel_clusters)} novel clusters to {novel_clusters_tsv}")

# Expand predicted sequences
novel_rows = []
for _, r in novel_clusters.iterrows():
    members = r["members"].split(";") if r["members"] else []
    for m in members:
        if mapping.get(m) == "PRED":
            novel_rows.append({
                "cluster_id": r["cluster_id"],
                "cluster_size": r["cluster_size"],
                "rep": r["rep"],
                "pred_seq_id": m
            })
novel_preds_df = pd.DataFrame(novel_rows)
novel_preds_tsv = outdir / "novel_predictions.tsv"
novel_preds_df.to_csv(novel_preds_tsv, sep="\t", index=False)
print(f"Wrote {len(novel_preds_df)} novel predicted sequences to {novel_preds_tsv}")

print("DONE. Results in:", outdir)
