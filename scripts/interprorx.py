#!/usr/bin/env python3
"""
InterProRX
----------
CLI tool to cluster protein sequences from InterProScan TSV results into families
based on shared domain/family identifiers using NetworkX.

Author: Mario Benítez-Prián
"""

import os
import argparse
import pandas as pd
import networkx as nx
from pathlib import Path


# ======================================================
# UTILS
# ======================================================

def read_fasta(fasta_file):
    """Read a FASTA file into a dictionary {seq_id: sequence}."""
    sequences = {}
    with open(fasta_file) as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq_lines)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id:
            sequences[seq_id] = "".join(seq_lines)
    return sequences


def write_clusters(components, sequences, outdir_clusters, outdir_fastas, df, add_annotation=False):
    """Write cluster ID lists and FASTA files for each cluster, optionally with InterPro annotations."""
    os.makedirs(outdir_clusters, exist_ok=True)
    os.makedirs(outdir_fastas, exist_ok=True)

    # Build mapping Protein_ID → best available description
    def best_description(subdf):
        interpro_desc = [
            v.strip() for v in subdf["InterPro_Description"].dropna().unique()
            if v.strip() != "" and v.strip() != "-"
        ]
        db_desc = [
            v.strip() for v in subdf["DB_Description"].dropna().unique()
            if v.strip() != "" and v.strip() != "-"
        ]
        if len(interpro_desc) > 0:
            desc = "|".join(sorted(interpro_desc))
        elif len(db_desc) > 0:
            desc = "|".join(sorted(db_desc))
        else:
            desc = "Unknown"
        return desc

    grouped = df.loc[:, df.columns != "Protein_ID"].groupby(df["Protein_ID"], group_keys=False)
    desc_map = grouped.apply(best_description).to_dict()

    # Write clusters
    for i, comp in enumerate(components, start=1):
        cluster_file = Path(outdir_clusters) / f"cluster_{i}.txt"
        fasta_file = Path(outdir_fastas) / f"cluster_{i}.fasta"

        with open(cluster_file, "w") as cf, open(fasta_file, "w") as ff:
            for seq_id in comp:
                seq = sequences.get(seq_id)
                desc = desc_map.get(seq_id, "Unknown")
                desc = desc.replace(" ", "_")

                # Only add annotation if the flag is enabled
                header = f"{seq_id}||{desc}" if add_annotation else seq_id

                # Write to cluster list
                cf.write(f"{header}\n")

                # Write to FASTA
                if seq:
                    ff.write(f">{header}\n{seq}\n")
                else:
                    ff.write(f">{header}\nSequence_not_found\n")
                    print(f"Warning: sequence {seq_id} not found in FASTA.")


# ======================================================
# CORE
# ======================================================

def read_and_filter_interpro(interpro_tsv, min_eval=None, min_len=None, max_len=None):
    """Read and filter InterProScan TSV according to user criteria."""
    usecols = [0, 2, 3, 4, 5, 8, 11, 12]
    colnames = [
        "Protein_ID", "Sequence_Length", "Database", "DB_ID",
        "DB_Description", "Evalue", "InterPro_ID", "InterPro_Description"
    ]

    df = pd.read_csv(
        interpro_tsv,
        sep="\t",
        header=None,
        usecols=usecols,
        names=colnames,
        dtype=str,
        na_values=["-"]
    )

    df = df[~df["Database"].isin(["MobiDBLite", "Coils"])]
    df["Evalue"] = pd.to_numeric(df["Evalue"], errors="coerce")
    df["Sequence_Length"] = pd.to_numeric(df["Sequence_Length"], errors="coerce")

    if min_eval is not None:
        mask_eval = (df["Evalue"] <= min_eval) | (df["Database"].isin(["ProSitePatterns", "ProSiteProfiles"]))
        df = df[mask_eval]

    if min_len is not None:
        df = df[df["Sequence_Length"] >= min_len]
    if max_len is not None:
        df = df[df["Sequence_Length"] <= max_len]

    return df


def build_graph_from_interpro(df):
    """Build a NetworkX graph connecting proteins and their InterPro annotations."""
    G = nx.Graph()
    connect_cols = [
        "Protein_ID", "DB_ID", "DB_Description",
        "InterPro_ID", "InterPro_Description"
    ]
    for _, row in df.iterrows():
        nodes = [row[c] for c in connect_cols if pd.notna(row[c])]
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                G.add_edge(nodes[i], nodes[j])
    return G


def cluster_sequences(df, min_cluster_size=2):
    """Build clusters (connected components) from filtered InterPro data."""
    G = build_graph_from_interpro(df)
    components = []
    for comp in nx.connected_components(G):
        proteins = [n for n in comp if n in df["Protein_ID"].values]
        if len(proteins) >= min_cluster_size:
            components.append(sorted(proteins))
    return components


# ======================================================
# CLI
# ======================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "InterProRX: Cluster protein sequences from InterProScan results using NetworkX.\n"
            "Generates FASTA files for each protein family identified by shared domain/family IDs.\n\n"
            "Author: Mario Benítez-Prián"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--interpro", required=True, help="InterProScan TSV results file")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file with all protein sequences")
    parser.add_argument("--outdir", default="interprorx_output", help="Output directory (default=interprorx_output)")
    parser.add_argument("--min-cluster-size", type=int, default=2, help="Minimum cluster size to output (default=2)")
    parser.add_argument("--min-evalue", type=float, default=None, help="Maximum allowed e-value (default=None)")
    parser.add_argument("--min-len", type=int, default=None, help="Minimum sequence length (default=None)")
    parser.add_argument("--max-len", type=int, default=None, help="Maximum sequence length (default=None)")
    parser.add_argument("--add-annotation", action="store_true",
                        help="Add InterPro annotation to sequence headers and cluster lists")

    args = parser.parse_args()

    print("Reading and filtering InterProScan results...")
    df = read_and_filter_interpro(
        args.interpro,
        min_eval=args.min_evalue,
        min_len=args.min_len,
        max_len=args.max_len
    )
    print(f"{len(df)} rows retained after filtering.")

    print("Building clusters...")
    components = cluster_sequences(df, min_cluster_size=args.min_cluster_size)
    print(f"{len(components)} clusters found.")

    print("Reading FASTA sequences...")
    sequences = read_fasta(args.fasta)

    out_clusters = Path(args.outdir) / "clusters"
    out_fastas = Path(args.outdir) / "fastas"

    print("Writing output files...")
    write_clusters(components, sequences, out_clusters, out_fastas, df, add_annotation=args.add_annotation)

    print(f"Done! Results written to {args.outdir}/")


if __name__ == "__main__":
    main()
