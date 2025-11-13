#!/usr/bin/env python3
"""
prepare_go_table.py
Extract, clean, and consolidate GO annotations from InterProScan results.
Keeps only columns 1 (sequence ID) and 14 (GO terms), removes duplicates,
and merges multiple GO annotations per sequence into unique sets.
"""

import pandas as pd
from pathlib import Path

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_FILE = Path("../results/interpro_results/without_domains_predictions_and_database_families_interpro.tsv")  # raw InterProScan TSV
OUTPUT_FILE = Path("../results/go_analysis/interpro_GO_unique.tsv")

# -----------------------------
# LOAD DATA
# -----------------------------
print("🧬 Loading InterProScan results...")
# InterProScan TSV files usually have no header, column 0 = sequence ID, column 13 = GO terms (0-indexed)
df = pd.read_csv(INPUT_FILE, sep="\t", header=None, dtype=str)

# Select only columns 1 (ID) and 14 (GO terms)
df = df.iloc[:, [0, 13]]
df.columns = ["seq_id", "go_terms"]

# Remove completely duplicated rows
df.drop_duplicates(inplace=True)

# Filter out rows with no GO term or with "-"
df = df[df["go_terms"].notna() & (df["go_terms"] != "-")]

# -----------------------------
# MERGE UNIQUE GO TERMS PER SEQUENCE
# -----------------------------
def merge_unique_gos(go_list):
    """Merge and deduplicate GO terms, stripping annotation sources (e.g., '(InterPro)')."""
    terms = set()
    for entry in go_list:
        for term in entry.split("|"):
            term = term.strip()
            if term.startswith("GO:"):
                term = term.split("(")[0]  # remove "(InterPro)" or similar
                terms.add(term)
    return "|".join(sorted(terms))

df_grouped = (
    df.groupby("seq_id")["go_terms"]
    .apply(merge_unique_gos)
    .reset_index()
)

# -----------------------------
# SAVE RESULTS
# -----------------------------
OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
df_grouped.to_csv(OUTPUT_FILE, sep="\t", index=False)

print(f"✅ Clean GO table created: {len(df_grouped)} unique sequences → {OUTPUT_FILE}")
