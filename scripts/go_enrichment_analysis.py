#!/usr/bin/env python3
"""
go_enrichment_predicted_vs_database.py
Perform GO term enrichment analysis for predicted novel peptides (PRED)
vs known database peptides (DATABASE) using GOATOOLS.
"""

import pandas as pd
from pathlib import Path
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

# -----------------------------
# CONFIGURATION
# -----------------------------
GO_FILE = Path("../results/go_analysis/interpro_GO_unique.tsv")
GO_OBO = Path("../data/go_database/go-basic.obo")
OUTPUT_FILE = Path("../results/go_analysis/go_enrichment_predicted.tsv")

# -----------------------------
# LOAD DATA
# -----------------------------
print("📖 Loading GO annotations file (predictions + database)...")
df = pd.read_csv(GO_FILE, sep="\t")

# Split into predicted and database sequences
pred_df = df[df["seq_id"].str.startswith("PRED|")]
db_df = df[df["seq_id"].str.startswith("DATABASE|")]

print(f"✅ Predicted sequences: {len(pred_df)}")
print(f"✅ Database sequences:  {len(db_df)}")

# Build GO mappings
def make_go_mapping(df):
    return {
        row["seq_id"]: [term for term in row["go_terms"].split("|") if term.startswith("GO:")]
        for _, row in df.iterrows()
    }

all_go = make_go_mapping(df)
pred_go = make_go_mapping(pred_df)

# Background set = all sequences
pop_ids = list(all_go.keys())

# Study set = predicted peptides
study_ids = list(pred_go.keys())

# GO term to gene mapping
gene2go = {k: v for k, v in all_go.items() if v}

print(gene2go)

# -----------------------------
# RUN ENRICHMENT ANALYSIS
# -----------------------------
print("🧠 Running GO enrichment analysis...")

obodag = GODag(str(GO_OBO))

goeaobj = GOEnrichmentStudy(
    pop_ids,          # background (all)
    gene2go,          # gene2go mapping
    obodag,           # ontology
    propagate_counts=False,
    alpha=0.05,
    methods=["fdr_bh"]
)

results = goeaobj.run_study(study_ids)

# -----------------------------
# SAVE RESULTS
# -----------------------------
print("💾 Saving results...")
sig_results = [r for r in results if r.p_fdr_bh < 0.05]

with open(OUTPUT_FILE, "w") as out:
    out.write("GO_ID\tNS\tName\tp_uncorrected\tp_fdr_bh\tenrichment\tstudy_count\tpop_count\n")
    for r in sig_results:
        enrich = "enriched" if r.enrichment == "e" else "depleted"
        out.write(f"{r.GO}\t{r.NS}\t{r.name}\t{r.p_uncorrected:.3e}\t{r.p_fdr_bh:.3e}\t{enrich}\t{r.study_count}\t{r.pop_count}\n")

print(f"✅ GO enrichment complete: {len(sig_results)} significant terms saved → {OUTPUT_FILE}")
