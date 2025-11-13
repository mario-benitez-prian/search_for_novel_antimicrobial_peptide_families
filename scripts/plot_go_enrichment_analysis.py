import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")   # Non-interactive backend for WSL or servers
import numpy as np
from pathlib import Path

# === Visualization section ===
print("📊 Generating GO enrichment plots...")

# Load GO enrichment results
results_path = Path("../results/go_analysis/go_enrichment_predicted.tsv")
df = pd.read_csv(results_path, sep="\t")

# Filter enriched terms only (not depleted)
df_enriched = df[df["enrichment"].isin(["enriched"])].copy()

# Add -log10(FDR) column
df_enriched["-log10(FDR)"] = -np.log10(df_enriched["p_fdr_bh"])

# Define namespaces
namespaces = ["BP", "CC", "MF"]
colors = {
    "BP": "#6baed6",
    "CC": "#9ecae1",
    "MF": "#c6dbef"
}

# --- Generate barplots for each namespace ---
for ns in namespaces:
    df_ns = df_enriched[df_enriched["NS"] == ns]
    if df_ns.empty:
        print(f"⚠️ No enriched terms found for namespace {ns}, skipping.")
        continue

    top20_ns = df_ns.nsmallest(20, "p_fdr_bh")  # Top 20 by FDR

    plt.figure(figsize=(8, 6))
    plt.barh(top20_ns["Name"], top20_ns["-log10(FDR)"], color=colors[ns])
    plt.xlabel("-log10(FDR)")
    plt.ylabel("GO term")
    plt.title(f"Top 20 Enriched GO Terms ({ns})")
    plt.gca().invert_yaxis()
    plt.tight_layout()

    # Save plot
    barplot_path = Path(f"../results/go_analysis/go_enrichment_top20_{ns}_predicted_vs_background.png")
    plt.savefig(barplot_path, dpi=300)
    plt.close()
    print(f"✅ Saved barplot for {ns} at {barplot_path}")

# --- Optional: Pie chart of GO namespaces ---
namespace_counts = df_enriched["NS"].value_counts()
plt.figure(figsize=(5, 5))
plt.pie(namespace_counts, labels=namespace_counts.index,
        autopct="%1.1f%%", startangle=90,
        colors=[colors.get(ns, "#cccccc") for ns in namespace_counts.index])
plt.title("Distribution of Enriched GO Terms by Namespace")
plt.tight_layout()

# Save pie chart
pie_path = Path("../results/go_analysis/go_enrichment_namespace_pie_predicted_vs_background.png")
plt.savefig(pie_path, dpi=300)
plt.close()
print(f"✅ Pie chart saved at {pie_path}")
