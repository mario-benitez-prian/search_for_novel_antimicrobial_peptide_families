from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np
import matplotlib
matplotlib.use("Agg")   # Non-interactive backend for WSL or servers


# === Input ===
FASTA_1 = "../data/database_amps/database_creation/database_precursors_training_plus_uniprot.fasta"
FASTA_2 = "../results/uncharacterized_shared_novel_seqs.fasta"
FASTA_3 = "../results/hypothetical_shared_novel_seqs.fasta"

LABEL_1 = "Database AMP Precursors"
LABEL_2 = "Predicted AMP Precursors Uncharacterized"
LABEL_3 = "Predicted AMP Precursors Hypothetical"

OUTPUT_FIG_PATH = "../results/physicochem_analyses/uncharacterized_AMP_precursors_physicochemical_PCA.png"
OUTPUT_TSV = "../results/physicochem_analyses/uncharacterized_AMP_precursors_physicochemical_properties.tsv"


def compute_hydrophobicity(seq):
    """Compute mean Kyte–Doolittle hydrophobicity."""
    hydrophobicity_scale = {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
        "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
        "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
        "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
    }
    values = [hydrophobicity_scale[a] for a in seq if a in hydrophobicity_scale]
    return np.mean(values) if values else np.nan


def compute_physicochemical_props(seq_record):
    seq = str(seq_record.seq).upper()
    allowed_aas = set("ACDEFGHIKLMNPQRSTVWY")
    seq = "".join([aa for aa in seq if aa in allowed_aas])

    if len(seq) < 5:
        return None

    analysed = ProteinAnalysis(seq)
    return {
        "ID": seq_record.id,
        "Length": len(seq),
        "GRAVY": analysed.gravy(),
        "Charge": analysed.charge_at_pH(7.0),
        "Aromaticity": analysed.aromaticity(),
        "InstabilityIndex": analysed.instability_index(),
        "IsoelectricPoint": analysed.isoelectric_point(),
        "Hydrophobicity": compute_hydrophobicity(seq)
    }


def fasta_to_dataframe(fasta_path, label):
    records = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        props = compute_physicochemical_props(rec)
        if props:
            props["Group"] = label
            records.append(props)
    return pd.DataFrame(records)


# --- Compute properties ---
df1 = fasta_to_dataframe(FASTA_1, LABEL_1)
df2 = fasta_to_dataframe(FASTA_2, LABEL_2)
df3 = fasta_to_dataframe(FASTA_3, LABEL_3)
df = pd.concat([df1, df2, df3], ignore_index=True)

print(f"Loaded {len(df1)} sequences from {LABEL_1}, "
      f"{len(df2)} from {LABEL_2}, and {len(df3)} from {LABEL_3}.")


# --- PCA ---
features = [
    "Length", "GRAVY", "Charge", "Aromaticity",
    "InstabilityIndex", "IsoelectricPoint", "Hydrophobicity"
]
X = df[features].values
X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=2)
components = pca.fit_transform(X_scaled)
df["PC1"] = components[:, 0]
df["PC2"] = components[:, 1]

explained_var = pca.explained_variance_ratio_

# --- Plot ---
plt.figure(figsize=(8, 6))
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]  # Blue, Orange, Green

for label, color in zip(df["Group"].unique(), colors):
    subset = df[df["Group"] == label]
    plt.scatter(subset["PC1"], subset["PC2"], label=label, alpha=0.7, s=60, c=color)

plt.xlabel(f"PC1 ({explained_var[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({explained_var[1]*100:.1f}% variance)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_FIG_PATH, dpi=300)
plt.show()

# Save results
df.to_csv(OUTPUT_TSV, sep="\t", index=False)
print("Saved results to AMPs_precursors_physicochemical_properties.tsv "
      "and PCA plot to AMPs_physicochemical_PCA.png")
