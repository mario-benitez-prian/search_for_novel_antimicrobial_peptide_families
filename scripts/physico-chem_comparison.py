from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend -> avoids Qt/xcb problems


# === Input ===
FASTA_1 = "../mature_candidates.faa"  # e.g., AMPs from species A
FASTA_2 = "../amp_dbs_mature.fasta"  # e.g., AMPs from species B
LABEL_1 = "Predicted AMPs"
LABEL_2 = "Database AMPs"

def compute_physicochemical_props(seq_record):
    seq = str(seq_record.seq).upper()

    # Keep only standard amino acids (remove U, O, B, Z, J, X, etc.)
    allowed_aas = set("ACDEFGHIKLMNPQRSTVWY")
    seq = "".join([aa for aa in seq if aa in allowed_aas])

    # Skip very short or empty sequences
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
        "IsoelectricPoint": analysed.isoelectric_point()
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
df = pd.concat([df1, df2], ignore_index=True)

print(f"Loaded {len(df1)} sequences from {LABEL_1} and {len(df2)} from {LABEL_2}")

# --- PCA ---
features = ["Length", "GRAVY", "Charge", "Aromaticity", "InstabilityIndex", "IsoelectricPoint"]
X = df[features].values
X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=2)
components = pca.fit_transform(X_scaled)
df["PC1"] = components[:, 0]
df["PC2"] = components[:, 1]

explained_var = pca.explained_variance_ratio_

# --- Plot ---
plt.figure(figsize=(8, 6))
for label, color in zip(df["Group"].unique(), ["#1f77b4", "#ff7f0e"]):
    subset = df[df["Group"] == label]
    plt.scatter(subset["PC1"], subset["PC2"], label=label, alpha=0.7, s=60)

plt.xlabel(f"PC1 ({explained_var[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({explained_var[1]*100:.1f}% variance)")
plt.title("PCA of Physicochemical Properties of AMPs")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("AMPs_physicochemical_PCA.png", dpi=300)
plt.show()

# Save results
df.to_csv("AMPs_physicochemical_properties.tsv", sep="\t", index=False)
print("Saved results to AMPs_physicochemical_properties.tsv and PCA plot to AMPs_physicochemical_PCA.png")
