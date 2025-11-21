from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np
import matplotlib
matplotlib.use("Agg")

# === Input for mature candidates ===
FASTA_1 = "../data/database_amps/database_creation/database_mature_APD_CAMP_DBAASP_no_dups.fasta"
FASTA_2 = "../data/predicted_amps/mature_candidates.faa"
FASTA_3 = "../data/physico-chem_data/mature_predicted_novel.fasta"

LABEL_1 = "Database Mature AMPs"
LABEL_2 = "Predicted Mature AMPs"
LABEL_3 = "Predicted Mature Novel AMPs"

OUTPUT_FIG_PATH = "../results/physicochem_analyses/AMP_mature_physicochemical_RDA.png"
OUTPUT_TSV = "../results/physicochem_analyses/AMP_mature_physicochemical_RDA.tsv"


def compute_hydrophobicity(seq):
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


# --- Load data ---
df1 = fasta_to_dataframe(FASTA_1, LABEL_1)
df2 = fasta_to_dataframe(FASTA_2, LABEL_2)
df3 = fasta_to_dataframe(FASTA_3, LABEL_3)
df = pd.concat([df1, df2, df3], ignore_index=True)

print(f"Loaded {len(df1)} sequences from {LABEL_1}, "
      f"{len(df2)} from {LABEL_2}, and {len(df3)} from {LABEL_3}.")


# === FEATURE MATRIX ===
features = [
    "Length", "GRAVY", "Charge", "Aromaticity",
    "InstabilityIndex", "IsoelectricPoint", "Hydrophobicity"
]

X = df[features].values
X = StandardScaler().fit_transform(X)


# === GROUP MATRIX ===
Y = pd.get_dummies(df["Group"], drop_first=False).values


# =======================================================
# ===============   RDA IMPLEMENTATION   ================
# =======================================================

# 1) Calculate fitted values X_hat = Y (Y'Y)^(-1) Y' X
XtY = Y.T @ X
inv = np.linalg.inv(Y.T @ Y)
B = inv @ XtY
X_hat = Y @ B

# 2) PCA over fitted values
pca = PCA(n_components=2)
rda_scores = pca.fit_transform(X_hat)

df["RDA1"] = rda_scores[:, 0]
df["RDA2"] = rda_scores[:, 1]

explained_var = pca.explained_variance_ratio_


# =======================================================
# ==========   PERMUTATION TEST (ANOVA.CCA)   ===========
# =======================================================

def permutation_test(X, Y, n_perm=999, random_state=42):
    rng = np.random.default_rng(random_state)

    # Original inertia
    pca_real = PCA(n_components=X.shape[1]).fit(Y @ np.linalg.inv(Y.T @ Y) @ Y.T @ X)
    inertia_real = sum(pca_real.explained_variance_)

    count = 0
    for _ in range(n_perm):
        Y_perm = rng.permutation(Y)
        Bp = np.linalg.inv(Y_perm.T @ Y_perm) @ (Y_perm.T @ X)
        X_hp = Y_perm @ Bp

        pca_p = PCA(n_components=X.shape[1]).fit(X_hp)
        inertia_perm = sum(pca_p.explained_variance_)

        if inertia_perm >= inertia_real:
            count += 1

    p_value = (count + 1) / (n_perm + 1)
    return p_value, inertia_real


p_value, inertia_real = permutation_test(X, Y, n_perm=999)
print(f"Permutation ANOVA-like test p-value: {p_value:.4f}")


# =======================================================
# ==============   BIPLOT LOADINGS (ARROWS)  ============
# =======================================================

loadings = pd.DataFrame(
    pca.components_.T,
    index=features,
    columns=["RDA1_loading", "RDA2_loading"]
)


# =======================================================
# ===============   RDA PLOT (POINTS + ARROWS) ==========
# =======================================================

plt.figure(figsize=(8, 6))

# Colors for groups
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

for label, color in zip(df["Group"].unique(), colors):
    subset = df[df["Group"] == label]
    plt.scatter(subset["RDA1"], subset["RDA2"], label=label, alpha=0.7, s=60, c=color)

# Draw arrows
for feat in features:
    x, y = loadings.loc[feat, "RDA1_loading"], loadings.loc[feat, "RDA2_loading"]
    plt.arrow(0, 0, x, y, head_width=0.03, color="black")
    plt.text(x * 1.1, y * 1.1, feat, fontsize=9)

plt.xlabel(f"RDA1 ({explained_var[0]*100:.1f}% variance)")
plt.ylabel(f"RDA2 ({explained_var[1]*100:.1f}% variance)")
plt.title(f"RDA (Permutation p = {p_value:.4f})")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_FIG_PATH, dpi=300)


# =======================================================
# =================   SAVE OUTPUT FILES  =================
# =======================================================

# Save dataframe with RDA scores
df.to_csv(OUTPUT_TSV, sep="\t", index=False)

# Save loadings
loadings.to_csv(OUTPUT_TSV.replace(".tsv", "_variable_loadings.tsv"),
                sep="\t")

print(f"Saved RDA scores to {OUTPUT_TSV}")
print(f"Saved variable loadings to {OUTPUT_TSV.replace('.tsv', '_variable_loadings.tsv')}")
print(f"Saved RDA plot to {OUTPUT_FIG_PATH}")
