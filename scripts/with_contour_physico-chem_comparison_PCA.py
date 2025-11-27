from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np
import matplotlib
matplotlib.use("Agg")   # Non-interactive backend for WSL or servers



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
        "Charge": analysed.charge_at_pH(7.0),
        "Aromaticity": analysed.aromaticity(),
        "InstabilityIndex": analysed.instability_index(),
        "IsoelectricPoint": analysed.isoelectric_point(),
        "Hydrophobicity": compute_hydrophobicity(seq)
    }

# === Remove extreme outliers in PC1 and PC2 using IQR ===

def remove_outliers_iqr(df, cols, k=1.5):
    """Remove rows where any of the specified cols are outside the IQR range."""
    mask = pd.Series([True] * len(df))
    
    for col in cols:
        Q1 = df[col].quantile(0.25)
        Q3 = df[col].quantile(0.75)
        IQR = Q3 - Q1
        
        lower_bound = Q1 - k * IQR
        upper_bound = Q3 + k * IQR
        
        mask &= df[col].between(lower_bound, upper_bound)
    
    removed = len(df) - mask.sum()
    print(f"Removed {removed} outliers using IQR.")
    
    return df[mask]


def get_unique_records(input_fasta):
    """
    Read a FASTA file and return only the unique SeqRecords.
    Duplicates are detected by identical sequence content.
    """

    seen = set()
    unique_records = []

    for rec in SeqIO.parse(input_fasta, "fasta"):
        seq = str(rec.seq)
        if seq not in seen:
            seen.add(seq)
            unique_records.append(rec)

    return unique_records

def records_to_dataframe(records, label, compute_physicochemical_props):
    """
    Convert a list of SeqRecords into a DataFrame of physicochemical properties.
    """

    rows = []

    for rec in records:
        props = compute_physicochemical_props(rec)
        if props:
            props["Group"] = label
            rows.append(props)

    return pd.DataFrame(rows)


import numpy as np
from matplotlib.patches import Ellipse

def plot_cov_ellipse(mean, cov, ax, n_std=2.0, edgecolor="black", linewidth=2, alpha=0.25):
    """
    Plots an ellipse representing the covariance matrix centered at mean.
    n_std: how many standard deviations the ellipse should extend (2 = 95%)
    """
    # Eigenvalues and eigenvectors
    vals, vecs = np.linalg.eigh(cov)

    # Sort by largest eigenvalue
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]

    # Angle of ellipse in degrees
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    # Width and height of ellipse (2*std)
    width, height = 2 * n_std * np.sqrt(vals)

    # Create ellipse
    ellipse = Ellipse(
        xy=mean,
        width=width,
        height=height,
        angle=theta,
        edgecolor=edgecolor,
        facecolor=edgecolor,
        alpha=alpha,
        linewidth=linewidth
    )
    ax.add_patch(ellipse)

def main(FASTA_1, FASTA_2, FASTA_3, LABEL_1, LABEL_2, LABEL_3, analysis):

   # --- Remove duplicates ---
    # --- Compute properties ---
    # --- Generate final dataframe ---
    unique_seqs_1 = get_unique_records(FASTA_1)
    unique_seqs_2 = get_unique_records(FASTA_2)
    unique_seqs_3 = get_unique_records(FASTA_3)

    df1 = records_to_dataframe(unique_seqs_1, label=LABEL_1, compute_physicochemical_props=compute_physicochemical_props)
    df2 = records_to_dataframe(unique_seqs_2, label=LABEL_2, compute_physicochemical_props=compute_physicochemical_props)
    df3 = records_to_dataframe(unique_seqs_3, label=LABEL_3, compute_physicochemical_props=compute_physicochemical_props)

    df = pd.concat([df1, df2, df3], ignore_index=True)

    print(f"Loaded {len(df1)} sequences from {LABEL_1}, "
          f"{len(df2)} from {LABEL_2}, and {len(df3)} from {LABEL_3}.")

    # --- PCA ---
    features = [
        "Length", "Charge", "Aromaticity",
        "InstabilityIndex", "IsoelectricPoint", "Hydrophobicity"
    ]
    X = df[features].values
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2)
    components = pca.fit_transform(X_scaled)
    df["PC1"] = components[:, 0]
    df["PC2"] = components[:, 1]

    # Apply outlier filtering
    df = remove_outliers_iqr(df, ["PC1", "PC2"], k=6)

    explained_var = pca.explained_variance_ratio_

    # --- Plot ---
    plt.figure(figsize=(8, 6))

    # --- Decide colours for each plot ---
    if analysis == "mature":
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    elif analysis == "precursors":
        colors = ["#ff7f0e", "#1f77b4", "#2ca02c"]
    elif analysis == "dark_proteome":
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

    # ======================================================================
    #   ADD KDE CONTOUR FOR EACH GROUP
    # ======================================================================
    from scipy.stats import gaussian_kde
    import numpy as np

    for label, color in zip(df["Group"].unique(), colors):

        subset = df[df["Group"] == label]

        x = subset["PC1"].values
        y = subset["PC2"].values

        # Scatter points
        plt.scatter(x, y, label=label, alpha=0.0, s=55, c=color, edgecolor="black", linewidth=0.4)

        # KDE only if enough points
        if len(subset) > 10:
            xmin, xmax = x.min(), x.max()
            ymin, ymax = y.min(), y.max()
            Xgrid, Ygrid = np.meshgrid(
                np.linspace(xmin, xmax, 200),
                np.linspace(ymin, ymax, 200)
            )
            grid = np.vstack([Xgrid.ravel(), Ygrid.ravel()])

            kde = gaussian_kde(np.vstack([x, y]))
            Z = kde(grid).reshape(Xgrid.shape)

            # --- NEW: Slight fill under contour ---
            plt.contourf(
                Xgrid, Ygrid, Z,
                levels=[np.percentile(Z, 80), Z.max()],  # between contour and max density
                colors=[color],
                alpha=0.15    # <-- soft transparency
            )

            # Draw only one clean contour per group — the outer shape
            plt.contour(
                Xgrid, Ygrid, Z,
                levels=[np.percentile(Z, 80)],  # adjust smoothness of contour
                colors=[color],
                linewidths=2
            )

    # --- Handmade legend for the ellipses
    from matplotlib.patches import Patch

    handles = []
    for label, color in zip(df["Group"].unique(), colors):
        handles.append(Patch(edgecolor=color, facecolor='none', linewidth=2, label=label))

    plt.legend(handles=handles, loc="best")

    plt.xlabel(f"PC1 ({explained_var[0]*100:.1f}% variance)")
    plt.ylabel(f"PC2 ({explained_var[1]*100:.1f}% variance)")
    #plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_PATH, dpi=300)



"""
    # === Plot loadings as horizontal barplots (positive vs negative) ===

    # Compute loadings dataframe
    loadings = pd.DataFrame(
        pca.components_.T,
        columns=["PC1", "PC2"],
        index=features
    )

    # Plot loadings
    for pc in ["PC1", "PC2"]:
        plt.figure(figsize=(7, 5))
        loadings[pc].sort_values().plot(
            kind="barh",
            color=["#d62728" if v < 0 else "#1f77b4" for v in loadings[pc].sort_values()]
        )
        plt.axvline(0, color="black", linewidth=1)
        plt.xlabel("Contribution (Loading)")
        plt.title(f"Variable Loadings for {pc}")
        plt.tight_layout()
        plt.savefig(
            OUTPUT_FIG_PATH.replace(".png", f"_{pc}_loadings.png"),
            dpi=300
        )
        plt.close()

    print("Saved loading barplots for PC1 and PC2.")
    """



if __name__ == "__main__":

    # === Input for mature candidates ===
    FASTA_1 = "../data/database_amps/database_creation/database_mature_APD_CAMP_DBAASP_no_dups.fasta"
    FASTA_2 = "../data/predicted_amps/mature_candidates.faa"
    FASTA_3 = "../data/physico-chem_data/mature_predicted_novel.fasta"

    LABEL_1 = "Database Mature AMPs"
    LABEL_2 = "Predicted Mature AMPs"
    LABEL_3 = "Predicted Mature Novel AMPs"

    OUTPUT_FIG_PATH = "../results/physicochem_analyses/AMP_mature_physicochemical_PCA_with_contour.png"
    OUTPUT_TSV = "../results/physicochem_analyses/AMP_mature_physicochemical_properties.tsv"

    main(FASTA_1, FASTA_2, FASTA_3, LABEL_1, LABEL_2, LABEL_3, analysis="mature")
    print("#################################")


    # === Input for predicted sequences (consensus) not all===

    FASTA_1 = "../data/predicted_amps/predicted_precursors_with_flag.fasta"
    FASTA_2 = "../data/database_amps/database_creation/database_precursors_training_plus_uniprot.fasta"
    FASTA_3 = "../data/physico-chem_data/shared_novel_seqs.fasta"

    LABEL_1 = "Predicted AMP Precursors"
    LABEL_2 = "Database AMP Precursors"
    LABEL_3 = "Predicted Novel AMP Precursors"

    OUTPUT_FIG_PATH = "../results/physicochem_analyses/AMP_precursors_physicochemical_PCA_with_contour.png"
    OUTPUT_TSV = "../results/physicochem_analyses/AMP_precursors_physicochemical_properties.tsv"

    main(FASTA_1, FASTA_2, FASTA_3, LABEL_1, LABEL_2, LABEL_3, analysis="precursors")
    print("#################################")

    # === Input for uncharacterized and hypothetical predicted (consensus not all) and novel sequences.

    FASTA_1 = "../data/database_amps/database_creation/database_precursors_training_plus_uniprot.fasta"
    FASTA_2 = "../data/physico-chem_data/uncharacterized_shared_novel_seqs.fasta"
    FASTA_3 = "../data/physico-chem_data/hypothetical_shared_novel_seqs.fasta"

    LABEL_1 = "Database AMP Precursors"
    LABEL_2 = "Predicted AMP Precursors Uncharacterized"
    LABEL_3 = "Predicted AMP Precursors Hypothetical"

    OUTPUT_FIG_PATH = "../results/physicochem_analyses/uncharacterized_AMP_precursors_physicochemical_PCA_with_contour.png"
    OUTPUT_TSV = "../results/physicochem_analyses/uncharacterized_AMP_precursors_physicochemical_properties.tsv"

    main(FASTA_1, FASTA_2, FASTA_3, LABEL_1, LABEL_2, LABEL_3, analysis="dark_proteome")
    print("#################################")






    # === Input for mature candidates (mmseqs & psi) ===
    FASTA_1 = "../data/database_amps/database_creation/database_mature_APD_CAMP_DBAASP_no_dups.fasta"
    FASTA_2 = "../data/predicted_amps/mature_candidates.faa"
    FASTA_3 = "../data/physico-chem_data/plus_psi_mmseqs_mature_novel_seqs.fasta"

    LABEL_1 = "Database Mature AMPs"
    LABEL_2 = "Predicted Mature AMPs"
    LABEL_3 = "Predicted Mature Novel AMPs (psi & mmseqs)"

    OUTPUT_FIG_PATH = "../results/physicochem_analyses/psi_mmseqs_AMP_mature_physicochemical_PCA_with_contour.png"
    OUTPUT_TSV = "../results/physicochem_analyses/psi_mmseqs_AMP_mature_physicochemical_properties.tsv"

    main(FASTA_1, FASTA_2, FASTA_3, LABEL_1, LABEL_2, LABEL_3, analysis="mature")
    print("#################################")


     # === Input for predicted sequences (consensus) not all===

    FASTA_1 = "../data/predicted_amps/predicted_precursors_with_flag.fasta"
    FASTA_2 = "../data/database_amps/database_creation/database_precursors_training_plus_uniprot.fasta"
    FASTA_3 = "../data/physico-chem_data/plus_psi_mmseqs_novel_seqs.fasta"

    LABEL_1 = "Predicted AMP Precursors"
    LABEL_2 = "Database AMP Precursors"
    LABEL_3 = "Predicted Novel AMP Precursors (psi & mmseqs)"

    OUTPUT_FIG_PATH = "../results/physicochem_analyses/psi_mmseqs_AMP_precursors_physicochemical_PCA_with_contour.png"
    OUTPUT_TSV = "../results/physicochem_analyses/psi_mmseqs_AMP_precursors_physicochemical_properties.tsv"

    main(FASTA_1, FASTA_2, FASTA_3, LABEL_1, LABEL_2, LABEL_3, analysis="precursors")
    print("#################################")


