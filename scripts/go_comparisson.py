#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# Ajusta estas rutas
GO_FILE = Path("../data/go_database/interpro_GO_unique_predicted_plus_database.tsv")
OUTPUT = Path("../results/go_analysis/go_enrichment_pred_vs_db_fisher.tsv")

# Cargar
df = pd.read_csv(GO_FILE, sep="\t", dtype=str)

# Asegúrate que columnas existen: seq_id, go_terms
if "seq_id" not in df.columns or "go_terms" not in df.columns:
    raise ValueError("El fichero debe contener 'seq_id' y 'go_terms'.")

# Normalizar go_terms: lista por fila
df["go_list"] = df["go_terms"].fillna("").map(lambda s: [t for t in s.split("|") if t.startswith("GO:")])

# Separar sets
pred_df = df[df["seq_id"].str.startswith("PRED|")].copy()
db_df   = df[df["seq_id"].str.startswith("DATABASE|")].copy()

N_pred = len(pred_df)
N_db   = len(db_df)
print(f"Predicted sequences: {N_pred}; Database sequences: {N_db}")

# Construir conteos por GO
from collections import Counter, defaultdict

pred_counter = Counter()
db_counter = Counter()
all_gos = set()

for lst in pred_df["go_list"]:
    pred_counter.update(lst)
    all_gos.update(lst)

for lst in db_df["go_list"]:
    db_counter.update(lst)
    all_gos.update(lst)

# Para cada GO calcular tabla 2x2: a,b,c,d
rows = []
for go in sorted(all_gos):
    a = pred_counter.get(go, 0)          # predicted with GO
    c = db_counter.get(go, 0)            # database with GO
    b = N_pred - a                       # predicted without GO
    d = N_db - c                         # database without GO

    # Seguridad: si alguno de los totales es negativo algo anda mal
    if b < 0 or d < 0:
        print("Warning negative counts for GO", go, "a,c,b,d:", a,c,b,d)
        continue

    # Fisher exact test (alternative='greater' para ver sobre-representación en predicted)
    # Use two-sided if you prefer to detect depletion also.
    # For enrichment (predicted > database) we use alternative='greater'
    try:
        oddsr, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
    except Exception as e:
        oddsr, pvalue = (np.nan, 1.0)

    rows.append({
        "GO": go,
        "pred_count": a,
        "db_count": c,
        "N_pred": N_pred,
        "N_db": N_db,
        "oddsratio": oddsr,
        "p_uncorrected": pvalue
    })

res_df = pd.DataFrame(rows)
# Ajuste BH
res_df["p_fdr_bh"] = multipletests(res_df["p_uncorrected"], method="fdr_bh")[1]

# Añadir nombre del GO si lo tienes en el dataframe original mapeado
# Si tu archivo contiene mapeo GO -> name puedes leerlo; si no, se puede extraer
# desde un fichero separado. Aquí intento usar df original como fuente:
go_name_map = {}
# construir mapeo GO -> any description si está presente en alguna fila
if "go_desc" in df.columns:
    # si tu tabla contiene descripciones por GO en cada fila
    for _, row in df.iterrows():
        for g in row["go_list"]:
            if g not in go_name_map and "go_desc" in row and pd.notna(row["go_desc"]):
                go_name_map[g] = row["go_desc"]

res_df["Name"] = res_df["GO"].map(go_name_map).fillna("NA")

# Ordenar por FDR y por odds ratio
res_df = res_df.sort_values(["p_fdr_bh", "p_uncorrected", "oddsratio"])

# Guardar
res_df.to_csv(OUTPUT, sep="\t", index=False, float_format="%.6g")
print(f"Saved results to {OUTPUT}")

# DIAGNOSTICOS rápidos: casos con pred_count == (pred_count + db_count) i.e. db_count==0
issue = res_df[res_df["db_count"] == 0].copy()
if not issue.empty:
    print("\n--- Diagnóstico: GO con db_count == 0 (solo aparecen en predicted) ---")
    print(issue[["GO", "pred_count", "db_count", "p_uncorrected", "p_fdr_bh"]].head(20).to_string(index=False))
    print("\nNota: si hay muchos de estos, indica que muchos GO solo están en tus predichos y no en la base de datos.")
else:
    print("\nNo GO exclusivos de predicted detectados (db_count > 0 para todos).")
