import pandas as pd
import re

# === Configuration ===
input_file = "CAMP_2025-11-03_17-07-12.txt"
output_file = "CAMP_2025-11-03_17-07-12.fasta"

# === Read input safely ===
df = pd.read_csv(input_file, sep='\t', dtype=str, encoding='latin-1', on_bad_lines='skip')

# Clean column names
df.columns = [c.strip() for c in df.columns]

# Replace non-breaking spaces and strip whitespace from all cells
df = df.applymap(lambda x: x.replace('\xa0', ' ').strip() if isinstance(x, str) else x)

# Combine Camp_ID and Title into ID
df["ID"] = (df["Camp_ID"].astype(str).str.strip() + "_" + df["Title"].astype(str).str.strip())

# Clean ID: remove spaces and strange chars, keep only letters, numbers, underscores, and hyphens
df["ID"] = df["ID"].apply(lambda x: re.sub(r'[^A-Za-z0-9_\-]', '', x))

# Clean sequence
df["Sequence"] = df["Seqence"].astype(str).str.strip()
df["Sequence"] = df["Sequence"].str.replace(r'[^A-Za-z]', '', regex=True)

# Remove rows with invalid, empty, or NaN sequences or IDs
df = df[
    df["Sequence"].notna() &
    df["ID"].notna() &
    df["Sequence"].str.len().gt(0) &
    df["ID"].str.len().gt(0) &
    ~df["ID"].str.lower().eq("nan") &
    ~df["Sequence"].str.lower().eq("nan")
]

# Drop duplicates if any
df = df.drop_duplicates(subset=["ID", "Sequence"])

# === Write FASTA ===
with open(output_file, "w") as f:
    for _, row in df.iterrows():
        f.write(f">{row['ID']}\n{row['Sequence']}\n")

print(f"✅ Clean FASTA file saved as: {output_file}")
print(f"✅ Total valid sequences: {len(df)}")
