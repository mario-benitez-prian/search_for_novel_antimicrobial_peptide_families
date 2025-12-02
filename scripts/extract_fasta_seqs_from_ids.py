# ============================
# Variable declarations
# ============================
id_list_file = "../data/database_amps/database_mature_amps_ids.txt"          # Archivo con IDs, uno por línea
input_fasta = "../data/database_amps/DATABASE_AMPs_clean.fasta"   # Archivo FASTA de entrada
output_fasta = "../data/database_amps/database_mature_amps.fasta"   # Archivo FASTA de salida

ids = set()
found_ids = set()
missing_ids = set()
records = []
current_id = None
current_seq = []

# ============================
# Load ID list
# ============================
with open(id_list_file, "r") as f:
    for line in f:
        clean = line.strip()
        if clean:
            ids.add(clean)

# ============================
# Parse FASTA and extract matches
# ============================
with open(input_fasta, "r") as f:
    for line in f:
        if line.startswith(">"):
            # Save previous record if matched
            if current_id is not None and current_id in ids:
                records.append((current_id, "".join(current_seq)))
                found_ids.add(current_id)

            # Start new record
            current_id = line[1:].strip().split()[0]   # ID before first space
            current_seq = []
        else:
            current_seq.append(line.strip())

# Handle last record
if current_id is not None and current_id in ids:
    records.append((current_id, "".join(current_seq)))
    found_ids.add(current_id)

# ============================
# Determine missing IDs
# ============================
missing_ids = ids - found_ids

# ============================
# Write output FASTA
# ============================
with open(output_fasta, "w") as f:
    for rid, seq in records:
        f.write(f">{rid}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

# ============================
# Report missing IDs
# ============================
if missing_ids:
    print("WARNING: Some IDs were not found:")
    for mid in missing_ids:
        print(" -", mid)
else:
    print("All IDs were found correctly.")

print(f"\nExtracted {len(records)} sequences.\nSaved to: {output_fasta}")