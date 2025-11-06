#!/usr/bin/env python3
import sys

# === Configuration ===
input_file = "DATABASE_AMPs.fasta"
output_file = "DATABASE_AMPs_clean.fasta"

sequences = {}
header = None
seq_lines = []

with open(input_file, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue  # Skip empty lines

        if line.startswith(">"):
            # Save previous sequence if exists
            if header and seq_lines:
                seq = "".join(seq_lines).upper()
                if seq not in sequences:
                    sequences[seq] = header

            # Replace spaces with underscores, keep all other chars
            header_text = line[1:].strip().replace(" ", "_")

            # Add DATABASE| prefix
            header = f"DATABASE|{header_text}"

            seq_lines = []

        else:
            # Keep all amino acid letters (including U, O, B, Z, X)
            clean_line = line.strip().upper()
            seq_lines.append(clean_line)

    # Save last sequence
    if header and seq_lines:
        seq = "".join(seq_lines).upper()
        if seq not in sequences:
            sequences[seq] = header

# === Write cleaned FASTA ===
with open(output_file, "w") as f:
    for seq, header in sequences.items():
        f.write(f">{header}\n{seq}\n")

print(f"✅ Clean FASTA saved as: {output_file}")
print(f"✅ Total unique sequences: {len(sequences)}")
