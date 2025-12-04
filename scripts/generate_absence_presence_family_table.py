#!/usr/bin/env python3

"""
This script loads a TSV summary file of species, then scans all .txt files
in a given folder. Each TXT file contains a list of species where genus
and species are joined with a hyphen (e.g. "Anolis-carolinensis").
The script adds one new column per TXT file to the TSV, marking 1 if the
species is present in that file and 0 otherwise. The output is a new TSV.
"""

import os
import pandas as pd

# ----------------------------------------------------------
# Variables declared at the top (required by the user)
# ----------------------------------------------------------
input_tsv = "../data/tree_data/AMPlocator_prediction_100_proteomes_analysis_summary.tsv"
species_lists_folder = "../data/novel_families_interpro/selected_families"     # folder containing the .txt files
output_tsv = "../data/tree_data/AMPlocator_prediction_100_proteomes_analysis_summary_with_presence.tsv"
species_column = "Species_Full_Name"
# ----------------------------------------------------------


# ----------------------------------------------------------
# Step 1: Load the main TSV table
# ----------------------------------------------------------
df = pd.read_csv(input_tsv, sep="\t", dtype=str)

# Standardize species names in the main table
# (strip whitespace and ensure consistent capitalization if needed)
df[species_column] = df[species_column].str.strip()

# ----------------------------------------------------------
# Step 2: Iterate over all TXT files in the folder
# ----------------------------------------------------------
for filename in os.listdir(species_lists_folder):

    # Process only .txt files
    if not filename.endswith(".txt"):
        continue

    file_path = os.path.join(species_lists_folder, filename)
    column_name = os.path.splitext(filename)[0]   # filename without extension

    # Read all species from the file
    with open(file_path, "r") as f:
        species_raw = [line.strip() for line in f if line.strip()]

    # Normalize species: replace hyphen by space
    # Example: "Anolis-carolinensis" → "Anolis carolinensis"
    species_normalized = [s.replace("-", " ") for s in species_raw]

    # Convert to a set for faster membership checks
    species_set = set(species_normalized)

    # Create the new column marking presence/absence
    df[column_name] = df[species_column].apply(lambda x: "1" if x in species_set else "0")


# ----------------------------------------------------------
# Step 3: Save the resulting table
# ----------------------------------------------------------
df.to_csv(output_tsv, sep="\t", index=False)

print(f"Done. Output written to: {output_tsv}")
