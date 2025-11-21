
#!/usr/bin/env python3
import os
from collections import defaultdict

# ==================================================================
# VARS
# ==================================================================
ID_FILE = "../data/shared_novel_seqs_ids_clean_no_header.txt"
PREDICTIONS_FOLDER = "/home/mario/projects/100_proteomes_2/1_amp_locator_prediction/predictions"
SUMMARY_TABLE = "../data/tree_data/AMPlocator_prediction_100_proteomes_analysis_summary.tsv"
NEW_COLUMN_NAME = "Novel_AMPs"   # Name of new column to add
OUTPUT_TABLE = "../data/tree_data/AMPlocator_prediction_100_proteomes_analysis_summary_updated.tsv"
# ==================================================================


def load_ids(id_file):
    """Load sequence IDs from file, one ID per line."""
    with open(id_file, "r") as f:
        return [line.strip() for line in f if line.strip()]


def parse_fasta_for_ids(fasta_file):
    """Parse a FASTA file and return a set of IDs found in header lines."""
    found = set()
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                found.add(seq_id)
    return found


def count_ids_per_species(ids, predictions_folder):
    """Search IDs in all FASTA files and count how many per species."""
    ids_to_species = defaultdict(list)

    fasta_files = [
        os.path.join(predictions_folder, f)
        for f in os.listdir(predictions_folder)
        if f.lower().endswith(".fasta") or f.lower().endswith(".fa")
    ]

    if not fasta_files:
        print("ERROR: No FASTA files found in predictions folder.")
        return {}

    for fasta in fasta_files:
        basename = os.path.basename(fasta)
        species = basename.split("_")[0]  # e.g. Acinetobacter-baumannii

        fasta_ids = parse_fasta_for_ids(fasta)

        for seq_id in ids:
            if seq_id in fasta_ids:
                ids_to_species[seq_id].append(species)

    # Aggregate counts
    species_counts = defaultdict(int)
    missing_ids = []
    duplicated_ids = []

    for seq_id in ids:
        if seq_id not in ids_to_species:
            missing_ids.append(seq_id)
        elif len(ids_to_species[seq_id]) > 1:
            duplicated_ids.append(seq_id)
        else:
            species_counts[ids_to_species[seq_id][0]] += 1

    # Warnings
    if missing_ids:
        print("\nWARNING – IDs NOT FOUND IN ANY FASTA:")
        for x in missing_ids:
            print(f"  {x}")

    if duplicated_ids:
        print("\nWARNING – IDs FOUND IN MULTIPLE SPECIES:")
        for x in duplicated_ids:
            print(f"  {x}")

    return species_counts


def add_column_to_summary(summary_file, species_counts, new_column_name, output_file):
    """
    Add a new column to the summary file based on counts.
    """
    output_lines = []

    with open(summary_file, "r") as f:
        header = f.readline().rstrip("\n")
        columns = header.split(",")

        # Add the new column
        columns.append(new_column_name)
        output_lines.append(",".join(columns))

        # Keep track of species found in summary
        species_seen = set()

        for line in f:
            line = line.rstrip("\n")
            fields = line.split(",")
            species_name = fields[0]   # e.g. "Acinetobacter baumannii"
            species_seen.add(species_name)

            # Convert species name from summary to FASTA-style for lookup
            fasta_style_name = species_name.replace(" ", "-")

            value = species_counts.get(fasta_style_name, 0)
            fields.append(str(value))
            output_lines.append(",".join(fields))

    # Check for species in counts that were not in the summary
    extra_species = set(species_counts.keys()) - {s.replace(" ", "-") for s in species_seen}
    if extra_species:
        print("\nWARNING – some species found in FASTA not present in summary table:")
        for sp in extra_species:
            print(f"  {sp}")

    # Write updated file
    with open(output_file, "w") as f:
        for line in output_lines:
            f.write(line + "\n")

    print(f"\nUpdated table written to: {output_file}")


if __name__ == "__main__":
    
    ids = load_ids(ID_FILE)
    species_counts = count_ids_per_species(ids, PREDICTIONS_FOLDER)
    add_column_to_summary(SUMMARY_TABLE, species_counts, NEW_COLUMN_NAME, OUTPUT_TABLE)
