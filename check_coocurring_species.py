import sys

# This script identifies and prints co-occurring organisms between two sets of FASTA sequences.
# It takes two FASTA files as input and extracts species names from the sequence IDs.
# The common species names between the two files are then printed.

# Function to iterate over each file, looks for FASTA index '>' and extracts species name.
def extract_organism_names(file_path):
    organism_names = set()
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line[1:].strip()
                organism_name_parts = header.split('[')
                if len(organism_name_parts) >= 2:
                    organism_name = organism_name_parts[1].split(']')[0]
                    organism_names.add(organism_name)
    return organism_names

# Provide file paths as command-line arguments.
file_path_1 = sys.argv[1]
file_path_2 = sys.argv[2]

# Extract organism names from each file.
organism_names_1 = extract_organism_names(file_path_1)
organism_names_2 = extract_organism_names(file_path_2)

# Find common organisms between the two files.
common_organisms = organism_names_1 & organism_names_2

# Print common organisms.
print("Common Organisms:")
for organism in common_organisms:
    print(organism)
