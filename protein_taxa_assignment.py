import os
import sys

# This script assigns the order of taxonomy for each protein by searching for it inside
# a proteome file and assigning a proteome ID, which is transferred to a taxa name.
# It takes a list of protein IDs and searches for corresponding proteome files in a specified directory.

# Path to the directory containing proteome files.
proteome_dir = "/home/proteome/"

# Input file containing a list of protein IDs separated by '\n'.
id_list_file = sys.argv[1]

# Search for proteome files in .faa format (.fasta) within the specified directory.
proteome_files = [file for file in os.listdir(proteome_dir) if file.endswith(".fasta")]

id_list = []

# Read IDs from the list in the file.
with open(id_list_file, 'r') as id_list_file:
    id_list = [line.strip() for line in id_list_file]

id_proteome_map = {}

# Initialize ID-proteome mapping.
for proteome_file in proteome_files:
    # Join paths to proteomes and proteome files with extensions specified in 'proteome_files' variable.
    proteome_path = os.path.join(proteome_dir, proteome_file)
    with open(proteome_path, 'r') as proteome:
        for line in proteome:
            if line.startswith('>'):
                id = line[1:].split()[0]
                if id in id_list:
                    id_proteome_map[id] = proteome_file

# Print ID and proteome for each sequence.
for id in id_list:
    # If no corresponding proteome is found, prints 'Not found'.
    proteome_file = id_proteome_map.get(id, "Not found")
    print(f"{id}\t{proteome_file}")
