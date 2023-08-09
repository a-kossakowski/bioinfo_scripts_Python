import sys

# This script retrieves FASTA sequences from a given directory based on a given pattern.
# It searches a .fasta file for sequences that contain a specified pattern in the header
# and prints the matching sequences in FASTA format.

# Function to search a .fasta file for sequences that contain a given pattern.
def retrieve_sequences_with_pattern(file_path, pattern):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        header = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # If a sequence and header are found and the pattern is in the header, add to sequences.
                if sequence and header and pattern in header:
                    sequences.append((header, sequence))
                header = line
                sequence = ''
            else:
                # Concatenate lines to form the sequence.
                sequence += line
        # Check the last sequence after the loop ends.
        if sequence and header and pattern in header:
            sequences.append((header, sequence))

    return sequences

# File path and user input pattern are provided as command-line arguments.
file_path = sys.argv[1]
pattern = sys.argv[2]

# Retrieve sequences with the specified pattern.
sequences = retrieve_sequences_with_pattern(file_path, pattern)

# Print sequences in FASTA format.
for header, sequence in sequences:
    print(header)
    print(sequence)
