import sys

#retrives FASTA sequences from given directory {file_path} based on given {pattern}

#this function will search a .fasta file for sequences that contain a given pattern
def retrieve_sequences_with_pattern(file_path, pattern):
    sequences = []
    with open(file_path, 'r') as file:
        sequence = ''
        header = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence and header and pattern in header:
                    sequences.append((header, sequence))
                header = line
                sequence = ''
            else:
                sequence += line
        #check last sequence after the loop ends
        if sequence and header and pattern in header:
            sequences.append((header, sequence))

    return sequences

#file path
file_path = sys.argv[1]

#user input pattern
pattern = sys.argv[2]

#retrieve sequences with the pattern
sequences = retrieve_sequences_with_pattern(file_path, pattern)

#print sequences in FASTA format
for header, sequence in sequences:
    print(header)
    print(sequence)
