import sys

#based on FASTA sequences IDs identifies species, if presents, and prints coocuring
#organisms between two sets of FASTA sequences


#function to iterate over each file, looks for FASTA incice '>' and extracts species name
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

#provide file paths
file_path_1 = sys.argv[1]
file_path_2 = sys.argv[2]

#extract organism names from each file
organism_names_1 = extract_organism_names(file_path_1)
organism_names_2 = extract_organism_names(file_path_2)

#find common organisms
common_organisms = organism_names_1.intersection(organism_names_2)

#print common organisms
print("Common Organisms:")
for organism in common_organisms:
    print(organism)
