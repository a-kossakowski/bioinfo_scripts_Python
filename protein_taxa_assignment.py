import os
import sys

#assigns order of taxonomy for each protein by searchin for it inside
#proteome file and assigning proteome ID, which is transfered to taxa name

#provide path and name of the file

#proteomes path
proteome_dir = "/home/proteome/"
#input file as a list of IDs separated by '\n'
id_list_file = sys.argv[1]

#search for proteome files in .faa format (.fasta)
proteome_files = [file for file in os.listdir(proteome_dir) if file.endswith(".fasta")]

id_list = []

#read ids from list in file
with open(id_list_file, 'r') as id_list_file:
    id_list = [line.strip() for line in id_list_file]

id_proteome_map = {}

#initialize ID-proteome mapping
for proteome_file in proteome_files:
    #join paths to proteomes and proteome files with axtensions specified in 'proteome_files' variable
    proteome_path = os.path.join(proteome_dir, proteome_file)
    with open(proteome_path, 'r') as proteome:
        for line in proteome:
            if line.startswith('>'):
                id = line[1:].split()[0] 
                if id in id_list:
                    id_proteome_map[id] = proteome_file

#print ID and proteome for each sequence
for id in id_list:
    #if no corresponting proteome, prints 'Not found'
    proteome_file = id_proteome_map.get(id, "Not found")
    print(f"{id}\t{proteome_file}")
