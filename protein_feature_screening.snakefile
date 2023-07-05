import os

#this workflow allows for prediction of fundamental protein features

# PfamScan: Program for protein domain analysis
rule pfamscan:
    input:
        input_file="{file}.fasta"
    output:
        output_file=temp("{file}/{file}.pfamscan")
    shell:
        "/home/tools/PfamScan/pfam_scan.pl -fasta {input.input_file} -outfile {output.output_file} -dir /home/db/PFAM/"
    run:
        touch(output.output_file)

# SEG: Program for detecting low complexity regions in proteins
rule seg:
    input:
        input_file="{file}.fasta"
    output:
        output_file=temp("{file}/{file}.lcr")
    shell:
        "/home/tools/seg/seg -x {input.input_file} >> {output.output_file}"
    run:
        touch(output.output_file)

# pepcoil: Program for predicting coiled coil regions in proteins
rule pepcoil:
    input:
        input_file="{file}.fasta"
    output:
        output_file=temp("{file}/{file}.pepcoil")
    shell:
        "pepcoil {input.input_file} -window 28 -outfile {output.output_file}"
    run:
        touch(output.output_file)

# iupred: Program for predicting intrinsically disordered regions in proteins
rule iupred:
    input:
        input_file="{file}.fasta"
    output:
        output_file=temp("{file}/{file}.diso")
    shell:
        "python3 /home/tools/iupred3/iupred3.py {input.input_file} long >> {output.output_file}"
    run:
        touch(output.output_file)

# WoLFPSort: Program for subcellular localization prediction of proteins in fungi
rule wolfpsort:
    input:
        input_file="{file}.fasta"
    output:
        output_file=temp("{file}/{file}.loc")
    shell:
        "/home/tools/WoLFPSort/bin/runWolfPsortSummary fungi < {input.input_file} >> {output.output_file}"
    run:
        touch(output.output_file)

# muscle: Program for multiple sequence alignment
rule muscle:
    input:
        input_file="{file}.fasta"
    output:
        output_file=temp("{file}/{file}.msa")
    shell:
        "muscle -in {input.input_file} -out {output.output_file}"
    run:
        touch(output.output_file)

# TMHMM: Program for transmembrane helix prediction in proteins
rule tmhmm:
    input:
        input_file="{file}.fasta"
    output:
        output_file=temp("{file}/{file}.tmhmm")
    shell:
        "/home/tools/TMHMM2/tmhmm-2.0c/bin/tmhmm < {input.input_file} >> {output.output_file}"
    run:
        touch(output.output_file)

# HHblits: Program for remote homology detection using profile-profile alignment
rule hhblits:
    input:
        input_file="{file}.msa"
    output:
        output_file=temp("{file}/{file}.hhr")
    shell:
        "/opt/hhsuite/build/bin/hhblits -i {input.input_file} -o {output.output_file} -d /home/db/hhsuitedb/pfam -M 50"
    run:
        touch(output.output_file)

# Rule to create separate folders and run all programs
rule process_files:
    input:
        "input_files/{file}.fasta"
    output:
        directory(temp("{file}"))
    shell:
        """
        mkdir -p {output.directory}
        """
    run:
        for file in os.listdir("input_files"):
            if file.endswith(".fasta"):
                file_prefix = os.path.splitext(file)[0]
                run_pfamscan(file_prefix)
                run_seg(file_prefix)
                run_pepcoil(file_prefix)
                run_iupred(file_prefix)
                run_wolfpsort(file_prefix)
                run_muscle(file_prefix)
                run_tmhmm(file_prefix)
                run_hhblits(file_prefix)
