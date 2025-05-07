#!/usr/bin/env python3
"""
SNIPz.py
--------
2025-05-07 13:39:26
Latest update. Now it also handles gzipped files.

This script processes FASTA files (plain or gzipped) to filter sequences based on specific criteria:
1. Removes sequences containing non-IUPAC characters.
2. Removes terminal stop codons (`*`) from sequences.
3. Outputs accepted and removed sequences into separate files.
4. Generates a report summarizing the results.

Features:
- Supports both single and multiple FASTA files.
- Handles gzipped input and output files.
- Generates a detailed TSV report for multiple files.

Usage:
- For a single file:
    python SNIPz.py -s input.fasta
    python SNIPz.py -s input.fasta.gz --gz
- For multiple files:
    python SNIPz.py -m "*.fasta"
    python SNIPz.py -m "*.gz" --gz

Dependencies:
- Biopython
"""

import os
import glob
import re
from Bio import SeqIO
import argparse
import gzip

# Regular expression to extract numbers for sorting
numbers = re.compile(r'(\d+)')

def numericalSort(value):
    """
    Sorts file names numerically to ensure proper order.
    """
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def input_files(ftype):
    """
    Retrieves a sorted list of files matching the given extension.
    """
    list_of_files = sorted((glob.glob(ftype)), key=numericalSort)
    return list_of_files

def seqs_extractor(fasta_file, gzipped=False):
    """
    Extracts sequence names and sequences from a FASTA file.
    Supports both plain and gzipped FASTA files.
    """
    names = []
    seqs = []
    
    # Use gzip.open if the file is gzipped
    if gzipped:
        with gzip.open(fasta_file, "rt") as handle:  # Open in text mode
            for seq_record in SeqIO.parse(handle, "fasta"):
                names.append(seq_record.description)
                seqs.append(seq_record.seq)
    else:
        with open(fasta_file, "r") as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                names.append(seq_record.description)
                seqs.append(seq_record.seq)
    
    return names, seqs

# Command-line argument parser
parser = argparse.ArgumentParser()

# Argument for a single FASTA file
parser.add_argument("-s",
                    help="Specify a single fasta file (target.faa). Accepts unaligned or aligned file.")

# Argument for multiple FASTA files
parser.add_argument("-m", default="",
                    help="Specify the extension for multiple fasta files (e.g. '*.fasta'). Accepts unaligned or aligned files.")

# Flag to indicate gzipped input files
parser.add_argument("--gz", action='store_true',
                    help="Specify if the input files are gzipped. This will be used for both single and multiple file inputs.")

# Parse all arguments
args = parser.parse_args()

# Redefine the parameters
single_file = args.s
ext = args.m
is_gzipped = args.gz

# List of accepted amino acids (IUPAC characters)
accepted_amino = ["A", "a", "c", "C", "d", "D", "e", "E", "f", "F", "g", "G",
                  "h", "H", "i", "I", "k", "K", "l", "L", "m", "M", "n", "N", "p", "P", "q",
                  "Q", "r", "R", "s", "S", "t", "T", "v", "V", "w", "W", "y", "Y", "-",
                  "."]

# Ensure only one of -s or -m is specified
if ext != "" and single_file != None:
    print("You can't specify both a single file and a list of files.")
    print("Please, choose only -s (single file) or -m (multiple files).")
    exit()

# Process a single FASTA file
if ext == "" and single_file != None:
    # Determine the base name of the file
    if is_gzipped:
        real_fasta = ".".join(single_file.split(".")[:-2])
    else:
        real_fasta = ".".join(single_file.split(".")[:-1])
    
    # Extract sequences from the file
    names, seqs = seqs_extractor(single_file, gzipped=is_gzipped)
    
    # Preprocess sequences: remove trailing whitespace and terminal stop codons
    original_seqs = seqs.copy()
    seqs = [seq.rstrip().rstrip('*') for seq in seqs]
    changes_made = any(original_seq != processed_seq for original_seq, processed_seq in zip(original_seqs, seqs))

    # Identify sequences with non-IUPAC characters
    bad_seqs_n = []
    bad_seqs_s = []
    for name, seq in zip(names, seqs):
        for amino in seq:
            if amino not in accepted_amino:
                bad_seqs_n.append(name)
                bad_seqs_s.append(seq)
                print("seq ", name, "has the following non-IUPAC character: ", amino)
                break

    # Separate good and bad sequences
    good_seqs_n = [seq for seq in names if seq not in bad_seqs_n]
    good_seqs_s = [seq for seq in seqs if seq not in bad_seqs_s]
    
    print(f'Results for file {single_file}')
    print("Total seqs: ", len(names))
    print("Accepted seqs: ", len(good_seqs_n))
    print("Removed seqs: ", len(bad_seqs_n))
   
    # Write results to output files
    if len(bad_seqs_n) > 0:
        # File names for removed and accepted sequences
        removed_file = str(real_fasta) + "_removed." + str((single_file.split(".")[-2]))
        accepted_file = str(real_fasta) + "_accepted." + str((single_file.split(".")[-2]))
        if is_gzipped:
            removed_file += ".gz"
            accepted_file += ".gz"            
            with gzip.open(removed_file, "wt") as f1:
                for name_b, seq_b in zip(bad_seqs_n, bad_seqs_s):
                    f1.write(">" + str(name_b) + "\n")
                    f1.write(str(seq_b) + "\n")
            with gzip.open(accepted_file, "wt") as f2:
                for name_g, seq_g in zip(good_seqs_n, good_seqs_s):
                    f2.write(">" + str(name_g) + "\n")
                    f2.write(str(seq_g) + "\n")
        else:
            removed_file = str(real_fasta) + "_removed." + str((single_file.split(".")[-1]))
            accepted_file = str(real_fasta) + "_accepted." + str((single_file.split(".")[-1]))
            with open(removed_file, "w") as f1:
                for name_b, seq_b in zip(bad_seqs_n, bad_seqs_s):
                    f1.write(">" + str(name_b) + "\n")
                    f1.write(str(seq_b) + "\n")
            with open(accepted_file, "w") as f2:
                for name_g, seq_g in zip(good_seqs_n, good_seqs_s):
                    f2.write(">" + str(name_g) + "\n")
                    f2.write(str(seq_g) + "\n")

    elif not changes_made: 
        print("There are no sequences with non-IUPAC characters nor terminal stop codons characters in your fasta file.")
    
    else: 
        print('There are no sequences with non-IUPAC characters in your fasta file, but terminal stop codons characters were found.')
        accep_file_prepr = str(real_fasta) + "_nsc." + str((single_file.split(".")[-2]))
        if is_gzipped:
            accep_file_prepr += ".gz"  
            with gzip.open(accep_file_prepr, "wt") as f3:
                for names_a, seqs_a in zip(names, seqs):
                    f3.write(">" + str(names_a) + "\n")
                    f3.write(str(seqs_a) + "\n")
        else:
            with open(accep_file_prepr, "w") as f3:
                for names_a, seqs_a in zip(names, seqs):
                    f3.write(">" + str(names_a) + "\n")
                    f3.write(str(seqs_a) + "\n")

# Process multiple FASTA files
if ext != "" and single_file == None:
    # Create directories for accepted and removed sequences
    if not os.path.exists('Seqs_Accepted'):
        os.mkdir('Seqs_Accepted') 
    if not os.path.exists('Seqs_Removed'):
        os.mkdir('Seqs_Removed')
    
    # Get the list of FASTA files
    list_of_fastas = input_files(ext)
    if list_of_fastas == []:
        print("No files found with the extension provided.")
        exit()
    
    # Open the report file
    with open("SNIP_report_multiple_files.tsv", "w") as out_f:
        out_f.write("File\tTotal_Seqs\tSeqs_accepted\tSeqs_removed\tterm_*\tOutput_File\n")
    
        for fasta in list_of_fastas:
            # Check if the file is gzipped
            is_gzipped = fasta.endswith(".gz")
            real_fasta = ".".join(fasta.split(".")[:-2]) if is_gzipped else ".".join(fasta.split(".")[:-1])
            
            names, seqs = seqs_extractor(fasta, gzipped=is_gzipped)
            # Preprocess sequences: remove trailing whitespace and asterisks
            original_seqs = seqs.copy()
            seqs = [seq.rstrip().rstrip('*') for seq in seqs]
            changes_made = any(original_seq != processed_seq for original_seq, processed_seq in zip(original_seqs, seqs))

            bad_seqs_n = []
            bad_seqs_s = []
            for name, seq in zip(names, seqs):
                for amino in seq:
                    if amino not in accepted_amino:
                        bad_seqs_n.append(name)
                        bad_seqs_s.append(seq)
                        print("seq ", name, "has the following non-IUPAC character: ", amino)
                        break

            good_seqs_n = [seq for seq in names if seq not in bad_seqs_n]
            good_seqs_s = [seq for seq in seqs if seq not in bad_seqs_s]

            print(f'Results for file {fasta}')
            print("Total seqs: ", len(names))
            print("Accepted seqs: ", len(good_seqs_n))
            print("Removed seqs: ", len(bad_seqs_n))

            if len(bad_seqs_n) > 0:
                if is_gzipped:
                    removed_file = "./Seqs_Removed/" + str(real_fasta) + "_removed." + str((fasta.split(".")[-2]))
                    accepted_file = "./Seqs_Accepted/" + str(real_fasta) + "_accepted." + str((fasta.split(".")[-2]))
                else:
                    removed_file = "./Seqs_Removed/" + str(real_fasta) + "_removed." + str((fasta.split(".")[-1]))
                    accepted_file = "./Seqs_Accepted/" + str(real_fasta) + "_accepted." + str((fasta.split(".")[-1]))
                # Write the bad sequences to a new file
                if is_gzipped:
                    removed_file += ".gz"
                    accepted_file += ".gz"
                    with gzip.open(removed_file, "wt") as f1:
                        for name_b, seq_b in zip(bad_seqs_n, bad_seqs_s):
                            f1.write(">" + str(name_b) + "\n")
                            f1.write(str(seq_b) + "\n")
                    with gzip.open(accepted_file, "wt") as f2:
                        for name_g, seq_g in zip(good_seqs_n, good_seqs_s):
                            f2.write(">" + str(name_g) + "\n")
                            f2.write(str(seq_g) + "\n")
                else:
                    with open(removed_file, "w") as f1:
                        for name_b, seq_b in zip(bad_seqs_n, bad_seqs_s):
                            f1.write(">" + str(name_b) + "\n")
                            f1.write(str(seq_b) + "\n")
                    with open(accepted_file, "w") as f2:
                        for name_g, seq_g in zip(good_seqs_n, good_seqs_s):
                            f2.write(">" + str(name_g) + "\n")
                            f2.write(str(seq_g) + "\n")
            
            elif not changes_made: 
                print("There are no sequences with non-IUPAC characters nor terminal stop codons characters in your fasta file.")
                os.system(f"cp {fasta} Seqs_Accepted/")
            else: 
                print('There are no sequences with non-IUPAC characters in your fasta file, but terminal stop codons characters were found.')
                if is_gzipped:
                    accep_file_prepr = "./Seqs_Accepted/" + str(real_fasta) + "_nsc." + str((fasta.split(".")[-2]))
                else:
                    accep_file_prepr = "./Seqs_Accepted/" + str(real_fasta) + "_nsc." + str((fasta.split(".")[-1]))
                if is_gzipped:
                    accep_file_prepr += ".gz"
                    with gzip.open(accep_file_prepr, "wt") as f3:
                        for names_a, seqs_a in zip(names, seqs):
                            f3.write(">" + str(names_a) + "\n")
                            f3.write(str(seqs_a) + "\n")
                else:
                    with open(accep_file_prepr, "w") as f3:
                        for names_a, seqs_a in zip(names, seqs):
                            f3.write(">" + str(names_a) + "\n")
                            f3.write(str(seqs_a) + "\n")
            
            term_SC = "yes" if changes_made else "no"
            # Determine which output file was actually created
            if len(bad_seqs_n) > 0:
                output_file = accepted_file  # _accepted file is created when sequences are removed
            elif changes_made:
                output_file = accep_file_prepr  # _nsc file is created when stop codons are removed
            else:
                output_file = f"./Seqs_Accepted/{os.path.basename(fasta)}"  # Original file is copied

            # Write to report
            out_f.write(f'{fasta}\t{len(names)}\t{len(good_seqs_n)}\t{len(bad_seqs_n)}\t{term_SC}\t{os.path.basename(output_file)}\n')
            out_f.flush()  # Ensure data is written to the file