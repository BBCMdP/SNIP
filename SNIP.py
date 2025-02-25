# -*- coding: utf-8 -*-

###############################################################################
#                         Script Description                                  #
# name:                                                                       #
#     SNIP.py                                                                 #
#                                                                             #
# description:                                                                #
#     detects and removes sequences with non-IUPAC characters from multifasta #
#     file(s)                                                                 #
#                                                                             #
# Authors:                                                                    #
#     Agustin Amalfitano, Nico Stocchi, Fernando Villarreal                   #
#                                                                             #
# Date:                                                                       #
#     2024-04-18                                                              #
###############################################################################

from Bio import SeqIO
import argparse
import re
import glob
import os

numbers = re.compile(r'(\d+)')

def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def input_files(ftype):
    list_of_files = sorted((glob.glob(ftype)), key=numericalSort)
    return list_of_files

def seqs_extractor(fasta_file):
    names = []
    seqs = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        names.append(seq_record.id)
        seqs.append(seq_record.seq)
    return names, seqs

# command line parameters
parser = argparse.ArgumentParser()

# input sinlge fasta file
parser.add_argument("-s",
                    help="Specify a single fasta file (target.faa). Accepts unaligned or aligned file.")

parser.add_argument("-m", default="",
                    help="Specify the extension for multiple fasta files (e.g. '*.fasta'). Accepts unaligned or aligned files.")

#all parameters in args
args = parser.parse_args()

# redefine the parameters names
single_file = args.s
ext = args.m

accepted_amino = ["A", "a", "c", "C", "d", "D", "e", "E", "f", "F", "g", "G",
                  "h", "H", "i", "I", "k", "K", "l", "L", "m", "M", "n", "N", "p", "P", "q",
                  "Q", "r", "R", "s", "S", "t", "T", "v", "V", "w", "W", "y", "Y", "-",
                  "."]

if ext != "" and single_file != None:
    print("You can't specify both a single file and a list of files.")
    print("Please, choose only -i (single file) or -m (multiple files).")
    exit()

if ext == "" and single_file != None:
    real_fasta = str((single_file.split(".")[:-1])[0])
    names, seqs = seqs_extractor(single_file)
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
    
    print(f'Results for file {single_file}')
    print("Total seqs: ", len(names))
    print("Accepted seqs: ", len(good_seqs_n))
    print("Removed seqs: ", len(bad_seqs_n))
   
    if len(bad_seqs_n) > 0:
        removed_file = str(real_fasta) + "_removed." + str((single_file.split(".")[-1]))
        accepted_file = str(real_fasta) + "_accepted." + str((single_file.split(".")[-1]))
       
        f1 = open(removed_file, "w")
        for name_b, seq_b in zip(bad_seqs_n, bad_seqs_s):
            f1.write(">" + str(name_b) + "\n")
            f1.write(str(seq_b) + "\n")
        f1.close()

        f2 = open(accepted_file, "w")
        for name_g, seq_g in zip(good_seqs_n, good_seqs_s):
            f2.write(">" + str(name_g) + "\n")
            f2.write(str(seq_g) + "\n")
        f2.close()

    elif not changes_made: 
        print("There are no sequences with non-IUPAC characters nor terminal stop codons characters in your fasta file.")
    
    else: 
        print('There are no sequences with non-IUPAC characters in your fasta file, but terminal stop codons characters were found.')
        accep_file_prepr = str(real_fasta) + "_nsc." + str((single_file.split(".")[-1]))
        f3 = open(accep_file_prepr, "w")
        for names_a, seqs_a in zip(names, seqs):
            f3.write(">" + str(names_a) + "\n")
            f3.write(str(seqs_a) + "\n")
        f3.close()

if ext != "" and single_file == None:
    if not os.path.exists('Seqs_Accepted'):
        os.mkdir('Seqs_Accepted') 
    if not os.path.exists('Seqs_Removed'):
        os.mkdir('Seqs_Removed')
    list_of_fastas = input_files(ext)
    if list_of_fastas == []:
        print("No files found with the extension provided.")
        exit()
    out_f = open("Report_multiple_files.tsv", "w")
    out_f.write("File\tTotal_Seqs\tSeqs_accepted\tSeqs_removed\tterm_*\tOutput_File\n")
    
    for fasta in list_of_fastas:
        real_fasta = str((fasta.split(".")[:-1])[0])
        names, seqs = seqs_extractor(fasta)
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
            removed_file = "./Seqs_Removed/" + str(real_fasta) + "_removed." + str((fasta.split(".")[-1]))
            accepted_file = "./Seqs_Accepted/" + str(real_fasta) + "_accepted." + str((fasta.split(".")[-1]))
        
            f1 = open(removed_file, "w")
            for name_b, seq_b in zip(bad_seqs_n, bad_seqs_s):
                f1.write(">" + str(name_b) + "\n")
                f1.write(str(seq_b) + "\n")
            f1.close()

            f2 = open(accepted_file, "w")
            for name_g, seq_g in zip(good_seqs_n, good_seqs_s):
                f2.write(">" + str(name_g) + "\n")
                f2.write(str(seq_g) + "\n")
            f2.close()

        elif not changes_made: 
            print("There are no sequences with non-IUPAC characters nor terminal stop codons characters in your fasta file.")
            os.system("cp " + fasta + " Seqs_Accepted/")
        else: 
            print('There are no sequences with non-IUPAC characters in your fasta file, but terminal stop codons characters were found.')
            accep_file_prepr = "./Seqs_Accepted/" + str(real_fasta) + "_nsc." + str((fasta.split(".")[-1]))
            f3 = open(accep_file_prepr, "w")
            for names_a, seqs_a in zip(names, seqs):
                f3.write(">" + str(names_a) + "\n")
                f3.write(str(seqs_a) + "\n")
            f3.close()
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
    
    out_f.close()
