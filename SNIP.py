# -*- coding: utf-8 -*-

###############################################################################
#                         Script Description                                  #
# name:                                                                       #
#     bad_seqs_v01.py                                                         #
#                                                                             #
# description:                                                                #
#     bad aminoacid detector , remove sequences with bad aminoacids           #
#                                                                             #
#                                                                             #
# Authors:                                                                    #
#     Agustin Amalfitano, Nico Stocchi                                        #
#                                                                             #
# Date:                                                                       #
#     2019-10-22                                                              #
# Modified for python 3 and output rename by Villarreal Fernando (04/11/2022) #
###############################################################################

from Bio import SeqIO
import argparse


# command line parameters
parser = argparse.ArgumentParser()

# input fastafile
parser.add_argument("-i", default="target.faa",
                    help="Specify the fasta file (-i target.fsa) by def")
#all parameters in args
args = parser.parse_args()

# redefine the parameters names
fasta_file = args.i
real_fasta = str((fasta_file.split(".")[:-1])[0])

accepted_amino = ["A", "a", "c", "C", "d", "D", "e", "E", "f", "F", "g", "G",
    "h", "H", "i", "I", "k", "K", "l", "L", "m", "M", "n", "N", "p", "P", "q",
    "Q", "r", "R", "s", "S", "t", "T", "v", "V", "w", "W", "y", "Y", "*", "-",
    "."]


def seqs_extractor(fasta_file):
    names = []
    seqs = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        names.append(seq_record.id)
        seqs.append(seq_record.seq)
    return names, seqs


def fetching(hits, target, output):

    with open(output, "w") as f:
        for seq_record in SeqIO.parse(target, "fasta"):
            if seq_record.id in hits:
                f.write(">" + str(seq_record.id) + "\n")
                f.write(str(seq_record.seq) + "\n")
    f.close()
    return()


names, seqs = seqs_extractor(fasta_file)
bad_seqs = list()
for name, seq in zip(names, seqs):
    for amino in seq:
        if amino not in accepted_amino:
            bad_seqs.append(name)
            print ("seq ", name, "has the following bad character: ", amino)
            break

print (bad_seqs)
good_seqs = [seq for seq in names if seq not in bad_seqs]


print ("total seqs: ", len(names))
print ("accepted seqs: ", len(good_seqs))
print ("bad seqs: ", len(bad_seqs))
fetching(bad_seqs, fasta_file, str(real_fasta) + "_removed.fsa")
fetching(good_seqs, fasta_file, str(real_fasta) + "_accepted.fsa")
