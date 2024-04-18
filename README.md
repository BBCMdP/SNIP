- [SNIP](#snip)
      - [Accesory script to remove Sequences with Non-IUPAC characters in multifasta file(s)](#accesory-script-to-remove-sequences-with-non-iupac-characters-in-multifasta-files)
    - [Requirements](#requirements)
    - [Single fasta file](#single-fasta-file)
    - [Multiple fasta files](#multiple-fasta-files)


#### Auxiliary script to remove Sequences with Non-IUPAC characters in multifasta file(s) 

A common issue derived from poorly annotated genomes is that resulting proteomes may come with characters that do not correspond to actual amino acids commonly found in natural proteins (what we define here as non-IUPAC characters). 
The most common scenario is the low definition of bases in the genome, typically annotated as nucleotide "N", for which translation of the predicted coding sequence may lead to an undefined residue, often times represented with symbol "X" in the protein sequence. 

Some downstream applications (even MSA generation, profile to sequence comparisons, and even phylogeny reconstruction) can't deal with for non-IUPAC characters, and results are compromised. Moreover, some proteomes show a high "contamination" with non-IUPAC annotations, making this a problem that can easily become a nuissance. In addition, stop codon character * is typically found in complete proteome annotations, and can also interfere with several downstream applications.  

In order to purge the protein datasets from sequences with non-IUPAC annotations, we have developed the Sequence Non-IUPAC Purge script, SNIP. SNIP is a python script that can deal with single or multiple fasta files, to automatically detect and remove sequences with non-IUPAC characters. SNIP do accept gap characters (such as "-" or "."), so either aligned or unaligned fasta files are accepted. The script will also process terminal stop codon character "*", and specifically remove it from the sequence. Note, however, that if a non terminal stop codon character is found in a sequence, it will be considered as non-IUPAC character, resulting in the removal of the sequence.

>**We recommend applying this script before to run MuFasA and Seqrutinator.** 

### Requirements

SNIP requires Biopython, and accepts by argument either single (`-s`) or multiple (`-m`) fasta files. For multiple file, it is required that all of them have the same extension. 

### Single fasta file
The output will depend on the results. If running a single fasta file:

`python3 SNIP.py -s input.fasta` 

and sequences with non-IUPAC characters are found, two output fasta files will be generated: `input_accepted.fasta` and `input_removed.fasta`. Of course, each file includes the sequences without and with non-IUPAC characters, respectively. The terminal will print out the amount of sequences (total, accepted and removed). It will also show which is the non-IUPAC character found for each removed sequence.

```
Results for file input.fasta
Total seqs:  65809
Accepted seqs:  65368
Removed seqs:  441
```
Note, in addition, that sequences in both files will not have terminal stop codon character "*".

Now, if no sequences with non-IUPAC characters are found, but terminal stop codons are detected, the resulting sequences, without these characters, will be written in a new file with extension `_nsc.fasta` (from no stop codons).

Finally, if neither non-IUPAC or terminal stop codon characters are found, only a printout in terminal is shown expliciting it (no further files are written). 

### Multiple fasta files
SNIP can be applied to multiple fasta files. The only requirements are that all files should (i) be in the same folder, and (ii) have the same extension. For example, can be run like this: 

`python3 SNIP.py -m '*.fa`  

in a folder with five fasta files with the provided extension:

```
SNIP.py
sfa.fa
crh.fa
smu.fa
tpl.fa
cri.fa
```

Will produce two folders, `/Seqs_Accepted` and `/Seqs_Removed`, plus a summary file `Report_multiple_files.tsv`. The files generated are the same and created conditionally as described before. The difference is that all files `*_removed.fa` will be moved to the `/Seqs_Removed` folder (if any). The `*_accepted.fa` and `*_ncs.fa` (if any) are moved to `/Seqs_Accepted`, as well as all files that were not modified by the script. This is to simplify the output recovery for the user. 

The `Report_multiple_files.tsv` summarizes the results, for example:

| File           | Total_Seqs | Seqs_accepted | Seqs_removed | term_* |
|----------------|------------|---------------|--------------|--------|
| crh.fa         | 19527      | 19514         | 13           | yes    |
| cri.fa | 75253      | 75253         | 0            | yes    |
| sfa.fa | 45611      | 45611         | 0            | no     |
| smu.fa | 27137      | 24567         | 2570         | no     |
| tpl.fa | 65809      | 65368         | 441          | yes    |
> The column term_* indicates if terminal stop codon characters were found

The user can identify here that for input file `crh.fa` has 19527 seqs, of which 13 have non-IUPAC characters (and are removed), and terminal stop codon characters were found. Thus, `crh_accepted.fa` and `crh_removed.fa` files are created and saved in the proper folder. The same goes for `smu.fa` (note that around 10% of the sequences have non-IUPAC characters) and `tpl.fa`. 
The case of `cri.fa` will result in a file `cri_nsc.fa`, since there were no sequences with non-IUPAC characters, but terminal stop codons were found. 
Finally, a copy of `sfa.fa` will be found in the `/Seqs_Accepted` folder, since neither non-IUPAC or terminal stop codons were found. 
