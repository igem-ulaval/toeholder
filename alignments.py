#!/usr/bin/env python3
# coding: utf-8

# # Specificity check
# 
# This script will read a list of toehold candidates and align them to a reference genome. This will allow us to see if there are any toeholds that could have promiscuous binding.


# Load libraries
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
import os
import csv
from collections import OrderedDict
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shutil import copyfile
import glob
import argparse
from toeholder_helper_functions import *

parser = argparse.ArgumentParser(description = 'This script will read a list of toehold candidates and align them to a reference genome. This will allow us to see if there are any toeholds that could have promiscuous binding.')

parser.add_argument('-w', type = str, help = 'The working directory', dest = 'working_directory')
parser.add_argument('-r', type = str, help = 'The file containing the list of reference genomes and tags for each', dest = 'reference_list')
parser.add_argument('-p', type = str, help = 'Minimum sequence identity percentage between the switch recognition sequence and any matches', dest = 'pct_ident', default = 0.95)
parser.add_argument('-e', type = str, help = 'E-value threshold for accepted matches', dest = 'evalue', default = 1e-06)

args = parser.parse_args()

# Load input variables 
working_directory = args.working_directory
reference_list = args.reference_list
pct_ident = args.pct_ident
evalue = args.evalue

# Load the toehold sequences 
handle = open(os.path.join(working_directory, 'all_toeholds_results.txt'), 'r')
reader = csv.reader(handle, delimiter = '\t')

toehold_list = []

for line in reader:
    # Skip header
    if line[1] == 'Structure':
        continue
        
    # Save toehold indices and sequences
    toehold_list.append([line[8], line[2].replace('U', 'T')])
    
df = pd.DataFrame(np.array(toehold_list), columns=['Index', 'Sequence'])
handle.close()

# Write a fasta file with all those sequences
records = []

for df_index, row in df.iterrows():
    index = row[0]
    sequence = row[1]
    
    record = SeqRecord(Seq(sequence), id = 'toehold_' + index, description = '')
    records.append(record)
    
toeholds_fasta = os.path.join(working_directory, 'toehold_seqs.fasta')
SeqIO.write(records, toeholds_fasta, 'fasta')

# Load the toehold sequences 
handle = open(os.path.join(working_directory, 'all_toeholds_results.txt'), 'r')
reader = csv.reader(handle, delimiter = '\t')

toehold_list = []

for line in reader:
    # Save header
    if line[1] == 'Structure':
        header = line
        continue
        
    
    # Save toehold indices and sequences
    toehold_list.append(line)
    
handle.close()
new_df = pd.DataFrame(np.array(toehold_list), columns=header)

#### Need to loop through the pairs of reference genomes and tags ####
handle = open(reference_list, 'r')
reader = csv.reader(handle, delimiter = '\t')

for entry in reader:
    reference_genome_file = entry[0]
    tag = entry[1]
    
    print('Working with the genome for  ' + tag)
    
    # Run the function to align and add the columns
    new_df = count_genome_matches(toeholds_fasta, working_directory, toehold_list, reference_genome_file, tag, new_df, pct_ident, evalue)
    
    print('--------')

# Save the data frame
new_df.to_csv(path_or_buf=os.path.join(working_directory, 'all_toeholds_results_genome_matches.txt'), sep = '\t', index = False)

