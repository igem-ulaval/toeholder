#!/usr/bin/env python
# coding: utf-8


#######################################################################################
####            Designing toeholds from a target sequence or genome                ####
####                                                                               #### 
#### This script will receive a target sequence or genome and look for candidates  ####
#### of sites for which we could design toeholds. Candidates will be selected      ####
#### as follows:                                                                   ####
#### - Regions having 2 weak pairs and 1 strong at the base of the candidate       ####
#### - User defines lengths a (unpaired part of the recognition sequence) and      ####
#### b (paired part)                                                               ####
#######################################################################################

# Load libraries
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
from Bio import SeqIO
import os
import csv
from collections import OrderedDict
import numpy as np
import pandas as pd
from shutil import copyfile
import glob
from toeholder_helper_functions import *
from input_variables import *

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Make a copy of the input parameter file in the new directory
copyfile('input_variables.py', os.path.join(output_folder, 'input_variables.py'))

selected = []

# ## 1.- Use a loop to find candidates
# Parse the input sequence
input_seq_record = next(SeqIO.parse(input_seq, 'fasta'))

# If input sequence is DNA, parse and transcribe to look at secondary structure of mRNA
if mol_type == 'DNA':   
    full_input_seq = str(input_seq_record.seq).upper().replace('T', 'U')   
# If it is RNA, prepare a NUPACK run to look at its secondary structure
elif mol_type == 'RNA':
    full_input_seq = str(input_seq_record.seq).upper()
    
# Calculate the secondary structure with NUPACK
output_file = os.path.join(output_folder, 'mRNA.in')
prepare_nupack_input(full_input_seq, output_file)
nupack_mfe(output_file[0:-3])

# Parse with NUPACK
mRNA_dict = parse_nupack(os.path.join(output_folder, 'mRNA.mfe'))
mRNA_structure = mRNA_dict['structure']
mRNA_sequence = mRNA_dict['sequence']

for hairpin_start_pos in range(length_unpaired, len(mRNA_structure) - length_paired):
    # Check the base of the hairpin
    hairpin_base = mRNA_sequence[hairpin_start_pos:hairpin_start_pos + 3]
    
    # Discard cases for which the base of the hairpin does not have two weak pairs and a strong one
    if hairpin_base.count('A') + hairpin_base.count('U') != 2:
        continue
    
    end_pos = hairpin_start_pos + length_paired
    start_pos = hairpin_start_pos - length_unpaired
    sub_seq = mRNA_sequence[start_pos:end_pos]
    sub_struc = mRNA_structure[start_pos:end_pos]
    
    ## Look at GC content at the least stable positions from MD ##

    # After the end of the trigger sequence, we always have AUA and the paired side of the hairpin
    # Disregard these parts of the sequence since the AUA does not change, and the rest is just complementary to the
    # bases to check
    for i in GC_weak:
        if i - 4 > len(sub_seq):
            GC_weak.remove(i)

    # Substract 4 from the positions to account for the GGG and the change from 1-based to 0-based
    nucleotides_GC_weak = [sub_seq[i - 4] for i in GC_weak]

    # Count the number of G's and C's in those nucleotides to check
    count_check_GC_weak = nucleotides_GC_weak.count('G') + nucleotides_GC_weak.count('C')

    ## Look at GC content at the most stable positions from MD ##

    # After the end of the trigger sequence, we always have AUA and the paired side of the hairpin
    # Disregard these parts of the sequence since the AUA does not change, and the rest is just complementary to the
    # bases to check
    for i in GC_strong:
        if i - 4 > len(sub_seq):
            GC_strong.remove(i)

    # Substract 4 from the positions to account for the GGG and the change from 1-based to 0-based
    nucleotides_GC_strong = [sub_seq[i - 4] for i in GC_strong]

    # Count the number of G's and C's in those nucleotides to check
    count_check_GC_strong = nucleotides_GC_strong.count('G') + nucleotides_GC_strong.count('C')

    # Count the number of unpaired bases
    unpaired = sub_struc.count('.')
    if unpaired >= min_unpaired:
        selected.append([unpaired, sub_struc, sub_seq, start_pos + 1, end_pos,
                         length_unpaired, length_paired, 
                         count_check_GC_weak, count_check_GC_strong])
    
tmp = pd.DataFrame(np.array(selected), columns=['Non paired count', 'Structure', 'Sequence', 'Start', 'End', 
                   'Length unpaired trigger', 'Length paired trigger', 'Count GC weak', 'Count GC strong'])    
df = tmp.sort_values(by = 'Non paired count', ascending = False)

df.to_csv(path_or_buf=os.path.join(output_folder, 'toehold_candidates.txt'), sep = '\t', index = False)
header_df = list(df.columns)

# Loop through the folders and gather the information in a table
results = []

# Retrieve the binding energy of the mRNA on its own
mRNA_binding_energy = mRNA_dict['energy']

# Prepare toeholds for each of the sequences
for index, row in df.iterrows():
    sequence = row[2]

    # ## 2.- Generate toeholds for each of the candidates
    toehold_folder = os.path.join(output_folder, str(index))
    generate_toehold(sequence, mol_type, reporter, toehold_folder, length_unpaired, length_paired)

    # ## 3.- Check how well the toeholds bind to the target in the mRNA
    if not os.path.exists(os.path.join(toehold_folder, 'switch1_python.mfe')):
        # This is the case in which the toehold generator function found a stop codon and stopped
        continue
    
    # For all other cases, write down the file that will work as the input for NUPACK
    # Read the file to extract the toehold's sequence
    toehold_dict = parse_nupack(os.path.join(toehold_folder, 'switch1_python.mfe'))
    toehold = toehold_dict['sequence']
    output_file = os.path.join(toehold_folder, 'toehold_mRNA.in')
    
    prepare_nupack_input_two(toehold, mRNA_sequence, output_file)
    
    # Run NUPACK
    nupack_mfe(output_file[0:-3])

    # Use the last toehold to save the length of generated switches
    len_toehold = len(toehold)

    # Retrieve the information from the table about this toehold
    toehold_num = int(index)
    toehold_data = [entry for entry in row]

    # Check if we have a result
    result_path = os.path.join(toehold_folder, 'toehold_mRNA.mfe')

    start_pos = int(toehold_data[3]) + len_toehold
    end_pos = int(toehold_data[4]) + len_toehold

    toehold_data.append(toehold_num)
    
    if os.path.exists(result_path):
        # Parse the output file
        toehold_mRNA_dict = parse_nupack(result_path)
        
        # Add the energy for the toehold binding to the mRNA
        binding_energy = toehold_mRNA_dict['energy']
        toehold_data.append(binding_energy)
        
        # Add to the current row the percentage of residues that are successfully bound
        pair_list = toehold_mRNA_dict['paired_bases']
        percentage_correct = check_matches(start_pos, end_pos, pair_list)
        toehold_data.append(percentage_correct)
        
        # Get the binding energy for the toehold on its own
        toehold_alone_path = os.path.join(toehold_folder, 'switch1_python.mfe')
        toehold_alone_dict = parse_nupack(toehold_alone_path)
        toehold_alone_binding_energy = toehold_alone_dict['energy']
        toehold_data.append(toehold_alone_binding_energy)
        
    else:
        # There are no values for binding energy or percentage of correct bases
        toehold_data.append('NA')
        toehold_data.append('NA')
        toehold_data.append('NA')
    
    # Add the binding energy of the mRNA on its own
    toehold_data.append(mRNA_binding_energy)
    
    # Look at the GC content
    sequence = toehold_data[2]
    gc_content = round(float(sequence.count('G') + sequence.count('C'))*100/len(sequence), 2)
    toehold_data.append(gc_content)
    
    ## Calculate the energy for the RBS linker ##
    
    toehold_mfe = toehold_dict['energy']
    toehold_sequence = toehold_dict['sequence']
    toehold_structure = toehold_dict['structure']
    pair_list = toehold_dict['paired_bases']
    
    ## Find the RBS inside the sequence
    rbs_start = str.find(toehold_sequence, 'AGGAGA')
    
    ## Look for the position that is paired to the position following the RBS
    ## Add 7 because of a shift of 6 positions and a change from 0-based to 1-based
    position_after_rbs = rbs_start + 7
 
    for bound_pair in pair_list:
        if str(position_after_rbs) in bound_pair:
            # Take the first position of the bound pair
            # Subtract 1 to go from 1-based to 0-based but add 1 to remove the final paired
            # base from the left side of the hairpin
            rbs_linker_start = int(bound_pair[0])
            break
    
    ## Use the RBS position to extract the RBS linker
    ## RBS linker goes from the last three paired bases at the top of the hairpin to
    ## the end of the last small hairpin formed with the start of the reporter gene
    rbs_linker_seq = toehold_sequence[(rbs_linker_start-3):-3]
    rbs_linker_structure = toehold_structure[(rbs_linker_start-3):-3]
                
    # Use NUPACK to calculate deltaG
    prepare_nupack_input(rbs_linker_seq, os.path.join(toehold_folder, 'RBS_linker.in'))
    
    # Run NUPACK with each one
    nupack_mfe(os.path.join(toehold_folder, 'RBS_linker'))
    
    # Parse NUPACK output
    rbs_linker_result = parse_nupack(os.path.join(toehold_folder, 'RBS_linker.mfe'))
    
    rbs_linker_mfe = rbs_linker_result['energy']
    
    toehold_data.append(rbs_linker_mfe)
    
    ### End calculation for RBS linker ###
    
    # Calculate the ddG of binding (bound - unbound)
    ddG_binding = binding_energy - (toehold_alone_binding_energy + mRNA_binding_energy)
    toehold_data.append(ddG_binding)
    
    # Add this row to the list of results
    results.append(toehold_data)
    
# Save as a dataframe #### Remove references no non-paired positions
results_df = pd.DataFrame(np.array(results), 
                          columns=['Non_paired_count', 'Structure', 'Sequence', 'Start', 'End',
                                   'Length_unpaired_trigger', 'Length_paired_trigger',
                                   'Count_GC_weak', 'Count_GC_strong', 'Index', 'dG_toehold_mRNA',
                                   'Percentage_correct_matches', 'dG_toehold',
                                   'dG_mRNA', 'GC_content', 'dG_RBS_linker', 
                                   'ddG_binding'])


results_df['dG_RBS_linker'] = pd.to_numeric(results_df['dG_RBS_linker'])
results_df['ddG_binding'] = pd.to_numeric(results_df['ddG_binding'])


sorted_results = results_df.sort_values(['dG_RBS_linker', 'ddG_binding'], ascending = [False, True])
sorted_results.to_csv(path_or_buf=os.path.join(output_folder, 'all_toeholds_results.txt'), sep = '\t', index = False)

#### Add call to the script that aligns to the genomes ####
part1 = 'python alignments.py'
part2 = ' -w ' + output_folder
part3 = ' -r ' + reference_list
part4 = ' -p ' + str(pct_ident)
part5 = ' -e ' + str(evalue)
os.system(part1 + part2 + part3 + part4 + part5)

