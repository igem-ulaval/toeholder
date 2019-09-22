# toehold_design
Team iGEM ULaval 2019 scripts for toehold design

## Dependencies

The following programs must be installed and added to the PATH:
- NUPACK: available at http://www.nupack.org/downloads
- BLAST+: available at https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

## Scripts

All scripts are written in Python 3 and depend on the following libraries

- master_script.py: Manages the other scripts to design toeholds for a particular target gene.

- generate_toehold.py: Sweeps through the sequence of the target gene looking for suitable candidate recognition sequences. All candidate recognition sequences are evaluated based on the following parameters:
	. Secondary structure on the mRNA
	. ddG of the bound (toehold + target) and unbound state (toehold and target, separately)

- toehold_functions.py: Contains several helper functions for the other scripts.

- specificity_check.py: Aligns the toeholds to a selected set of reference genomes to identify matches. Toeholds matching more than one sequence or matching sequences from other genomes would not be completely specific to the target.

## Examples

In order to run the examples, reference genomes must be downloaded from: 
- E_coli_DH5alpha: https://www.ncbi.nlm.nih.gov/nuccore/CP026085.1
- Human: https://www.ncbi.nlm.nih.gov/nuccore/?term=human+genome
