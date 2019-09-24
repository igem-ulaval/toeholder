# toehold_design
Team iGEM ULaval 2019 scripts for toehold design

## Dependencies

The following programs must be installed and added to the PATH:
- NUPACK (Zadeh et al. 2011. Journal of Computational Chemistry): available at http://www.nupack.org/downloads
- BLAST+ (Camacho et al. 2009. BMC Bioinformatics): available at https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

BLAST+ can also be installed with the following command line:

```
sudo apt-get install ncbi-blast+
```

It also depends on the following Python libraries:
- Biopython (Cock et al. 2009. Bioinformatics)
- Numpy
- Pandas

These libraries can be installed with the following commands:

```
pip3 install biopython
pip3 install numpy
pip3 install pandas
```

## Scripts

All scripts are written in Python 3 and depend on the following libraries

- toehold_generator.py: Sweeps through the sequence of the target gene looking for suitable candidate recognition sequences. All candidate recognition sequences are evaluated based on the following parameters:
	- Secondary structure on the mRNA
	- ddG of the bound (toehold + target) and unbound state (toehold and target, separately)

- input_variables.py: Defines tunable parameters and input files for the toehold_generator script.

- toehold_helper_functions.py: Contains several helper functions for the other scripts.

- alignments.py: Aligns the toeholds to a selected set of reference genomes to identify matches. Toeholds matching more than one sequence or matching sequences from other genomes would not be completely specific to the target.

## Input

The workflow is fully automated to be executed as follows as long as all the scripts are in the current directory:

```
python toehold_generator.py
```

The input_variables.py script contains all the necessary adjustable variables for the script, including:
- Input sequence (FASTA formatted)
- Length of the unpaired region of the recognition sequence (α in the diagram)
- Length of the paired region of the recognition sequence (β in the diagram)
- Path to the output folder
- Reporter gene or tag to be added to the end of the toehold
- Input sequence molecule type (DNA or RNA)
- List of reference genomes formatted as follows:

path_to_genome_1	Tag_1
path_to_genome_2	Tag_2
path_to_genome_3	Tag_3

- Percentage identity threshold for hits to retain from the alignments
- Evalue threshold for hits to retain from the alignments
- Minimum number of unpaired bases in the secondary structure of the target mRNA for a candidate trigger to be considered

![](Figures/toehold_diagram.png)

## Output

The toehold_generator.py script generates an output folder with a subfolder for each of the candidate toeholds generated. When the candidate toehold contains a stop codon, its corresponding subfolder is empty. When it does not contain a stop codon, there are four files inside the subfolder:
- switch1_python.in: NUPACK-formatted input to test the toehold's secondary structure.
- switch1_python.mfe: NUPACK-formatted output with the most favorable secondary structure.
- toehold_mRNA.in: NUPACK-formatted input to test the toehold's ability to bind to the target mRNA.
- toehold_mRNA.mfe: NUPACK-formatted output with the most favorable structure of the toehold-mRNA complex.

Outside those folders, the rest of the files are:
- input_variables.py: Copy of the input variables used for this run.
- mRNA.in: NUPACK-formatted input to test the mRNA's secondary structure.
- mRNA.mfe: NUPACK-formatted output with the most favorable secondary structure for the mRNA.
- toehold_candidates.txt: First list of toeholds ranked by the number of non-paired positions in the secondary structure of the trigger region of the target mRNA.
- toehold_seqs.fasta: FASTA-formatted file containing the recognition sequences for each of the toeholds.
- all_toeholds_results.txt: Results of the tests performed on the toeholds. It adds the following columns to the toehold_candidates.txt file:
	- Toehold index
	- Free energy of the toehold-mRNA complex
	- Percentage of paired bases of the toehold-mRNA complex that correspond the intended base pairing
	- Free energy of the toehold secondary structure
	- Free energy of the mRNA secondary structure
	- GC content of the toehold recognition sequence
- \<tag\>_toeholds_alignment.aln: Output of the BLAST alignment of the library of toeholds to the genome referenced with the corresponding tag in the genome list.
- all_toeholds_results_genome_matches.txt: Adds the counts of matches for each toehold in each of the genomes referenced in the genome list.



