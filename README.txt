Scripts for "The human BUBR1 pseudokinase domain is required for efficient Spindle Checkpoint silencing and chromosome alignment via allosteric activation of KARD phosphorylation" by Gama-Braga, L., Cisneros, A. F., Mathieu, M., Clerc, M., Garcia, P., Lottin, B., Garand, C., Thebault, P., Landry, C. R., & Elowe, S.

#########################################
#	Dependencies			#
#########################################

The analyses performed by these scripts depend on the following programs and libraries:

- Evol (from the ProDy Python library, version 1.10.8, Bakan et al. Bioinformatics 2014): Available at: http://prody.csb.pitt.edu/evol/

- BLASTP (version 2.6.0+, Camacho et al. BMC Bioinformatics 2009): Available at: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

- MDAT (version 1.0, Kemena et al. BMC Bioinformatics 2015): Available at: http://www.bornberglab.org/pages/mdat

- Jackhmmer (from the HMMER suite, version 3.1b2, Johnson, et al. BMC Bioinformatics 2010): Available at: http://hmmer.org/download.html

#########################################
#	Data folder			#
#########################################

Contains the starting files for the analysis on domain co-evolution, as well as all the files created by the scripts.

Starting with the Sequences subfolder and the following files, the scripts generate the rest of the files and folders:

- Human_BUB1B_orthologues_2018_12_05.fa: The list of orthologs for the human BUB1B protein as downloaded from Ensembl Orthologs on December 5th, 2018.
- Human_BUB1_orthologues_2018_12_05.fa: The list of orthologs for the human BUB1 protein as downloaded from Ensembl Orthologs on December 5th, 2018.
- rsob160315supp3.fasta: Supplementary materials containing sequences used in Tromer et al. Open Biol. 2016, available at http://dx.doi.org/10.1098/rsob.160315.
- rsob160315supp4.fasta: Supplementary materials containing domain annotations curated in Tromer et al. Open Biol. 2016, available at http://dx.doi.org/10.1098/rsob.160315.

#########################################
#	Figures folder			#
#########################################

Contains all the figures produced for the paper by the scripts.

#########################################
#	Scripts folder			#
#########################################

Contains the scripts used to perform the analyses on domain co-evolution for BUB1 and BuBR1. Unless otherwise stated, they can be run without arguments.

- 001_organize_annotations.sh
Produces a folder organization that will be used to store the domain annotations, as well as lists of sequences to separate BUB1 and BuBR1 sequences from Tromer's published lists.

- 002_extract_sequences.py
Extracts the BUB1 and BuBR1 sequences in a FASTA formatted file. It is used as follows:
> python 002_extract_sequences.py -l ../Data/Sequences/Tromer_BUB1_seqs_list.txt -f ../Data/Sequences/rsob160315supp3.fasta -o ../Data/Sequences/Tromer_BUB1_seqs.fasta
> python 002_extract_sequences.py -l ../Data/Sequences/Tromer_BuBR1_seqs_list.txt -f ../Data/Sequences/rsob160315supp3.fasta -o ../Data/Sequences/Tromer_BuBR1_seqs.fasta

- 003_align_merge_sets.sh
Runs the alignment of the sequences from both data sets to find duplicates. Stores the results in the BLAST_Alignments folder.

- 004_remove_sequences.py
Removes the duplicate sequences from the merged data set to obtain the final data sets with 173 sequences for BUB1 and 176 for BUBR1. It is used as follows:
> python 004_remove_sequences.py -l ../Data/BLAST_Alignments/duplicate_sequences_BUB1.txt -f ../Data/BLAST_Alignments/all_Ensembl_Tromer_BUB1_seqs.fasta -o ../Data/Sequences/Ensembl_Tromer_all_BUB1_no_duplicates.fa
> python 004_remove_sequences.py -l ../Data/BLAST_Alignments/duplicate_sequences_BuBR1.txt -f ../Data/BLAST_Alignments/all_Ensembl_Tromer_BuBR1_seqs.fasta -o ../Data/Sequences/Ensembl_Tromer_all_BuBR1_no_duplicates.fa

- 005_jackhmmer_runs.sh
Runs jackhmmer to have the matching regions of each protein for each domain.

- 006_annotate_domains.py
This script looks at the merged output from the jackhmmer runs and uses it to obtain an annotation with the best hits and a table indicating which domains are present in each sequence.

- 007_domain_copresence.R
This script looks at table of presence of domains in each sequence and produces figures Suppl1A and Suppl1F on domain co-presence in BUBR1 and BUB1, respectively.

- 008_domain_alignments.sh
Performs the multiple domain alignments with MDAT. Stores the results in the MDAT_Alignments folder. 

- 009_coevolution_analysis.py
This is the script that looks at co-evolution between domains based on variation observed in the multiple domain alignments. Stores the resulting files in the Evol_results folder.

- 010_coevolution_figures.R
Uses the output from the co-evolution analysis performed in the previous script and produces figures 1B, 1C, Suppl1A, Suppl1E, and Suppl6D.



