Scripts for "The human BUBR1 pseudokinase domain is required for efficient Spindle Checkpoint silencing and chromosome alignment via allosteric activation of KARD phosphorylation" by Gama-Braga, L., Cisneros, A. F., Mathieu, M., Clerc, M., Garcia, P., Lottin, B., Garand, C., Thebault, P., Landry, C. R., & Elowe, S.

#########################################
#	Dependencies			#
#########################################

The analyses performed by these scripts depend on the following programs and libraries:

- Evol (from the ProDy Python library, version 1.10.11, Bakan et al. Bioinformatics 2014): Available at: http://prody.csb.pitt.edu/evol/

- BLASTP (version 2.6.0+, Camacho et al. BMC Bioinformatics 2009): Available at: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

- MAFFT (version 7.450, Katoh, et al. Nucleic Acids Research 2002; Katoh, et al. Molecular Biology and Evolution 2013): Available at: https://mafft.cbrc.jp/alignment/software/

- Jackhmmer (from the HMMER suite, version 3.1b2, Johnson, et al. BMC Bioinformatics 2010): Available at: http://hmmer.org/download.html

- MrBayes (version 3.2.7a,  Huelsenbeck & Ronquist. Bioinformatics 2001; Ronquist & Huelsenbeck. Bioinformatics 2003): Available at: https://nbisweden.github.io/MrBayes/download.html

- UpSetR (version 1.4.0, Lex et al. IEEE Transctions on Visualization and Computer Graphics, 2014): Available at: https://cran.r-project.org/web/packages/UpSetR/index.html

- ggtree (version 2.0.4, Yu. Current Protocols in Bioinformatics, 2020): Available at: https://bioconductor.org/packages/release/bioc/html/ggtree.html

#########################################
#	Data folder			#
#########################################

Contains the starting files for the analysis on domain co-evolution, as well as all the files created by the scripts.

Starting with the Sequences subfolder and the following files, the scripts generate the rest of the files and folders:

- Human_BUB1B_orthologues_2020_04_27.fa: The list of orthologs for the human BUB1B protein as downloaded from Ensembl Orthologs on April 27th, 2020.
- Human_BUB1_orthologues_2020_04_27.fa: The list of orthologs for the human BUB1 protein as downloaded from Ensembl Orthologs on April 27th, 2020.
- rsob160315supp3.fasta: Supplementary materials containing sequences used in Tromer et al. Open Biol. 2016, available at http://dx.doi.org/10.1098/rsob.160315.
- rsob160315supp4.fasta: Supplementary materials containing domain annotations curated in Tromer et al. Open Biol. 2016, available at http://dx.doi.org/10.1098/rsob.160315.
- Tromer_domain_annotation_plk1.fasta: A copy of the rsob160315supp4.fasta file that contains the annotation for the PLK binding domain for both hBUB1 and hBUBR1.

#########################################
#	Figures folder			#
#########################################

Contains all the figures produced for the paper by the scripts.

#########################################
#	Scripts folder			#
#########################################

Contains the scripts used to perform the analyses on domain co-evolution for BUB1 and BuBR1. Unless otherwise stated, they can be run without arguments. For the R scripts, it is necessary to set the working directory as indicated within each script.

- 001_organize_annotations.sh
Produces a folder organization that will be used to store the domain annotations, as well as lists of sequences to separate MADBUB, BUB1, and BuBR1 sequences from Tromer's published lists.

- 002_extract_sequences.py
Extracts the BUB1 and BuBR1 sequences in a FASTA formatted file. It is a helper script that is automatically called by the other scripts when needed.

- 003_jackhmmer_runs.sh
Runs jackhmmer to have the matching regions of each protein for each domain.

- 004_annotate_domains.py
This script looks at the merged output from the jackhmmer runs and uses it to obtain an annotation with the best hits and a table indicating which domains are present in each sequence.

- 005_domain_copresence.R
This script looks at table of presence of domains in each sequence and produces figures on domain co-presence in BUBR1 and BUB1, respectively.

- 006_MAFFT_alignments.sh
Performs the multiple domain alignments with MAFFT. Stores the results in the Data/MAFFT_Alignments folder. It also aligns the sequences of domains/motifs of hBUB1 to those of hBUBR1.

- 007_filter_sequences.py
This script uses the Evol suite from the Prody package to filter the sequences and refine the multiple sequence alignment. The alignment is refined by removing sequences that have more than 20% gaps and taking the human sequence as a reference. Finally, the filter makes sure that only sequences from vertebrate species that have exactly one copy of BUB1 and one copy of BUBR1 are retained. The refined alignment is then used to calculate entropy for each position as a metric of sequence variation.

- 008_entropy_table.py
This script uses the domain alignments of hBUB1 and hBUBR1 as well as the table with entropy values. It outputs a table that indicates which positions of the hBUB1 and the hBUBR1 are matched in the alignment with their corresponding entropy values.

- 009_entropy_analysis.py
This is the script reads the table with the entropy values and draws a figure to visualize it.

- 010_align_filtered_sequences.R
This script takes the BUB1 and BUBR1 sequences from the refined alignments as well as 5 MADBUB sequences (outgroups) and realigns them. This new alignment will be used to generate a phylogeny of these sequences that were used for the entropy analysis.

- 011_rename_sequences.py
This script replaces the sequence IDs with their corresponding species and whether it is the BUB1 or the BUBR1 protein.

- 012_build_tree.sh 
This script contains the series of commands used with MrBayes to build the phylogenetic tree. It also explains the conversion procedure used to convert the alignment from a fasta format to a nexus format.

- 013_tree_visualization.R
This script reads the phylogenetic tree and produces the figures.



