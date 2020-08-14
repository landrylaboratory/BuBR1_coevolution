#! /bin/bash

#### 010_align_filtered_sequences ####

#### This script aligns the filtered BUB1 and BUBR1 sequences to build the phylogenetic tree. ####

# Create a folder for the tree files
mkdir -p ../Results/Species_tree_files

# Retrieve the sequences
# BUBR1
python3 002_extract_sequences.py -l ../Results/MAFFT/BuBR1_IDs.txt -f ../Data/Sequences/Human_BUB1B_orthologues_2020-04-27.fa -o ../Results/Species_tree_files/BuBR1_common_seqs.fasta

# BUB1
python3 002_extract_sequences.py -l ../Results/MAFFT/BUB1_IDs.txt -f ../Data/Sequences/Human_BUB1_orthologues_2020-04-27.fa -o ../Results/Species_tree_files/BUB1_common_seqs.fasta

# MADBUB
cat > ../Results/Species_tree_files/MADBUB_IDs.txt << EOF
LAMPREY_MADBUB
BFLO_MADBUB
CINT_MADBUB
SPUR_MADBUB
SMIM_MADBUB
EOF

python3 002_extract_sequences.py -l ../Results/Species_tree_files/MADBUB_IDs.txt -f ../Data/Sequences/rsob160315supp3.fasta -o ../Results/Species_tree_files/MADBUB_outgroups.fasta

# Concatenate the sequences in a single file
cat ../Results/Species_tree_files/BUB1_common_seqs.fasta ../Results/Species_tree_files/BuBR1_common_seqs.fasta ../Results/Species_tree_files/MADBUB_outgroups.fasta > ../Results/Species_tree_files/Tree_sequences.fasta

# Align the sequences
mafft --ep 0 --genafpair --maxiterate 1000 ../Results/Species_tree_files/Tree_sequences.fasta > ../Results/Species_tree_files/Tree_sequences_aln_mafft.fasta
