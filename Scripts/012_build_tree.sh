#!/bin/bash

#### 012_build_tree ####

# This script uses the alignment of the filtered BUB1, BUBR1, and MADBUB sequences to build a phylogenetic tree with MrBayes
# Before running the script: Convert the alignment from the fasta format to a Nexus format. This can be done with the following web tool: http://phylogeny.lirmm.fr/phylo_cgi/data_converter.cgi

# File to convert: ../Results/Species_tree_files/Tree_sequences_aln_mafft_newIDs.fasta
# Save the converted file as: ../Results/Species_tree_files/Tree_sequences_aln_mafft_newIDs.nex

cd ../Results/Species_tree_files

# Open MrBayes
mb

# Load file
execute Tree_sequences_aln_mafft_newIDs.nex 

# Set evolutionary model for proteins
prset aamodelpr=mixed 	

# Set outgroup
outgroup Strongylocentrotus_purpuratus_MADBUB 

# Start the sampling
mcmc ngen=60000 samplefreq=10 

# Summarise results
sump burnin=50 
sumt burnin=50


