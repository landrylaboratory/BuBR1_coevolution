#!/bin/bash
# -*- coding: utf-8 -*-

###############################################################################
#### This script aligns the sequences from Tromer et al. Open Biol 2016    ####
#### to the sequences from Ensembl Orthologs to create the complete data   ####
#### set without duplicated sequences.                                     ####
###############################################################################

start_dir=$PWD

# Prepare the BLAST DB files for the alignments
cd ../Data/Sequences/
makeblastdb -in Tromer_BUB1_seqs.fasta -dbtype prot
makeblastdb -in Tromer_BuBR1_seqs.fasta -dbtype prot

# Run BLASTP to get the alignments
mkdir ../BLAST_Alignments
blastp -query Human_BUB1_orthologues_2018_12_05.fa -db Tromer_BUB1_seqs.fasta -evalue 0.000001 -outfmt 6 > ../BLAST_Alignments/Ensembl_orthologs_vs_Tromer_BUB1.aln
blastp -query Human_BUB1B_orthologues_2018_12_05.fa -db Tromer_BuBR1_seqs.fasta -evalue 0.000001 -outfmt 6 > ../BLAST_Alignments/Ensembl_orthologs_vs_Tromer_BuBR1.aln

# Look for the duplicate sequences (above 95% sequence identity and belonging to the same organism)
awk '$3 + 0 >= 95' Ensembl_orthologs_vs_Tromer_BUB1.aln | cut -f 2 
awk '$3 + 0 >= 95' Ensembl_orthologs_vs_Tromer_BuBR1.aln | cut -f 2

# Retrieve for the duplicate sequences (98.5% identity allows removal of the ones that are annotated to the same organism)
awk '$3 + 0 >= 98.5' Ensembl_orthologs_vs_Tromer_BUB1.aln | cut -f 2 > duplicate_sequences_BUB1.txt
awk '$3 + 0 >= 98.5' Ensembl_orthologs_vs_Tromer_BuBR1.aln | cut -f 2 > duplicate_sequences_BuBR1.txt

# Merge the two datasets prior to removing the duplicates
cat ../Sequences/Human_BUB1_orthologues_2018_12_05.fa ../Sequences/Tromer_BUB1_seqs.fasta > all_Ensembl_Tromer_BUB1_seqs.fasta
cat ../Sequences/Human_BUB1B_orthologues_2018_12_05.fa ../Sequences/Tromer_BuBR1_seqs.fasta > all_Ensembl_Tromer_BuBR1_seqs.fasta
