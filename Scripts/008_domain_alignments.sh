#!/bin/bash
# -*- coding: utf-8 -*-

###############################################################################
#### This script performs the multiple domain alignments with MDAT using   ####
#### the full sets of sequences and the best-hit domain annotations.       ####
###############################################################################

mkdir ../Data/MDAT_Alignments

mdat -s ../Data/Sequences/Ensembl_Tromer_all_BUB1_no_duplicates.fa \
-d ../Data/Domain_annotation/BUB1/domain_table_BUB1_Pfam_scan_format.txt \
-o ../Data/MDAT_Alignments/BUB1_mdat_alignment.fasta

mdat -s ../Data/Sequences/Ensembl_Tromer_all_BuBR1_no_duplicates.fa \
-d ../Data/Domain_annotation/BuBR1/domain_table_BuBR1_Pfam_scan_format.txt \
-o ../Data/MDAT_Alignments/BuBR1_mdat_alignment.fasta