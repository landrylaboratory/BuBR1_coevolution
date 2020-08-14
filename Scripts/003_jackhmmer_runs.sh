#!/bin/bash
# -*- coding: utf-8 -*-

###############################################################################
#### This script executes the jackhmmer runs and merges their output to    ####
#### have matching coordinates for each domain.                            ####
###############################################################################

start_dir=$PWD

#### BUB1 ####

# Move to the main folder for the BUB1 domain annotation
cd ../Data/Domain_annotation/BUB1/

# Run jackhmmer with the Ensembl Orthologs BUB1 sequences
for domain in $(ls)
do
        echo Running jackhmmer with BUB1, ${domain}
        cd $domain
        jackhmmer --tblout EnsemblBUB1_${domain}_jackhmmer.tblout --domtblout EnsemblBUB1_${domain}_jackhmmer.domtblout --qformat fasta --tformat fasta TromerBUB1_${domain}.fasta ../../../Sequences/Human_BUB1_orthologues_2020-04-27.fa > EnsemblBUB1_${domain}_jackhmmer.log
        echo ------
        cd ..

done

# Concatenate all the files to have the final output
cat <(head -n 3 kinase/EnsemblBUB1_kinase_jackhmmer.domtblout) <(cat */*BUB1_*_jackhmmer.domtblout | grep -v '^#' ) > all_Ensembl_BUB1_domains_jackhmmer.domtblout

cd $start_dir

#### BuBR1 ####

# Move to the main folder for the BuBR1 domain annotation
cd ../Data/Domain_annotation/BuBR1/

# Run jackhmmer with the Ensembl Orthologs BuBR1 sequences
for domain in $(ls)
do
        echo Running jackhmmer with BuBR1, ${domain}
        cd $domain
        jackhmmer --tblout EnsemblBuBR1_${domain}_jackhmmer.tblout --domtblout EnsemblBuBR1_${domain}_jackhmmer.domtblout --qformat fasta --tformat fasta TromerBuBR1_${domain}.fasta ../../../Sequences/Human_BUB1B_orthologues_2020-04-27.fa > EnsemblBuBR1_${domain}_jackhmmer.log
        echo ------
        cd ..

done

# Concatenate all the files to have the final output
cat <(head -n 3 kinase/EnsemblBuBR1_kinase_jackhmmer.domtblout) <(cat */*BuBR1_*_jackhmmer.domtblout | grep -v '^#' ) > all_Ensembl_BuBR1_domains_jackhmmer.domtblout

cd $start_dir

#### MADBUB ####

# Move to the main folder for the MADBUB domain annotation
cd ../Data/Domain_annotation/MADBUB/

# Run jackhmmer with the Tromer MADBUB sequences
for domain in $(ls)
do
        echo Running jackhmmer with MADBUB, ${domain}
        cd $domain
        jackhmmer --tblout TromerMADBUB_${domain}_jackhmmer.tblout --domtblout TromerMADBUB_${domain}_jackhmmer.domtblout --qformat fasta --tformat fasta TromerMADBUB_${domain}.fasta ../../../Sequences/Tromer_MADBUB_seqs.fasta > TromerMADBUB_${domain}_jackhmmer.log
        echo ------
        cd ..

done

cat <(head -n 3 kinase/TromerMADBUB_kinase_jackhmmer.domtblout) <(cat */*MADBUB_*_jackhmmer.domtblout | grep -v '^#' ) > all_Tromer_MADBUB_domains_jackhmmer.domtblout



