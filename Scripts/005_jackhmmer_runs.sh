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
        jackhmmer --tblout EnsemblBUB1_${domain}_jackhmmer.tblout --domtblout EnsemblBUB1_${domain}_jackhmmer.domtblout --qformat fasta --tformat fasta TromerBUB1_${domain}.fasta ../../../Sequences/Human_BUB1_orthologues_2018_12_05.fa > EnsemblBUB1_${domain}_jackhmmer.log
        echo ------
        cd ..

done

# Run jackhmmer with the Tromer BUB1 sequences to have the same format

for domain in $(ls)
do
        echo Running jackhmmer with BUB1, ${domain}
        cd $domain
        jackhmmer --tblout TromerBUB1_${domain}_jackhmmer.tblout --domtblout TromerBUB1_${domain}_jackhmmer.domtblout --qformat fasta --tformat fasta TromerBUB1_${domain}.fasta ../../../Sequences/rsob160315supp3.fasta > TromerBUB1_${domain}_jackhmmer.log
        echo ------
        cd ..

done

# Concatenate all the files to have the final output
cat <(head -n 3 kinase/EnsemblBUB1_kinase_jackhmmer.domtblout) <(cat */*BUB1_*_jackhmmer.domtblout | grep -v '^#' ) > all_Ensembl_Tromer_BUB1_domains_jackhmmer.domtblout

cd $start_dir

#### BuBR1 ####

# Move to the main folder for the BuBR1 domain annotation
cd ../Data/Domain_annotation/BuBR1/

# Run jackhmmer with the Ensembl Orthologs BuBR1 sequences
for domain in $(ls)
do
        echo Running jackhmmer with BuBR1, ${domain}
        cd $domain
        jackhmmer --tblout EnsemblBuBR1_${domain}_jackhmmer.tblout --domtblout EnsemblBuBR1_${domain}_jackhmmer.domtblout --qformat fasta --tformat fasta TromerBuBR1_${domain}.fasta ../../../Sequences/Human_BUB1B_orthologues_2018_12_05.fa > EnsemblBuBR1_${domain}_jackhmmer.log
        echo ------
        cd ..

done

# Run jackhmmer with the Tromer BuBR1 sequences to have the same format

for domain in $(ls)
do
        echo Running jackhmmer with BuBR1, ${domain}
        cd $domain
        jackhmmer --tblout TromerBuBR1_${domain}_jackhmmer.tblout --domtblout TromerBuBR1_${domain}_jackhmmer.domtblout --qformat fasta --tformat fasta TromerBuBR1_${domain}.fasta ../../../Sequences/rsob160315supp3.fasta > TromerBuBR1_${domain}_jackhmmer.log
        echo ------
        cd ..

done

# Concatenate all the files to have the final output
cat <(head -n 3 kinase/EnsemblBuBR1_kinase_jackhmmer.domtblout) <(cat */*BuBR1_*_jackhmmer.domtblout | grep -v '^#' ) > all_Ensembl_Tromer_BuBR1_domains_jackhmmer.domtblout


