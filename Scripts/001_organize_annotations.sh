#!/bin/bash
# -*- coding: utf-8 -*-

###############################################################################
#### This script processes the sequences from Tromer et al. Open Biol 2016 ####
#### to prepare the folders for domain annotation                          ####
###############################################################################

start_dir=$PWD

# Extract the list of sequences for BUB1
grep '^>' ../Data/Sequences/rsob160315supp3.fasta | grep 'BUB' | grep -v 'MAD' | cut -c 2- > ../Data/Sequences/Tromer_BUB1_seqs_list.txt

# Extract the list of sequences for BUBR1
grep '^>' ../Data/Sequences/rsob160315supp3.fasta | grep 'MAD' | grep -v 'BUB' | cut -c 2- > ../Data/Sequences/Tromer_BuBR1_seqs_list.txt

#  Extract the domains for BUB1 and create folders for each of them
grep -v 'MAD' ../Data/Sequences/rsob160315supp4.fasta | grep -o '||.*\/' | sort | uniq | cut -d '/' -f 1 | cut -c 3- | xargs -I {} mkdir -p ../Data/Domain_annotation/BUB1/{}

# Extract the domains for BUBR1 and create folders for each of them
grep -v 'BUB' ../Data/Sequences/rsob160315supp4.fasta | grep -o '||.*\/' | sort | uniq | cut -d '/' -f 1 | cut -c 3- | xargs -I {} mkdir -p ../Data/Domain_annotation/BuBR1/{}

# Extract all the domain annotations for BUB1 domains to fill in the folders
cd ../Data/Domain_annotation/BUB1/
for domain in $(ls)
do
	cd $domain
	echo Extracting sequences of BUB1, $domain
	grep -v 'MAD' ../../../Sequences/rsob160315supp4.fasta | grep '^>' | grep $domain | cut -c 2- > list_TromerBUB1_${domain}_ids.txt
	python ../../../../Scripts/002_extract_sequences.py -l list_TromerBUB1_${domain}_ids.txt -f ../../../Sequences/rsob160315supp4.fasta -o TromerBUB1_${domain}.fasta
	cd ..
	echo --------
done

cd $start_dir

# Extract all the domain annotations for BuBR1 domains to fill in the folders
cd ../Data/Domain_annotation/BuBR1/
for domain in $(ls)
do
        cd $domain
        echo Extracting sequences of BuBR1, $domain
        grep -v 'BUB' ../../../Sequences/rsob160315supp4.fasta | grep '^>' | grep $domain | cut -c 2- > list_TromerBuBR1_${domain}_ids.txt
        python ../../../../Scripts/002_extract_sequences.py -l list_TromerBuBR1_${domain}_ids.txt -f ../../../Sequences/rsob160315supp4.fasta -o TromerBuBR1_${domain}.fasta
        cd ..
        echo --------
done