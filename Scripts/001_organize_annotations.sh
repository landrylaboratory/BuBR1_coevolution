#!/bin/bash
# -*- coding: utf-8 -*-

###############################################################################
#### This script processes the sequences from Tromer et al. Open Biol 2016 ####
#### to prepare the folders for domain annotation                          ####
###############################################################################

start_dir=$PWD

#### BUB1 ####

# Extract the list of sequences for BUB1
grep '^>' ../Data/Sequences/rsob160315supp3.fasta | grep 'BUB' | grep -v 'MAD' | cut -c 2- > ../Data/Sequences/Tromer_BUB1_seqs_list.txt

#  Extract the domains for BUB1 and create folders for each of them
grep -v 'MAD' ../Data/Sequences/Tromer_domain_annotation_plk1.fasta | grep -o '||.*\/' | sort | uniq | cut -d '/' -f 1 | cut -c 3- | xargs -I {} mkdir -p ../Data/Domain_annotation/BUB1/{}

# Extract all the domain annotations for BUB1 domains to fill in the folders
cd ../Data/Domain_annotation/BUB1/
for domain in $(ls)
do
	cd $domain
	echo Extracting sequences of BUB1, $domain
	grep -v 'MAD' ../../../Sequences/Tromer_domain_annotation_plk1.fasta | grep '^>' | grep $domain | cut -c 2- > list_TromerBUB1_${domain}_ids.txt
	python3 ../../../../Scripts/002_extract_sequences.py -l list_TromerBUB1_${domain}_ids.txt -f ../../../Sequences/Tromer_domain_annotation_plk1.fasta -o TromerBUB1_${domain}.fasta
	cd ..
	echo --------
done

cd $start_dir

#### BUBR1 ####

# Extract the list of sequences for BUBR1
grep '^>' ../Data/Sequences/rsob160315supp3.fasta | grep 'MAD' | grep -v 'BUB' | cut -c 2- > ../Data/Sequences/Tromer_BuBR1_seqs_list.txt

# Extract the domains for BUBR1 and create folders for each of them
grep -v 'BUB' ../Data/Sequences/Tromer_domain_annotation_plk1.fasta | grep -o '||.*\/' | sort | uniq | cut -d '/' -f 1 | cut -c 3- | xargs -I {} mkdir -p ../Data/Domain_annotation/BuBR1/{}

# Extract all the domain annotations for BuBR1 domains to fill in the folders
cd ../Data/Domain_annotation/BuBR1/
for domain in $(ls)
do
        cd $domain
        echo Extracting sequences of BuBR1, $domain
        grep -v 'BUB' ../../../Sequences/Tromer_domain_annotation_plk1.fasta | grep '^>' | grep $domain | cut -c 2- > list_TromerBuBR1_${domain}_ids.txt
        python3 ../../../../Scripts/002_extract_sequences.py -l list_TromerBuBR1_${domain}_ids.txt -f ../../../Sequences/Tromer_domain_annotation_plk1.fasta -o TromerBuBR1_${domain}.fasta
        cd ..
        echo --------
done

cd $start_dir

#### MADBUB ####

grep '^>' ../Data/Sequences/rsob160315supp3.fasta | grep 'MADBUB' | cut -c 2- > ../Data/Sequences/Tromer_MADBUB_seqs_list.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Tromer_MADBUB_seqs_list.txt -f ../Data/Sequences/rsob160315supp3.fasta -o ../Data/Sequences/Tromer_MADBUB_seqs.fasta

grep -E '(MADBUB|PLKBD)' ../Data/Sequences/Tromer_domain_annotation_plk1.fasta | grep -o '||.*\/' | sort | uniq | cut -d '/' -f 1 | cut -c 3- | xargs -I {} mkdir -p ../Data/Domain_annotation/MADBUB/{}

cd ../Data/Domain_annotation/MADBUB/
for domain in $(ls)
do
        cd $domain
        echo Extracting sequences of MADBUB, $domain
        grep -E '(MADBUB|PLKBD)' ../../../Sequences/Tromer_domain_annotation_plk1.fasta | grep '^>' | grep $domain | cut -c 2- > list_TromerMADBUB_${domain}_ids.txt
        python3 ../../../../Scripts/002_extract_sequences.py -l list_TromerMADBUB_${domain}_ids.txt -f ../../../Sequences/Tromer_domain_annotation_plk1.fasta -o TromerMADBUB_${domain}.fasta
        cd ..
        echo --------
done

