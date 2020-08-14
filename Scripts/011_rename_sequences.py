#!/usr/bin/env python
# coding: utf-8

#### 011_rename_sequences ####

# This script will take the file with the alignments and rename the sequence IDs. The new IDs will indicate the species name and whether it is a BUB1 (BUB), BUBR1 (MAD), or a MADBUB (MADBUB) protein.

# Load libraries
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import OrderedDict
import csv
import re

# Load files
fasta_file = SeqIO.parse('../Results/Species_tree_files/Tree_sequences_aln_mafft.fasta', 'fasta')

handle = open('../Data/Ensembl_species_names_no_dups.txt', 'r')
id_table = csv.reader(handle, delimiter = '\t')

# Load the lists of BUB1 and BUBR1 sequence IDs to be able to distinguish them
BUBR1_handle = open('../Results/MAFFT/BuBR1_IDs.txt', 'r')
BUBR1_IDs = csv.reader(BUBR1_handle, delimiter = '\t')

BUB1_handle = open('../Results/MAFFT/BUB1_IDs.txt', 'r')
BUB1_IDs = csv.reader(BUB1_handle, delimiter = '\t')

# Save the BUB1 and BUBR1 IDs as a list
BUB1_list = []
for entry in BUB1_IDs:
    BUB1_list.append(entry[0])

BUBR1_list = []
for entry in BUBR1_IDs:
    BUBR1_list.append(entry[0])

# Fill the dictionary for the IDs
id_dict = OrderedDict()

header = next(id_table)

for entry in id_table:
    original_id = entry[0]
    species_id = entry[1]
    
    # Remove the common names and use underscores instead of spaces
    new_id = species_id.split(' (')[0].replace(' ', '_')
    
    id_dict[original_id] = new_id

# Loop through the sequence IDs
outfile = '../Results/Species_tree_files/Tree_sequences_aln_mafft_newIDs.fasta'
records = []

for record in fasta_file:
    
    record.description = ''
    
    # For MADBUB sequences, I can just rename and add "_MADBUB"
    if record.id.endswith('MADBUB'):
        record.id = id_dict[record.id] + '_MADBUB'
    # Work with BUB1 and BUBR1 sequences
    else:    

        # Extract the prefixes
        prefix = re.search('(.*)P[0-9]+', record.id).group(1)

        # Use the prefixes to map to the species name
        new_record_id = id_dict[prefix]

        # Check if this ID corresponds to a BUB1 or a BUBR1 sequence
        if record.id in BUB1_list:
            new_record_id = new_record_id + '_BUB'
        elif record.id in BUBR1_list:
            new_record_id = new_record_id + '_MAD'

        # Replace the sequence IDs
        record.id = new_record_id
    
    records.append(record)

# Write a new fasta file
SeqIO.write(records, outfile, 'fasta')



