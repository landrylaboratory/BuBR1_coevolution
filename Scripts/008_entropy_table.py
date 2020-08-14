#!/usr/bin/env python
# coding: utf-8

#### 008_entropy_table ####
# 
# This script will look at the alignments of the domains the hBUB1 and hBUBR1 sequences have in common. It will use them to match the entropy values to return a table that will simplify the analysis.

# Load libraries
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
import os

# Create the output directory
os.makedirs('../Data/Sequences/Human_domains/Alignments/Tables', exist_ok = True)

# Load the entropy data
entropy_bub1_file = open('../Data/Evol_results/MAFFT/BUB1_entropy_2020_05_08.txt', 'r')
entropy_bub1 = [entry for entry in csv.reader(entropy_bub1_file, delimiter = '\t')]

entropy_bubr1_file = open('../Data/Evol_results/MAFFT/BUBR1_entropy_2020_05_08.txt', 'r')
entropy_bubr1 = [entry for entry in csv.reader(entropy_bubr1_file, delimiter = '\t')]


# Loop through the domains to build the table containing entropy values and the positions of hBUB1 and hBUBR1 that are matched in the alignment

domains = ['TPR', 'GLEBS', 'ABBA3', 'PLKBD', 'KARD', 'CDII', 'kinase']

outfile_path = os.path.join('../Data/Sequences/Human_domains/Alignments/Tables/', 'Human_domains_Tromer_table.txt')
outfile = open(outfile_path, 'w')

writer = csv.writer(outfile, delimiter = '\t')
header = ['Domain', 'Protein', 'Position', 'Residue', 'Group']
writer.writerow(header)

accumulated = 0

for domain in domains:

    print('Working with', domain)
    
    aln_file = os.path.join('../Data/Sequences/Human_domains/Alignments/', 'Human_' + domain + '_Tromer_sequences_aln.fasta')
    domain_aln = SeqIO.parse(aln_file, 'fasta')

    for entry in domain_aln:
        # Check the ID to get the data
        if('BUB' in entry.id):
            protein = 'BUB1'

        elif('MAD' in entry.id):
            protein = 'BUBR1'

        other_info = entry.id.split('||')[1].split('/')
        positions = [int(pos) for pos in other_info[1].split('-')]

        # Adjust numbering after the TPR to account for the 14-residue insertion in Tromer's sequence
        for i in range(len(positions)):        
            if domain != 'TPR' and positions[i] > 93 and protein == 'BUBR1':
                positions[i] = positions[i] - 14

        # Adjust the numbering of the BUB1 kinase, which was offset by 13 residues
        if domain == 'kinase' and protein == 'BUB1':
            for i in range(len(positions)):
                positions[i] = positions[i] + 13
                
        # Match the domains to their positions in the entropy vector
        # Use groups to pair the matching residues in the alignment
        curr_group = 1
        curr_pos = positions[0]
        for residue in entry.seq:
            # Remove the 14-residue insertion that is only present in Tromer's BUBR1 sequence
            if domain == 'TPR' and (curr_group >= 21 and curr_group <= 34):
                curr_group = curr_group + 1
                continue

            if residue == '-':
                # Outfile columns will be domain, protein, position, residue, group
                new_line = [domain, protein, -1, '-', curr_group  + accumulated]           
            else:
                new_line = [domain, protein, curr_pos, residue, curr_group  + accumulated]
                curr_pos = curr_pos + 1
            curr_group = curr_group + 1
            writer.writerow(new_line)

    accumulated = accumulated + curr_group -1
    print('------')
outfile.close()




