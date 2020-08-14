#!/usr/bin/env python
# coding: utf-8

#####################################################
####		007_filter_sequences		 ####
#### This script refines the alignments of BUB1  ####
#### and BUBR1 sequences. To do this, it removes ####
#### sequences that have over 20% gaps in the	 ####
#### alignment.					 ####
#####################################################


from prody import *
from matplotlib.pylab import *
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os
import matplotlib.colors as mcol
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import OrderedDict

# Write a function to save the refined alignments
def write_prody_aln(ref_alignment, outfile):
    '''This function receives a refined alignment variable (ProDy's MSA class) and a path to an output file.
    It will save the alignment in the fasta format.
    '''
    # Start the output variable
    records = []
    
    id_list = ref_alignment.getLabels()
    
    # Loop through each of the records
    for i in range(len(id_list)):
        new_id = id_list[i]
        full_sequence = ''
        seq_array = ref_alignment.getArray()[i]
        for j in range(len(seq_array)):
            full_sequence = full_sequence + seq_array[j].decode('UTF-8')
    
        new_record = SeqRecord(Seq(full_sequence), id = new_id)
        records.append(new_record)
    # Write the sequences
    SeqIO.write(records, outfile, 'fasta')

# Make output directories
os.makedirs('../Results', exist_ok=True)
os.makedirs('../Results/Final_refined_alignments', exist_ok=True)
os.makedirs('../Results/MAFFT', exist_ok = True)
os.makedirs('../Results/MAFFT/Final_refined_alignments', exist_ok = True)
os.makedirs('../Data/Evol_results', exist_ok = True)
os.makedirs('../Data/Evol_results/MAFFT', exist_ok = True)


# ## Load BUBR1 alignment

msa_bubr1 = parseMSA('../Data/MAFFT_Alignments/Human_BUB1B_orthologues_2020-04-27_mafft.fasta')


# ## Load BUB1 alignment


msa_bub1 = parseMSA('../Data/MAFFT_Alignments/Human_BUB1_orthologues_2020-04-27_mafft.fasta')

# ## Filter so that my BUB1 and BUBR1 alignments have sequences from the same species

msa_bub1_refined = refineMSA(msa_bub1, label = 'ENSP00000302530', rowocc=0.80)
msa_bubr1_refined = refineMSA(msa_bubr1, label = 'ENSP00000287598', rowocc=0.80)

# Save the refined alignments
write_prody_aln(msa_bubr1_refined, '../Results/BUBR1_refined_aln_mafft.fasta')
write_prody_aln(msa_bub1_refined, '../Results/BUB1_refined_aln_mafft.fasta')

# Filter for organisms that pass the threshold and only have one sequence

# BuBR1
bubr1_sequences = OrderedDict()
for i in range(msa_bubr1_refined.numSequences()):
    # Reformat the labels to have a consistent formatting
    sequence_id = msa_bubr1_refined[i].getLabel()
    
    # Use regular expressions to get the prefix for each sequence
    matches = re.search(pattern = '(.*)P[0-9]+', string = sequence_id)
    new_sequence_id = matches.group(1)
    
    if bubr1_sequences.get(new_sequence_id, -1) == -1:
        bubr1_sequences[new_sequence_id] = msa_bubr1_refined[i].getLabel()
    else:
        # Remove cases of multiple duplications
        bubr1_sequences.pop(new_sequence_id)

# BUB1
bub1_sequences = OrderedDict()
for i in range(msa_bub1_refined.numSequences()):
    # Reformat the labels to have a consistent formatting
    sequence_id = msa_bub1_refined[i].getLabel()
    
    # Use regular expressions to get the prefix for each sequence
    matches = re.search(pattern = '(.*)P[0-9]+', string = sequence_id)
    new_sequence_id = matches.group(1)
    
    if bub1_sequences.get(new_sequence_id, -1) == -1:
        bub1_sequences[new_sequence_id] = msa_bub1_refined[i].getLabel()
    else:
        # Remove cases of multiple duplications
        bub1_sequences.pop(new_sequence_id)

# Make new dictionaries for the sequences that are present in both the BUB1 and the BuBR1 alignments
new_bubr1_dict = OrderedDict()
new_bub1_dict = OrderedDict()
for entry, value in bubr1_sequences.items():
    if entry in bub1_sequences.keys():
        # Then it is in both and I should save them
        new_bubr1_dict[entry] = value
        new_bub1_dict[entry] = bub1_sequences[entry]

print(len(new_bubr1_dict))
print(len(new_bub1_dict))


# Save the common sequences from each set to a file to extract them from the MSA
outfile_BuBR1 = open('../Results/MAFFT/BuBR1_IDs.txt', 'w')

for entry in new_bubr1_dict.values():
    outfile_BuBR1.write(entry + '\n')
outfile_BuBR1.close()

outfile_BUB1 = open('../Results/MAFFT/BUB1_IDs.txt', 'w')

for entry in new_bub1_dict.values():
    outfile_BUB1.write(entry + '\n')
outfile_BUB1.close()

# Extract the sequences from the MSA
# BUBR1
part0 = 'python3 002_extract_sequences.py'
part1 = ' -l ../Results/MAFFT/BuBR1_IDs.txt'
part2 = ' -f ../Data/MAFFT_Alignments/Human_BUB1B_orthologues_2020-04-27_mafft.fasta'
part3 = ' -o ../Data/Evol_results/MAFFT/BuBR1_common_seqs.fasta'
os.system(part0 + part1 + part2 + part3)

# BUB1
part0 = 'python3 002_extract_sequences.py'
part1 = ' -l ../Results/MAFFT/BUB1_IDs.txt'
part2 = ' -f ../Data/MAFFT_Alignments/Human_BUB1_orthologues_2020-04-27_mafft.fasta'
part3 = ' -o ../Data/Evol_results/MAFFT/BUB1_common_seqs.fasta'
os.system(part0 + part1 + part2 + part3)


# ## Write the final alignments for the filtered sequences and extract data on entropy

msa_bubr1 = parseMSA('../Data/Evol_results/MAFFT/BuBR1_common_seqs.fasta')
msa_bubr1_refined = refineMSA(msa_bubr1, label = 'ENSP00000287598', rowocc=0.80)
write_prody_aln(msa_bubr1_refined, '../Results/MAFFT/Final_refined_alignments/BUBR1_refined_aln_mafft.fasta')

# Save entropy data
entropy = calcShannonEntropy(msa_bubr1_refined)
np.savetxt('../Data/Evol_results/MAFFT/BUBR1_entropy_2020_05_08.txt', entropy)


# Repeat for the BUB1 data

msa_bub1 = parseMSA('../Data/Evol_results/MAFFT/BUB1_common_seqs.fasta')
msa_bub1_refined = refineMSA(msa_bub1, label = 'ENSP00000302530', rowocc=0.80)
write_prody_aln(msa_bub1_refined, '../Results/MAFFT/Final_refined_alignments/BUB1_refined_aln_mafft.fasta')

# Save entropy data
entropy = calcShannonEntropy(msa_bub1_refined)
np.savetxt('../Data/Evol_results/MAFFT/BUB1_entropy_2020_05_08.txt', entropy)


