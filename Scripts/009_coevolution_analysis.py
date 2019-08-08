#!/usr/bin/env python
# coding: utf-8

# # Coevolution Evol MDAT
# 
# This is the script that looks at coevolution between domains based on variation observed in the multiple domain alignments.


from prody import *
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib.colors as mcol
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import OrderedDict


# ## 1.- Look at how the number of BuBR1 sequences decays when filtering by occupancy
# Set working directory to the scripts folder
# cd <path_to_scripts_folder> 


# Load the alignment
msa_bubr1 = parseMSA('../Data/MDAT_Alignments/BuBR1_mdat_alignment.fasta')


# Use the refine function with the human sequence and different occupancy thresholds
checks = []

for i in range(0,100,5):
    threshold = float(i)/100
    msa_bubr1_refined = refineMSA(msa_bubr1, label = 'ENSP00000287598', rowocc=threshold)
    checks.append(len(msa_bubr1_refined))

y_pos = range(0,100,5)

msa_bubr1_refined = refineMSA(msa_bubr1, label = 'ENSP00000287598', rowocc=0.80)
occupancy_data = pd.DataFrame(np.column_stack((y_pos, checks)), columns=['Occupancy', 'Sequences'])
occupancy_data.to_csv(path_or_buf=os.path.join('../Data/Evol_results/', 'BuBR1_occupancy_filter.txt'), sep = '\t', index = False)


# ## 2.- Look at how the number of BUB1 sequences decays when filtering by occupancy
# Load the alignment
msa_bub1 = parseMSA('../Data/MDAT_Alignments/BUB1_mdat_alignment.fasta')


# Use the refine function with the human sequence and different occupancy thresholds
checks = []

for i in range(0,100,5):
    threshold = float(i)/100
    msa_bub1_refined = refineMSA(msa_bub1, label = 'ENSP00000302530', rowocc=threshold)
    checks.append(len(msa_bub1_refined))

y_pos = range(0,100,5)

msa_bub1_refined = refineMSA(msa_bub1, label = 'ENSP00000302530', rowocc=0.80)
occupancy_data = pd.DataFrame(np.column_stack((y_pos, checks)), columns=['Occupancy', 'Sequences'])
occupancy_data.to_csv(path_or_buf=os.path.join('../Data/Evol_results/', 'BUB1_occupancy_filter.txt'), sep = '\t', index = False)

# ## 3.- Look at the BuBR1 and BUB1 sequences to filter for those that belong to the same organism

# BuBR1
bubr1_sequences = OrderedDict()
for i in range(msa_bubr1_refined.numSequences()):
    # Reformat the labels to have a consistent formatting
    sequence_id = msa_bubr1_refined[i].getLabel()
    
    # First case: Tromer's sequences
    if sequence_id.endswith('MAD'):
        new_sequence_id = sequence_id[0:4]
    # Second case: Ensembl sequences (not mouse genome project)
    elif sequence_id.startswith('ENS'):
        new_sequence_id = sequence_id[0:7]
    # Third case: Mouse Genome Project sequences
    elif sequence_id.startswith('MGP'):
        new_sequence_id = sequence_id.split('_')[1]
    
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
    
    # First case: Tromer's sequences
    if sequence_id.endswith('BUB'):
        new_sequence_id = sequence_id[0:4]
    # Second case: Ensembl sequences (not mouse genome project)
    elif sequence_id.startswith('ENS'):
        new_sequence_id = sequence_id[0:7]
    # Third case: Mouse Genome Project sequences
    elif sequence_id.startswith('MGP'):
        new_sequence_id = sequence_id.split('_')[1]
    
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
        # Save sequences that are present in both
        new_bubr1_dict[entry] = value
        new_bub1_dict[entry] = bub1_sequences[entry]

# Save the common sequences from each set to a file to extract them from the MSA
outfile_BuBR1 = open('../Data/Evol_results/BuBR1_IDs.txt', 'w')
for entry in new_bubr1_dict.values():
    outfile_BuBR1.write(entry + '\n')
outfile_BuBR1.close()

outfile_BUB1 = open('../Data/Evol_results/BUB1_IDs.txt', 'w')
for entry in new_bub1_dict.values():
    outfile_BUB1.write(entry + '\n')
outfile_BUB1.close()

# Extract the sequences from the MSA
# BUBR1
part0 = 'python 002_extract_sequences.py'
part1 = ' -l ../Data/Evol_results/BuBR1_IDs.txt'
part2 = ' -f ../Data/MDAT_Alignments/BuBR1_mdat_alignment.fasta'
part3 = ' -o ../Data/Evol_results/BuBR1_common_seqs.fasta'
os.system(part0 + part1 + part2 + part3)

# BUB1
part0 = 'python 002_extract_sequences.py'
part1 = ' -l ../Data/Evol_results/BUB1_IDs.txt'
part2 = ' -f ../Data/MDAT_Alignments/BUB1_mdat_alignment.fasta'
part3 = ' -o ../Data/Evol_results/BUB1_common_seqs.fasta'
os.system(part0 + part1 + part2 + part3)


# ## 4.- Perform the coevolution analyses with the extracted sequences
msa_bubr1 = parseMSA('../Data/Evol_results/BuBR1_common_seqs.fasta')
msa_bubr1_refined = refineMSA(msa_bubr1, label = 'ENSP00000287598', rowocc=0.80)

indices = range(len(msa_bubr1_refined[0]))
entropy = calcShannonEntropy(msa_bubr1_refined)

mutinfo = buildMutinfoMatrix(msa_bubr1_refined)
mutinfo_corr = applyMutinfoCorr(mutinfo, corr='apc')
np.savetxt('../Data/Evol_results/BuBR1_common_seqs_APC.txt', mutinfo_corr)
np.savetxt('../Data/Evol_results/BuBR1_common_seqs_entropy.txt', entropy)

msa_bub1 = parseMSA('../Data/Evol_results/BUB1_common_seqs.fasta')

msa_bub1_refined = refineMSA(msa_bub1, label = 'ENSP00000302530', rowocc=0.80)

indices = range(len(msa_bub1_refined[0]))
entropy = calcShannonEntropy(msa_bub1_refined)

mutinfo = buildMutinfoMatrix(msa_bub1_refined)
mutinfo_corr = applyMutinfoCorr(mutinfo, corr='apc')
np.savetxt('../Data/Evol_results/BUB1_common_seqs_APC.txt', mutinfo_corr)
np.savetxt('../Data/Evol_results/BUB1_common_seqs_entropy.txt', entropy)


# ## 5.- Merge the sequences from the same organism to look at coevolution between domains of the different proteins.

# Parse the sequences and save them to dictionaries
bub1_seqs = SeqIO.parse('../Data/Evol_results/BUB1_common_seqs.fasta', 'fasta')
bubr1_seqs = SeqIO.parse('../Data/Evol_results/BuBR1_common_seqs.fasta', 'fasta')

bub1_dict = OrderedDict()
bubr1_dict = OrderedDict()

for record in bub1_seqs:
    # I could reformat the labels to have a consistent formatting
    sequence_id = record.id
    # First case: Tromer's sequences
    if sequence_id.endswith('BUB'):
        new_sequence_id = sequence_id[0:4]
    # Second case: Ensembl sequences (not mouse genome project)
    elif sequence_id.startswith('ENS'):
        new_sequence_id = sequence_id[0:7]
    # Third case: Mouse Genome Project sequences
    elif sequence_id.startswith('MGP'):
        new_sequence_id = sequence_id.split('_')[1]
    
    bub1_dict[new_sequence_id] = record
    
for record in bubr1_seqs:
    # I could reformat the labels to have a consistent formatting
    sequence_id = record.id
    # First case: Tromer's sequences
    if sequence_id.endswith('MAD'):
        new_sequence_id = sequence_id[0:4]
    # Second case: Ensembl sequences (not mouse genome project)
    elif sequence_id.startswith('ENS'):
        new_sequence_id = sequence_id[0:7]
    # Third case: Mouse Genome Project sequences
    elif sequence_id.startswith('MGP'):
        new_sequence_id = sequence_id.split('_')[1]
    
    bubr1_dict[new_sequence_id] = record


merged_records = []

for key in bub1_dict.keys():
    # Get the corresponding two records
    record_bub1 = bub1_dict[key]
    record_bubr1 = bubr1_dict[key]
    
    # Prepare the record for the merged sequences
    new_sequence = Seq(str(record_bub1.seq) + str(record_bubr1.seq))
    new_record = SeqRecord(id=record_bub1.id + '-' + record_bubr1.id, seq=new_sequence, description = '')

    merged_records.append(new_record)
    
# Save the records to an output fasta file
SeqIO.write(merged_records, '../Data/Evol_results/merged_BUB1-BuBR1_common_seqs.fasta', 'fasta')

# ## 6.- Repeat the coevolution analyses with the merged sequences
msa_merged = parseMSA('../Data/Evol_results/merged_BUB1-BuBR1_common_seqs.fasta')

msa_merged_refined = refineMSA(msa_merged, label = 'ENSP00000302530-ENSP00000287598', rowocc=0.80)

indices = range(len(msa_merged_refined[0]))
entropy = calcShannonEntropy(msa_merged_refined)

mutinfo = buildMutinfoMatrix(msa_merged_refined)
mutinfo_corr = applyMutinfoCorr(mutinfo, corr='apc')
np.savetxt('/media/axelle/Angel_backup/Dropbox/BuBR1/BuBR1_coevolution/Data/Evol_results/merged_common_seqs_APC.txt', mutinfo_corr)
np.savetxt('/media/axelle/Angel_backup/Dropbox/BuBR1/BuBR1_coevolution/Data/Evol_results/merged_common_seqs_entropy.txt', entropy)



