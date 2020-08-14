#!/usr/bin/env python
# coding: utf-8

# # 004_annotate_domains.py
# 
# This script will use the jackhmmer searches and save the best hits of each of the domains as the annotation.
# 


# Set current directory to the scripts folder
# cd <path_to_scripts_folder>


# Import libraries
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
import os
import numpy as np

# Define some helper functions
def parse_jackhmmer(infile):
    '''This function parses jackhmmer files into generators.'''
    with open(infile, 'r') as handle:
        for line in handle:
            # I am not interested in the lines that have comments
            if not line.startswith('#'):
                yield re.split('\s+', line)

# Define the thresholds used to accept annotations to merge them
min_coverage = 0.6
min_overlap = 0.5
shortest_domain_threshold = 0.7

# Define the paths to the BUB1 and BuBR1 folders
bub1_folder = '../Data/Domain_annotation/BUB1/'
bubr1_folder = '../Data/Domain_annotation/BuBR1/'
madbub_folder = '../Data/Domain_annotation/MADBUB/'


# ## 1.- Look at the BUBR1 sequences 

# Use the jackhmmer files to produce Pfam_scan formatted file as detailed in https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Protein+Functional+Analysis+Output+Examples.

jackhmmer_domtblout = os.path.join(bubr1_folder, 'all_Ensembl_BuBR1_domains_jackhmmer.domtblout')

domains_BuBR1 = parse_jackhmmer(jackhmmer_domtblout)

# Initialize a dictionary that will save the information about the domains
domains_BuBR1_dict = OrderedDict()

# A list that will keep track of all the domains used as a reference
domain_list = []

proteins = []

# Open output file
# handle_out = open(os.path.join(bubr1_folder, 'new_annotation', 'domain_table_BuBR1_Pfam_scan_format.txt'), 'w')
handle_out = open(os.path.join(bubr1_folder, 'domain_table_BuBR1_Pfam_scan_format.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

# Loop through all the lines in the jackhmmer output file
for entry in domains_BuBR1:
    # Save the needed data
    sequence_name = entry[0] # Column 1 of the output
    ali_start = entry[17] # Column 2 of the output
    ali_end = entry[18] # Column 3 of the output
    env_start = entry[19] # Column 4 of the output
    env_end = entry[20] # Column 5 of the output
    query_name = entry[3] # Column 6 of the output
    query_prot_name = query_name.split('||')[0] # Column 7 of the output
    query_prot_domain = query_name.split('||')[1].split('/')[0] # Column 8 of the output
    hmm_start = entry[15] # Column 9 of the output
    hmm_end = entry[16] # Column 10 of the output
    hmm_length = entry[5] # Column 11 of the output
    bitscore = entry[7] # Column 12 of the output
    evalue = entry[6] # Column 13 of the output
    post_prob = float(entry[21]) # Column 14 of the output
    clan = 'No_clan' # Column 15 of the output, not working with clans

    # Process the query name to extract the name of the domain
    domain_name = query_name.split('||')[1].split('/')[0]
    
    # Skip sequences that correspond to Tromer's BUB or MADBUB sequences
    if '_MADBUB' in sequence_name or '_BUB' in sequence_name:
        continue
        
    # Only save matches with a posterior probability of a least 0.95
    if post_prob >= 0.95:
        # Save the line to the output
        out_line = [sequence_name, ali_start, ali_end, env_start, env_end, query_name, query_prot_name, query_prot_domain, hmm_start, hmm_end, hmm_length, bitscore, evalue, post_prob, clan]
        writer.writerow(out_line)
                    
        # Add the domain to the list of domains
        if not domain_name in domain_list:
            domain_list.append(domain_name)
            
        # Add the sequence and the domain matches to the dictionary
        if not sequence_name in domains_BuBR1_dict.keys():
            domains_BuBR1_dict[sequence_name] = OrderedDict()
            domains_BuBR1_dict[sequence_name][domain_name] = 1
        else:
            domains_BuBR1_dict[sequence_name][domain_name] = 1

handle_out.close()           


# ## Select the best matches from the table as the annotations

# Save the best hits as a Pfam_scan formatted table.

# in_file = open(os.path.join(bubr1_folder, 'new_annotation', 'domain_table_BuBR1_Pfam_scan_format.txt'), 'r')
in_file = open(os.path.join(bubr1_folder, 'domain_table_BuBR1_Pfam_scan_format.txt'), 'r')
reader = csv.reader(in_file, delimiter = '\t')

domain_consensus_dict = OrderedDict()

for entry in reader:
    # Retrieve the data
    ali_start = float(entry[1])
    ali_end = float(entry[2])
    query_length = float(entry[10])
    protein = entry[0]
    score = float(entry[11]) # Bitscore
    
    query_protein = entry[6]
    query_start = entry[8]
    query_end = entry[9]
    evalue = entry[12]
    prob = entry[13]
    
    # Merge the three kinase-like domain annotations, as well as the ABBA and KEN annotations
    if entry[7] in ['kinase', 'pseudo_kinase', '_pseudo_kinase']:
        domain = '(Pseudo)kinase'
    elif entry[7] in ['ABBA1', 'ABBA2', 'ABBA_other']:
        domain = 'ABBA'
    elif entry[7] in ['KEN1', 'KEN2', 'KEN_other']:
        domain = 'KEN'
    elif entry[7] == 'Dbox':
        domain = 'D-Box'
    else:
        domain = entry[7]
    
    # Get the length of the alignment (1-based positions)
    ali_length = ali_end - ali_start + 1
    
    # Check if the length of the alignment has a good coverage of the query
    if ali_length/query_length >= min_coverage:
        if domain_consensus_dict.get(protein, -1) == -1:
            domain_consensus_dict[protein] = OrderedDict()
        
        # If this is the first time this domain appears, save it
        if domain_consensus_dict[protein].get(domain, -1) == -1:
            domain_consensus_dict[protein][domain] = OrderedDict()
            domain_consensus_dict[protein][domain][1] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
        # If this domain has already appeared, look at the overlap with all the instances
        else:
            new_domain = True
            for domain_num, position_list in domain_consensus_dict[protein][domain].items():
                # Calculate overlap as intersection over union (the percentage of the union that
                # is comprised by the intersection)
                current_start = position_list[0]
                current_end = position_list[1]
                current_best_score = position_list[2]
                         
                intersection = float(min(ali_end, current_end) - max(ali_start, current_start)) + 1
                union = float(max(ali_end, current_end) - min(ali_start, current_start)) + 1
                overlap = intersection / union
                
                # Check how much of the shorter sequence is in the intersection
                short_seq_pct = intersection / (min((ali_end - ali_start), (current_end - current_start)) + 1)

                # If the overlap threshold is cleared or if the intersection contains most of the 
                # shorter sequence, keep only the one with the higher score
                if overlap >= min_overlap or short_seq_pct >= shortest_domain_threshold:
                    # Then look at the bitscore to decide which to maintain
                    if score > current_best_score:
                        domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
                    new_domain = False
                    break
                # else:
                elif intersection > 0:
                    # Otherwise, discard the proposed domain annotation
                    print('Discarded domain: ', protein, domain)
                    print('New observation:', ali_start, ali_end)
                    print('Current annotation:', current_start, current_end)
                    print(intersection)
                    print('------')
                    new_domain = False
 
            # Instances of domains that did not overlap with previous ones will be added as new instances
            if new_domain:
                domain_num = max(domain_consensus_dict[protein][domain].keys()) + 1
                domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
                
in_file.close()                 


# Need to make sure in the domain_consensus_dict all the domains are sorted

temp_dict = OrderedDict()

# Loop through each protein's domains
for protein, domain_dict in domain_consensus_dict.items():
    temp_dict[protein] = OrderedDict()
    
    for domain, instance_dict in domain_dict.items():
        temp_dict[protein][domain] = OrderedDict()
        instance_list = []
        # Build the list of entries
        for num, entry in instance_dict.items():
            entry[0] = float(entry[0])
            instance_list.append(entry)
            
        # Use the list of entries to build a numpy table
        domain_matrix = np.array(instance_list)
        
        # Sort the matrix by the start position
        sorted_matrix = domain_matrix[domain_matrix[:,0].astype(np.float).argsort()]
        
        # Loop through the rows to assign them in the dictionary
        num_rows, num_cols = sorted_matrix.shape
        for i in range(num_rows):
            temp_dict[protein][domain][i + 1] = sorted_matrix[i,:].tolist()


# Replace the domain_consensus_dict with the temp_dict
domain_consensus_dict = temp_dict


# Save the best hits as a Pfam_scan formatted table.

# outfile = open(os.path.join(bubr1_folder, 'Final_annotation', 'BuBR1_best_hit_domain_annotation_table_new.txt'), 'w')
outfile = open(os.path.join(bubr1_folder, 'BuBR1_best_hit_domain_annotation_table_new.txt'), 'w')
writer = csv.writer(outfile, delimiter = '\t')

for protein, domain_dict in domain_consensus_dict.items():
    for domain, instance_list in domain_dict.items():
        
        for domain_num, positions in instance_list.items():
            domain_start = positions[0]
            domain_end = positions[1]
            score = positions[2]
            query_prot = positions[3]
            query_start = positions[4]
            query_end = positions[5]
            query_length = int(float(positions[6]))
            evalue = positions[7]
            prob = positions[8]
            
            new_row = [protein, domain_start, domain_end, domain_start, domain_end, 'Best_hit', query_prot, domain + '_' + str(domain_num), query_start, query_end, query_length, score, evalue, prob, 'No_clan']
            writer.writerow(new_row)

outfile.close()


# Save the best hits with their sequences as a FASTA formatted file.


records = []

sequences = SeqIO.parse('../Data/Sequences/Human_BUB1B_orthologues_2020-04-27.fa', 'fasta')

# Save a dictionary to simplify access to the sequences
seq_dict = OrderedDict()
for seq in sequences:
    seq_dict[seq.id] = seq.seq
    
for protein, domain_dict in domain_consensus_dict.items():
    for domain, instance_list in domain_dict.items():
        for domain_num, positions in instance_list.items():
            domain_start = int(float(positions[0]))
            domain_end = int(float(positions[1]))
            
            # Subtract 1 to the start to change from 1-based notation to 0-based notation
            sequence = seq_dict[protein][domain_start-1:domain_end]
            
            # Produce the output in Tromer's notation
            new_id = protein + '||' + domain + '_' + str(domain_num) + '/' + str(domain_start) + '-' + str(domain_end)
            
            record = SeqRecord(sequence, id = new_id, description = '')
            
            records.append(record)
            
SeqIO.write(records, os.path.join(bubr1_folder, 'BuBR1_best_hit_domain_annotation.fasta'), 'fasta')

# Write a table that indicates the positions of each domain

# Prepare a new copy of the dictionary with the positions of the best hit domains.
# infile = open(os.path.join(bubr1_folder, 'Final_annotation', 'BuBR1_best_hit_domain_annotation_table_new.txt'), 'r')
infile = open(os.path.join(bubr1_folder, 'BuBR1_best_hit_domain_annotation_table_new.txt'), 'r')
reader = csv.reader(infile, delimiter = '\t')

bubr1_position_dict = OrderedDict()
domain_list = []

for entry in reader:
    sequence_id = entry[0]
    domain_start = int(float(entry[1]))
    domain_end = int(float(entry[2]))
    domain_id = entry[7]
    
    if bubr1_position_dict.get(sequence_id, -1) == -1:
        bubr1_position_dict[sequence_id] = OrderedDict()
        
    if not domain_id in domain_list:
        domain_list.append(domain_id)
        
    bubr1_position_dict[sequence_id][domain_id] = str(domain_start) + '-' + str(domain_end)

infile.close()

# Open output file
# handle_out = open(os.path.join(bubr1_folder, 'Final_annotation', 'domain_table_BuBR1_new.txt'), 'w')
handle_out = open(os.path.join(bubr1_folder, 'domain_table_BuBR1_new.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

# Write headers
writer.writerow(['Sequence'] + domain_list)

# Look at each of the sequences
for sequence in bubr1_position_dict.keys():
    new_line = [sequence]
    
    # Loop through the list of domains to get the table on which sequences have which domains
    for domain in domain_list:
        # Get a 1 if that domain was found in that sequence or a 0 otherwise
        new_line.append(bubr1_position_dict[sequence].get(domain, 0))
        
    # Write the new line
    writer.writerow(new_line)
    
handle_out.close()


# The following command saves a discretized table of presence (1) or absence (0) of each domain in each sequence.

# handle_in = open(os.path.join(bubr1_folder, 'Final_annotation', 'domain_table_BuBR1_new.txt'), 'r')
handle_in = open(os.path.join(bubr1_folder, 'domain_table_BuBR1_new.txt'), 'r')
reader = csv.reader(handle_in, delimiter = '\t')

# handle_out = open(os.path.join(bubr1_folder, 'Final_annotation', 'domain_table_BuBR1_discrete_new.txt'), 'w')
handle_out = open(os.path.join(bubr1_folder, 'domain_table_BuBR1_discrete_new.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

header = next(reader)
writer.writerow(header)

# Loop through the entries
for entry in reader:
    # If an entry contains positions for a domain change them to a 1
    for i in range(1, len(entry)):
        if entry[i] != '0':
            entry[i] = 1
    
    writer.writerow(entry)

handle_out.close()


# ## 2.- Work with the BUB1 sequences.

jackhmmer_domtblout = os.path.join(bub1_folder, 'all_Ensembl_BUB1_domains_jackhmmer.domtblout')

domains_BUB1 = parse_jackhmmer(jackhmmer_domtblout)

# Initialize a dictionary that will save the information about the domains
domains_BUB1_dict = OrderedDict()

# A list that will keep track of all the domains used as a reference
domain_list = []

# Open output file
# handle_out = open(os.path.join(bub1_folder, 'new_annotation', 'domain_table_BUB1_Pfam_scan_format.txt'), 'w')
handle_out = open(os.path.join(bub1_folder, 'domain_table_BUB1_Pfam_scan_format.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

proteins = []

# Loop through all the lines in the jackhmmer output file
for entry in domains_BUB1:
    # Save the data I am interested in
    sequence_name = entry[0] # Column 1 of the output
    ali_start = entry[17] # Column 2 of the output
    ali_end = entry[18] # Column 3 of the output
    env_start = entry[19] # Column 4 of the output
    env_end = entry[20] # Column 5 of the output
    query_name = entry[3] # Column 6 of the output
    query_prot_name = query_name.split('||')[0] # Column 7 of the output
    query_prot_domain = query_name.split('||')[1].split('/')[0] # Column 8 of the output
    hmm_start = entry[15] # Column 9 of the output
    hmm_end = entry[16] # Column 10 of the output
    hmm_length = entry[5] # Column 11 of the output
    bitscore = entry[7] # Column 12 of the output
    evalue = entry[6] # Column 13 of the output
    post_prob = float(entry[21]) # Column 14 of the output
    clan = 'No_clan' # Column 15 of the output, since I am not working with clans

    # I can process the query name to extract the name of the domain
    domain_name = query_name.split('||')[1].split('/')[0]
    
    # Skip sequences that correspond to Tromer's BUB or MADBUB sequences
    if '_MADBUB' in sequence_name or '_MAD' in sequence_name:
        continue
    
    # Only save matches with a posterior probability of a least 0.95
    if post_prob >= 0.95:
        # Save the line to the output
        out_line = [sequence_name, ali_start, ali_end, env_start, env_end, query_name, query_prot_name, query_prot_domain, hmm_start, hmm_end, hmm_length, bitscore, evalue, post_prob, clan]
        writer.writerow(out_line)
                    
        # Add the domain to the list of domains
        if not domain_name in domain_list:
            domain_list.append(domain_name)
            
        # Add the sequence and the domain matches to the dictionary
        if not sequence_name in domains_BUB1_dict.keys():
            domains_BUB1_dict[sequence_name] = OrderedDict()
            domains_BUB1_dict[sequence_name][domain_name] = 1
        else:
            domains_BUB1_dict[sequence_name][domain_name] = 1

handle_out.close()

# ## Select the best matches from the table as the annotations

# in_file = open(os.path.join(bub1_folder, 'new_annotation', 'domain_table_BUB1_Pfam_scan_format.txt'), 'r')
in_file = open(os.path.join(bub1_folder, 'domain_table_BUB1_Pfam_scan_format.txt'), 'r')
reader = csv.reader(in_file, delimiter = '\t')

domain_consensus_dict = OrderedDict()

for entry in reader:
    # Retrieve the needed data
    ali_start = float(entry[1])
    ali_end = float(entry[2])
    query_length = float(entry[10])
    protein = entry[0]
    score = float(entry[11]) # Bitscore
    
    query_protein = entry[6]
    query_start = entry[8]
    query_end = entry[9]
    evalue = entry[12]
    prob = entry[13]
    
    # Merge the three kinase-like domain annotations, as well as the ABBA and KEN annotations
    if entry[7] in ['kinase', 'pseudo_kinase', '_pseudo_kinase']:
        domain = 'kinase'
    elif entry[7] in ['ABBA1', 'ABBA2', 'ABBA_other']:
        domain = 'ABBA'
    elif entry[7] in ['KEN1', 'KEN2', 'KEN_other']:
        domain = 'KEN'
    elif entry[7] == 'Dbox':
        domain = 'D-Box'
    else:
        domain = entry[7]
    
    # Get the length of the alignment (1-based positions)
    ali_length = ali_end - ali_start + 1
    
    # Check if the length of the alignment has a good coverage of the query
    if ali_length/query_length >= min_coverage:
        if domain_consensus_dict.get(protein, -1) == -1:
            domain_consensus_dict[protein] = OrderedDict()
        
        # If this is the first time this domain appears, save it
        if domain_consensus_dict[protein].get(domain, -1) == -1:
            domain_consensus_dict[protein][domain] = OrderedDict()
            if domain == 'ABBA':
                # ABBA for BUB1 should be initialized as ABBA3 because it does not have the first two
                domain_consensus_dict[protein][domain][3] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
            else:
                domain_consensus_dict[protein][domain][1] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]

        # If this domain has already appeared, look at the overlap with all the instances
        else:
            new_domain = True
            for domain_num, position_list in domain_consensus_dict[protein][domain].items():
                # Calculate overlap as intersection over union (the percentage of the union that
                # is comprised by the intersection)
                current_start = position_list[0]
                current_end = position_list[1]
                         
                intersection = float(min(ali_end, current_end) - max(ali_start, current_start)) + 1
                union = float(max(ali_end, current_end) - min(ali_start, current_start)) + 1
                overlap = intersection / union

                # Check how much of the shorter sequence is in the intersection
                short_seq_pct = intersection / (min((ali_end - ali_start), (current_end - current_start)) + 1)

                # If the overlap threshold is cleared or if the intersection contains most of the 
                # shorter sequence, keep only the one with the higher scorer
                if overlap >= min_overlap or short_seq_pct >= shortest_domain_threshold:
                    # Then look at the bitscore to decide which to maintain
                    if score > current_best_score:
                        domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
                    new_domain = False
                    break
                # else:
                elif intersection > 0:
                    # Otherwise, discard the proposed domain annotation
                    print('Discarded domain: ', protein, domain)
                    print('New observation:', ali_start, ali_end)
                    print('Current annotation:', current_start, current_end)
                    print(intersection, check)
                    print('------')
                    new_domain = False
 
            # Instances of domains that did not overlap with previous ones will be added as new instances
            if new_domain:
                domain_num = max(domain_consensus_dict[protein][domain].keys()) + 1
                domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]

in_file.close()    


# Need to make sure in the domain_consensus_dict all the domains are sorted

temp_dict = OrderedDict()

# Loop through each protein's domains
for protein, domain_dict in domain_consensus_dict.items():
    temp_dict[protein] = OrderedDict()
    
    for domain, instance_dict in domain_dict.items():
        temp_dict[protein][domain] = OrderedDict()
        instance_list = []
        # Build the list of entries
        for num, entry in instance_dict.items():
            entry[0] = float(entry[0])
            instance_list.append(entry)
            
        # Use the list of entries to build a numpy table
        domain_matrix = np.array(instance_list)
        
        # Sort the matrix by the start position
        sorted_matrix = domain_matrix[domain_matrix[:,0].astype(np.float).argsort()]
        
        # Loop through the rows to assign them in the dictionary
        num_rows, num_cols = sorted_matrix.shape
        for i in range(num_rows):
            # ABBA domains for BUB1 should start with ABBA3
            if domain == 'ABBA':
                temp_dict[protein][domain][i + 3] = sorted_matrix[i,:].tolist()
            else:
                temp_dict[protein][domain][i + 1] = sorted_matrix[i,:].tolist()

# Replace the domain_consensus_dict with the temp_dict
domain_consensus_dict = temp_dict

# Save the best hits as a Pfam_scan formatted table.
# outfile = open(os.path.join(bub1_folder, 'Final_annotation', 'BUB1_best_hit_domain_annotation_table_new.txt'), 'w')
outfile = open(os.path.join(bub1_folder, 'BUB1_best_hit_domain_annotation_table_new.txt'), 'w')
writer = csv.writer(outfile, delimiter = '\t')

for protein, domain_dict in domain_consensus_dict.items():
    
    for domain, instance_list in domain_dict.items():
        for domain_num, positions in instance_list.items():
            domain_start = positions[0]
            domain_end = positions[1]
            score = positions[2]
            query_prot = positions[3]
            query_start = positions[4]
            query_end = positions[5]
            query_length = int(float(positions[6]))
            evalue = positions[7]
            prob = positions[8]
            
            new_row = [protein, domain_start, domain_end, domain_start, domain_end, 'Best_hit', query_prot, domain + '_' + str(domain_num), query_start, query_end, query_length, score, evalue, prob, 'No_clan']
            writer.writerow(new_row)
            
outfile.close()


# Save the best hits with their sequences as a FASTA formatted file.

records = []

sequences = SeqIO.parse('../Data/Sequences/Human_BUB1_orthologues_2020-04-27.fa', 'fasta')

# Save a dictionary to simplify access to the sequences
seq_dict = OrderedDict()
for seq in sequences:
    seq_dict[seq.id] = seq.seq
    
for protein, domain_dict in domain_consensus_dict.items():
    for domain, instance_list in domain_dict.items():
        for domain_num, positions in instance_list.items():
            domain_start = int(float(positions[0]))
            domain_end = int(float(positions[1]))
            
            # Subtract 1 to the start to change from 1-based notation to 0-based notation
            sequence = seq_dict[protein][domain_start-1:domain_end]
            
            # Produce the output in Tromer's notation
            new_id = protein + '||' + domain + '_' + str(domain_num) + '/' + str(domain_start) + '-' + str(domain_end)
            
            record = SeqRecord(sequence, id = new_id, description = '')
            
            records.append(record)
            
# SeqIO.write(records, os.path.join(bub1_folder, 'Final_annotation', 'BUB1_best_hit_domain_annotation.fasta'), 'fasta')
SeqIO.write(records, os.path.join(bub1_folder, 'BUB1_best_hit_domain_annotation.fasta'), 'fasta')

## Write a table with the positions of each domain

# Prepare a new copy of the dictionary with the positions of the best hit domains.
# infile = open(os.path.join(bub1_folder, 'Final_annotation', 'BUB1_best_hit_domain_annotation_table_new.txt'), 'r')
infile = open(os.path.join(bub1_folder, 'BUB1_best_hit_domain_annotation_table_new.txt'), 'r')
reader = csv.reader(infile, delimiter = '\t')

bub1_position_dict = OrderedDict()
domain_list = []

for entry in reader:
    sequence_id = entry[0]
    domain_start = int(float(entry[1]))
    domain_end = int(float(entry[2]))
    domain_id = entry[7]
    
    if bub1_position_dict.get(sequence_id, -1) == -1:
        bub1_position_dict[sequence_id] = OrderedDict()
        
    if not domain_id in domain_list:
        domain_list.append(domain_id)
        
    bub1_position_dict[sequence_id][domain_id] = str(domain_start) + '-' + str(domain_end)

infile.close()

# Open output file
# handle_out = open(os.path.join(bub1_folder, 'Final_annotation', 'domain_table_BUB1_new.txt'), 'w')
handle_out = open(os.path.join(bub1_folder, 'domain_table_BUB1_new.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

# Write headers
writer.writerow(['Sequence'] + domain_list)

# Look at each of the sequences
for sequence in bub1_position_dict.keys():
    new_line = [sequence]
    
    # Loop through the list of domains to get the table on which sequences have which domains
    for domain in domain_list:
        # Get a 1 if that domain was found in that sequence or a 0 otherwise
        new_line.append(bub1_position_dict[sequence].get(domain, 0))
        
    # Write the new line
    writer.writerow(new_line)
    
handle_out.close()

# Write a table with 0 for absence of a domain or 1 for presence

# handle_in = open(os.path.join(bub1_folder, 'Final_annotation', 'domain_table_BUB1_new.txt'), 'r')
handle_in = open(os.path.join(bub1_folder, 'domain_table_BUB1_new.txt'), 'r')
reader = csv.reader(handle_in, delimiter = '\t')

# handle_out = open(os.path.join(bub1_folder, 'Final_annotation', 'domain_table_BUB1_discrete_new.txt'), 'w')
handle_out = open(os.path.join(bub1_folder, 'domain_table_BUB1_discrete_new.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

header = next(reader)
writer.writerow(header)

# Loop through the entries
for entry in reader:
    # If an entry contains positions for a domain change them to a 1
    for i in range(1, len(entry)):
        if entry[i] != '0':
            entry[i] = 1
    
    writer.writerow(entry)

handle_out.close()


# ## 3.- Work with the MADBUB annotations

jackhmmer_domtblout = os.path.join(madbub_folder, 'all_Tromer_MADBUB_domains_jackhmmer.domtblout')

domains_MADBUB = parse_jackhmmer(jackhmmer_domtblout)

# Initialize a dictionary that will save the information about the domains
domains_MADBUB_dict = OrderedDict()

# A list that will keep track of all the domains used as a reference
domain_list = []

# Open output file
handle_out = open(os.path.join(madbub_folder, 'domain_table_MADBUB_Pfam_scan_format.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

proteins = []

# Loop through all the lines in the jackhmmer output file
for entry in domains_MADBUB:
    # Save the data I am interested in
    sequence_name = entry[0] # Column 1 of the output
    ali_start = entry[17] # Column 2 of the output
    ali_end = entry[18] # Column 3 of the output
    env_start = entry[19] # Column 4 of the output
    env_end = entry[20] # Column 5 of the output
    query_name = entry[3] # Column 6 of the output
    query_prot_name = query_name.split('||')[0] # Column 7 of the output
    query_prot_domain = query_name.split('||')[1].split('/')[0] # Column 8 of the output
    hmm_start = entry[15] # Column 9 of the output
    hmm_end = entry[16] # Column 10 of the output
    hmm_length = entry[5] # Column 11 of the output
    bitscore = entry[7] # Column 12 of the output
    evalue = entry[6] # Column 13 of the output
    post_prob = float(entry[21]) # Column 14 of the output
    clan = 'No_clan' # Column 15 of the output, since I am not working with clans

    # I can process the query name to extract the name of the domain
    domain_name = query_name.split('||')[1].split('/')[0]
    
    # Only save matches with a posterior probability of a least 0.95
    if post_prob >= 0.95:
        # Save the line to the output
        out_line = [sequence_name, ali_start, ali_end, env_start, env_end, query_name, query_prot_name, query_prot_domain, hmm_start, hmm_end, hmm_length, bitscore, evalue, post_prob, clan]
        writer.writerow(out_line)
                    
        # Add the domain to the list of domains
        if not domain_name in domain_list:
            domain_list.append(domain_name)
            
        # Add the sequence and the domain matches to the dictionary
        if not sequence_name in domains_MADBUB_dict.keys():
            domains_MADBUB_dict[sequence_name] = OrderedDict()
            domains_MADBUB_dict[sequence_name][domain_name] = 1
        else:
            domains_MADBUB_dict[sequence_name][domain_name] = 1

handle_out.close()


# ## Select the best matches from the table as the annotations

in_file = open(os.path.join(madbub_folder, 'domain_table_MADBUB_Pfam_scan_format.txt'), 'r')
reader = csv.reader(in_file, delimiter = '\t')

domain_consensus_dict = OrderedDict()

for entry in reader:
    # Retrieve the needed data
    ali_start = float(entry[1])
    ali_end = float(entry[2])
    query_length = float(entry[10])
    protein = entry[0]
    score = float(entry[11]) # Bitscore
    
    query_protein = entry[6]
    query_start = entry[8]
    query_end = entry[9]
    evalue = entry[12]
    prob = entry[13]
    
    # Merge the three kinase-like domain annotations, as well as the ABBA and KEN annotations
    if entry[7] in ['kinase', 'pseudo_kinase', '_pseudo_kinase']:
        domain = 'kinase'
    elif entry[7] in ['ABBA1', 'ABBA2', 'ABBA_other']:
        domain = 'ABBA'
    elif entry[7] in ['KEN1', 'KEN2', 'KEN_other']:
        domain = 'KEN'
    elif entry[7] == 'Dbox':
        domain = 'D-Box'
    else:
        domain = entry[7]
    
    # Get the length of the alignment (1-based positions)
    ali_length = ali_end - ali_start + 1
    
    # Check if the length of the alignment has a good coverage of the query
    if ali_length/query_length >= min_coverage:
        if domain_consensus_dict.get(protein, -1) == -1:
            domain_consensus_dict[protein] = OrderedDict()
        
        # If this is the first time this domain appears, save it
        if domain_consensus_dict[protein].get(domain, -1) == -1:
            domain_consensus_dict[protein][domain] = OrderedDict()
            domain_consensus_dict[protein][domain][1] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
        
        # If this domain has already appeared, look at the overlap with all the instances
        else:
            new_domain = True
            for domain_num, position_list in domain_consensus_dict[protein][domain].items():
                # Calculate overlap as intersection over union (the percentage of the union that
                # is comprised by the intersection)
                current_start = position_list[0]
                current_end = position_list[1]
                current_best_score = position_list[2]
                         
                intersection = float(min(ali_end, current_end) - max(ali_start, current_start)) + 1
                union = float(max(ali_end, current_end) - min(ali_start, current_start)) + 1
                overlap = intersection / union

                # Check how much of the shorter sequence is in the intersection
                short_seq_pct = intersection / (min((ali_end - ali_start), (current_end - current_start)) + 1)

                # If the overlap threshold is cleared or if the intersection contains most of the 
                # shorter sequence, keep only the one with the higher scorer
                if overlap >= min_overlap or short_seq_pct >= shortest_domain_threshold:
                    # Then look at the bitscore to decide which to maintain
                    if score > current_best_score:
                        domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
                    new_domain = False
                    break
                # else:
                elif intersection > 0:
                    # Otherwise, discard the proposed domain annotation
                    print('Discarded domain: ', protein, domain)
                    print('New observation:', ali_start, ali_end)
                    print('Current annotation:', current_start, current_end)
                    print(intersection)
                    print('------')
                    new_domain = False
 
            # Instances of domains that did not overlap with previous ones will be added as new instances
            if new_domain:
                domain_num = max(domain_consensus_dict[protein][domain].keys()) + 1
                domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]

in_file.close() 


# Need to make sure in the domain_consensus_dict all the domains are sorted

temp_dict = OrderedDict()

# Loop through each protein's domains
for protein, domain_dict in domain_consensus_dict.items():
    temp_dict[protein] = OrderedDict()
    
    for domain, instance_dict in domain_dict.items():
        temp_dict[protein][domain] = OrderedDict()
        instance_list = []
        # Build the list of entries
        for num, entry in instance_dict.items():
            entry[0] = float(entry[0])
            instance_list.append(entry)
            
        # Use the list of entries to build a numpy table
        domain_matrix = np.array(instance_list)
        
        # Sort the matrix by the start position
        sorted_matrix = domain_matrix[domain_matrix[:,0].astype(np.float).argsort()]
        
        # Loop through the rows to assign them in the dictionary
        num_rows, num_cols = sorted_matrix.shape
        for i in range(num_rows):
            temp_dict[protein][domain][i + 1] = sorted_matrix[i,:].tolist()

# Replace the domain_consensus_dict with the temp_dict
domain_consensus_dict = temp_dict

# Save the best hits as a Pfam_scan formatted table.

# outfile = open(os.path.join(madbub_folder, 'Final_annotation', 'MADBUB_best_hit_domain_annotation_table.txt'), 'w')
outfile = open(os.path.join(madbub_folder, 'MADBUB_best_hit_domain_annotation_table.txt'), 'w')
writer = csv.writer(outfile, delimiter = '\t')

for protein, domain_dict in domain_consensus_dict.items():
    
    for domain, instance_list in domain_dict.items():
        for domain_num, positions in instance_list.items():
            domain_start = positions[0]
            domain_end = positions[1]
            score = positions[2]
            query_prot = positions[3]
            query_start = positions[4]
            query_end = positions[5]
            query_length = int(float(positions[6]))
            evalue = positions[7]
            prob = positions[8]
            
            
            new_row = [protein, domain_start, domain_end, domain_start, domain_end, 'Best_hit', query_prot, domain + '_' + str(domain_num), query_start, query_end, query_length, score, evalue, prob, 'No_clan']
            writer.writerow(new_row)

outfile.close()

# Save the best hits with their sequences as a FASTA formatted file.

records = []

sequences = SeqIO.parse('../Data/Sequences/Tromer_MADBUB_seqs.fasta', 'fasta')

# Save a dictionary to simplify access to the sequences
seq_dict = OrderedDict()
for seq in sequences:
    seq_dict[seq.id] = seq.seq
    
for protein, domain_dict in domain_consensus_dict.items():
    for domain, instance_list in domain_dict.items():
        for domain_num, positions in instance_list.items():
            domain_start = int(float(positions[0]))
            domain_end = int(float(positions[1]))
            
            # Subtract 1 to the start to change from 1-based notation to 0-based notation
            sequence = seq_dict[protein][domain_start-1:domain_end]
            
            # Produce the output in Tromer's notation
            new_id = protein + '||' + domain + '_' + str(domain_num) + '/' + str(domain_start) + '-' + str(domain_end)
            
            record = SeqRecord(sequence, id = new_id, description = '')
            
            records.append(record)
            
# SeqIO.write(records, os.path.join(madbub_folder, 'Final_annotation', 'MADBUB_best_hit_domain_annotation.fasta'), 'fasta')
SeqIO.write(records, os.path.join(madbub_folder, 'MADBUB_best_hit_domain_annotation.fasta'), 'fasta')

# Write a table indicating the positions of the domains

# Prepare a new copy of the dictionary with the positions of the best hit domains.
# infile = open(os.path.join(madbub_folder, 'Final_annotation', 'MADBUB_best_hit_domain_annotation_table.txt'), 'r')
infile = open(os.path.join(madbub_folder, 'MADBUB_best_hit_domain_annotation_table.txt'), 'r')
reader = csv.reader(infile, delimiter = '\t')

madbub_position_dict = OrderedDict()
domain_list = []

for entry in reader:
    sequence_id = entry[0]
    domain_start = int(float(entry[1]))
    domain_end = int(float(entry[2]))
    domain_id = entry[7]
    
    if madbub_position_dict.get(sequence_id, -1) == -1:
        madbub_position_dict[sequence_id] = OrderedDict()
        
    if not domain_id in domain_list:
        domain_list.append(domain_id)
        
    madbub_position_dict[sequence_id][domain_id] = str(domain_start) + '-' + str(domain_end)

infile.close()


# Open output file
# handle_out = open(os.path.join(madbub_folder, 'Final_annotation', 'domain_table_MADBUB_new.txt'), 'w')
handle_out = open(os.path.join(madbub_folder, 'domain_table_MADBUB_new.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

# Write headers
writer.writerow(['Sequence'] + domain_list)

# Look at each of the sequences
for sequence in madbub_position_dict.keys():
    new_line = [sequence]
    
    # Loop through the list of domains to get the table on which sequences have which domains
    for domain in domain_list:
        # Get a 1 if that domain was found in that sequence or a 0 otherwise
        new_line.append(madbub_position_dict[sequence].get(domain, 0))
        
    # Write the new line
    writer.writerow(new_line)
    
handle_out.close()


# Write a table that indicates presence (1) or absence (0) of a domain

# handle_in = open(os.path.join(madbub_folder,'Final_annotation', 'domain_table_MADBUB_new.txt'), 'r')
handle_in = open(os.path.join(madbub_folder, 'domain_table_MADBUB_new.txt'), 'r')
reader = csv.reader(handle_in, delimiter = '\t')

# handle_out = open(os.path.join(madbub_folder,'Final_annotation', 'domain_table_MADBUB_discrete_new.txt'), 'w')
handle_out = open(os.path.join(madbub_folder, 'domain_table_MADBUB_discrete_new.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

header = next(reader)
writer.writerow(header)

# Loop through the entries
for entry in reader:
    # If an entry contains positions for a domain change them to a 1
    for i in range(1, len(entry)):
        if entry[i] != '0':
            entry[i] = 1
    
    writer.writerow(entry)

handle_out.close()




