#!/usr/bin/env python
# coding: utf-8

# # Domain annotation Jackhmmer
# 
# This script will use the jackhmmer searches and save the best hits of each of the domains as the annotation.

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
bubr1_folder = '../Data/Domain_annotation/BuBR1/'
bub1_folder = '../Data/Domain_annotation/BUB1/'


# ## 1.- Look at the BUBR1 sequences 

# Use the jackhmmer files to produce Pfam_scan formatted file as detailed in https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Protein+Functional+Analysis+Output+Examples.


jackhmmer_domtblout = os.path.join(bubr1_folder, 'all_Ensembl_Tromer_BuBR1_domains_jackhmmer.domtblout')

domains_BuBR1 = parse_jackhmmer(jackhmmer_domtblout)

# Initialize a dictionary that will save the information about the domains
domains_BuBR1_dict = OrderedDict()

# List the sequences that were duplicated in the Ensembl Orthologs set and Tromer's set
duplicates = ['FALB_MAD', 'LOCU_MAD', 'ACAR_MAD', 'TRUB_MAD', 'MMUS_MAD', 'HSAP_MAD', 'DRER_MAD_B']

# A list that will keep track of all the domains used as a reference
domain_list = []

proteins = []

# Open output file
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
    if '_MADBUB' in sequence_name or '_BUB' in sequence_name or sequence_name in duplicates:
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
        domain = 'kinase-like'
    elif entry[7] in ['ABBA1', 'ABBA2', 'ABBA_other']:
        domain = 'ABBA'
    elif entry[7] in ['KEN1', 'KEN2', 'Other_KEN']:
        domain = 'KEN'
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
                elif overlap > 0:
                    # Otherwise, discard the proposed domain annotation
                    print('Discarded domain: ', protein, domain)
                    print('New observation:', ali_start, ali_end)
                    print('Current annotation:', current_start, current_end)
                    print('------')
                    new_domain = False
 
            # Instances of domains that did not overlap with previous ones will be added as new instances
            if new_domain:
                domain_num = max(domain_consensus_dict[protein][domain].keys()) + 1
                domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
                
in_file.close()                 


# Save the best hits as a Pfam_scan formatted table.

outfile = open(os.path.join(bubr1_folder, 'BuBR1_best_hit_domain_annotation_table.txt'), 'w')
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
            query_length = int(positions[6])
            evalue = positions[7]
            prob = positions[8]
            
            # Change domain num to other if greater than 2
            if domain_num > 2:
                domain_num = 'other'
            
            new_row = [protein, domain_start, domain_end, domain_start, domain_end, 'Best_hit', query_prot, domain + '_' + str(domain_num), query_start, query_end, query_length, score, evalue, prob, 'No_clan']
            writer.writerow(new_row)

outfile.close()


# Save the best hits with their sequences as a FASTA formatted file.

records = []

sequences = SeqIO.parse('../Data/Sequences/Ensembl_Tromer_all_BuBR1_no_duplicates.fa', 'fasta')

# Save a dictionary to simplify access to the sequences
seq_dict = OrderedDict()
for seq in sequences:
    seq_dict[seq.id] = seq.seq

# Update the domain list
domain_list = []

for protein, domain_dict in domain_consensus_dict.items():
    for domain, instance_list in domain_dict.items():
        for domain_num, positions in instance_list.items():
            domain_start = int(positions[0])
            domain_end = int(positions[1])
            
            # Change domain num to other if greater than 2
            if domain_num > 2:
                domain_num = 'other'
            
            # Subtract 1 to the start to change from 1-based notation to 0-based notation
            sequence = seq_dict[protein][domain_start-1:domain_end]
       
            domain_name = domain + '_' + str(domain_num)
            
            # Produce the output in Tromer's notation
            new_id = protein + '||' + domain_name + '/' + str(domain_start) + '-' + str(domain_end)
            
            if not domain_name in domain_list:
                domain_list.append(domain_name)
            
            record = SeqRecord(sequence, id = new_id, description = '')
            
            records.append(record)
            
SeqIO.write(records, os.path.join(bubr1_folder, 'BuBR1_best_hit_domain_annotation.fasta'), 'fasta')


# The following command saves a discretized table of presence (1) or absence (0) of each domain in each sequence.
# Open output file
handle_out = open(os.path.join(bubr1_folder, 'domain_table_BuBR1.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

writer.writerow(['Sequence'] + domain_list)

# Look at each of the sequences
for protein, domain_dict in domain_consensus_dict.items():
    new_line = [protein]
    
    # Loop through the list of domains to get the table on which sequences have which domains
    for domain_name in domain_list:
        # Retrieve the domain from the list
        domain = domain_name.split('_')[0]
        domain_num = domain_name.split('_')[1]
        
        if domain_num == 'other':
            domain_num = 3
        else:
            domain_num = int(domain_num)
        
        domain_check = domain_consensus_dict[protein].get(domain, 0)
        
        if domain_check != 0:
            domain_check = domain_check.get(domain_num, 0)
            if domain_check != 0:
                domain_check = 1
        
        # Get a 1 if that domain was found in that sequence or a 0 otherwise
        new_line.append(domain_check)
        
    # Write the new line
    writer.writerow(new_line)
    
handle_out.close()

# ## 2.- Work with the BUB1 sequences.

jackhmmer_domtblout = os.path.join(bub1_folder, 'all_Ensembl_Tromer_BUB1_domains_jackhmmer.domtblout')

domains_BUB1 = parse_jackhmmer(jackhmmer_domtblout)

# Initialize a dictionary that will save the information about the domains
domains_BUB1_dict = OrderedDict()

# List the sequences that were duplicated in the Ensembl Orthologs set and Tromer's set
duplicates = ['LOCU_BUB', 'FALB_BUB', 'XTRO_BUB', 'ACAR_BUB', 'DRER_BUB', 'HSAP_BUB', 'MMUS_BUB']

# A list that will keep track of all the domains used as a reference
domain_list = []

# Open output file
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
    if '_MADBUB' in sequence_name or '_MAD' in sequence_name or sequence_name in duplicates:
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
                elif overlap > 0:
                    # Otherwise, discard the proposed domain annotation
                    print('Discarded domain: ', protein, domain)
                    print('New observation:', ali_start, ali_end)
                    print('Current annotation:', current_start, current_end)
                    print('------')
                    new_domain = False
 
            # Instances of domains that did not overlap with previous ones will be added as new instances
            if new_domain:
                domain_num = max(domain_consensus_dict[protein][domain].keys()) + 1
                domain_consensus_dict[protein][domain][domain_num] = [ali_start, ali_end, score, query_protein, query_start, query_end, query_length, evalue, prob]
                
in_file.close()    

# Save the best hits as a Pfam_scan formatted table.
outfile = open(os.path.join(bub1_folder, 'BUB1_best_hit_domain_annotation_table.txt'), 'w')
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
            query_length = int(positions[6])
            evalue = positions[7]
            prob = positions[8]
            
            # Change domain num to other if greater than 2
            if domain_num > 2:
                domain_num = 'other'
            
            new_row = [protein, domain_start, domain_end, domain_start, domain_end, 'Best_hit', query_prot, domain + '_' + str(domain_num), query_start, query_end, query_length, score, evalue, prob, 'No_clan']
            writer.writerow(new_row)

outfile.close()


# Save the best hits with their sequences as a FASTA formatted file.
records = []

sequences = SeqIO.parse('../Data/Sequences/Ensembl_Tromer_all_BUB1_no_duplicates.fa', 'fasta')

# Save a dictionary to simplify access to the sequences
seq_dict = OrderedDict()
for seq in sequences:
    seq_dict[seq.id] = seq.seq

# Update the domain list
domain_list = []
    
for protein, domain_dict in domain_consensus_dict.items():
    for domain, instance_list in domain_dict.items():
        for domain_num, positions in instance_list.items():
            domain_start = int(positions[0])
            domain_end = int(positions[1])
            
            # Change domain num to other if greater than 2
            if domain_num > 2:
                domain_num = 'other'
            
            # Subtract 1 to the start to change from 1-based notation to 0-based notation
            sequence = seq_dict[protein][domain_start-1:domain_end]
            
            domain_name = domain + '_' + str(domain_num)
            
            # Produce the output in Tromer's notation
            new_id = protein + '||' + domain_name + '/' + str(domain_start) + '-' + str(domain_end)
            
            if not domain_name in domain_list:
                domain_list.append(domain_name)
            
            record = SeqRecord(sequence, id = new_id, description = '')
            
            records.append(record)
            
SeqIO.write(records, os.path.join(bub1_folder, 'BUB1_best_hit_domain_annotation.fasta'), 'fasta')

# The next cell discretizes domain presence or absence.
# Open output file
handle_out = open(os.path.join(bub1_folder, 'domain_table_BUB1.txt'), 'w')
writer = csv.writer(handle_out, delimiter = '\t')

writer.writerow(['Sequence'] + domain_list)

# Look at each of the sequences
for protein, domain_dict in domain_consensus_dict.items():
    new_line = [protein]
    
    # Loop through the list of domains to get the table on which sequences have which domains
    for domain_name in domain_list:
        # Retrieve the domain from the list
        domain = domain_name.split('_')[0]
        domain_num = domain_name.split('_')[1]
        
        if domain_num == 'other':
            domain_num = 3
        else:
            domain_num = int(domain_num)
        
        domain_check = domain_consensus_dict[protein].get(domain, 0)
        
        if domain_check != 0:
            domain_check = domain_check.get(domain_num, 0)
            if domain_check != 0:
                domain_check = 1
        
        # Get a 1 if that domain was found in that sequence or a 0 otherwise
        new_line.append(domain_check)
        
    # Write the new line
    writer.writerow(new_line)
    
handle_out.close()

