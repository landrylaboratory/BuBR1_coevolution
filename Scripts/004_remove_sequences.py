#### This script will receive a reference fasta file and a list of IDs. It will
#### then remove the sequences for such IDs from the fasta file.

# Import libraries
from Bio.PDB import *
from Bio import SeqIO
from Bio.Seq import *
from Bio.SeqRecord import *
from Bio.Align.Applications import MuscleCommandline
import csv
import re
import os
import glob
from collections import OrderedDict
from shutil import copyfile
import argparse

parser = argparse.ArgumentParser(description = 'This script will remove a list of sequences from a reference fasta file.')

parser.add_argument('-l', type = str, help = 'The path to the list of IDs', dest = 'ID_list')
parser.add_argument('-f', type = str, help = 'The reference fasta file', dest = 'ref_fasta')
parser.add_argument('-o', type = str, help = 'The path to the output file', dest = 'out_file')

args = parser.parse_args()
paralog_list_file = args.ID_list
ref_fasta = args.ref_fasta
out_file = args.out_file

# Parse the file
proteome = SeqIO.parse(ref_fasta, "fasta")
proteome_dict = {}
for record in proteome:
    proteome_dict[record.id] = record

# Start reading the list of IDs
handle = open(paralog_list_file, 'r')
reader = csv.reader(handle, delimiter = "\t")

# Lists to save everything
records = []
seqs_to_remove = []

# Get the list of sequences that will be removed
for entry in reader:
    seqs_to_remove.append(entry[0])


# Loop through the list of IDs
for protein_id, record in proteome_dict.items():

    if protein_id in seqs_to_remove:
        continue
    
    records.append(record)

# Write records to the output file
SeqIO.write(records, out_file, 'fasta')




