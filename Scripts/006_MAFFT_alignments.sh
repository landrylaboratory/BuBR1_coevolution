#####################################################
####		006_MAFFT_alignments		 ####
#### This script performs aligns the BUB1 and	 ####
#### BUBR1 sequences separately. It then aligns	 ####
#### the domains that the hBUB1 and hBUBR1	 ####
#### sequences have in common.			 ####
#####################################################

mkdir ../Data/MAFFT_Alignments
mkdir ../Data/Sequences/Human_domains
mkdir ../Data/Sequences/Human_domains/Alignments

#### Aligning the two complete sets of sequences ####

mafft --ep 0 --genafpair --maxiterate 1000 ../Data/Sequences/Human_BUB1B_orthologues_2020-04-27.fa > ../Data/MAFFT_Alignments/Human_BUB1B_orthologues_2020-04-27_mafft.fasta

mafft --ep 0 --genafpair --maxiterate 1000 ../Data/Sequences/Human_BUB1_orthologues_2020-04-27.fa > ../Data/MAFFT_Alignments/Human_BUB1_orthologues_2020-04-27_mafft.fasta

#### Aligning the common human BUB1 and BUBR1 domains ####

# Save the sequence IDs we need

grep 'HSAP' ../Data/Sequences/Tromer_domain_annotation_plk1.fasta | cut -c 2- > ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt

# Start pairing the corresponding domains
# TPR
grep 'TPR' ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt > ../Data/Sequences/Human_domains/Human_TPR_Tromer_IDs.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Human_domains/Human_TPR_Tromer_IDs.txt -f ../Data/Sequences/Tromer_domain_annotation_plk1.fasta -o ../Data/Sequences/Human_domains/Human_TPR_Tromer_sequences.fasta

mafft ../Data/Sequences/Human_domains/Human_TPR_Tromer_sequences.fasta > ../Data/Sequences/Human_domains/Alignments/Human_TPR_Tromer_sequences_aln.fasta

# GLEBS
grep 'GLEBS' ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt > ../Data/Sequences/Human_domains/Human_GLEBS_Tromer_IDs.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Human_domains/Human_GLEBS_Tromer_IDs.txt -f ../Data/Sequences/Tromer_domain_annotation_plk1.fasta -o ../Data/Sequences/Human_domains/Human_GLEBS_Tromer_sequences.fasta

mafft ../Data/Sequences/Human_domains/Human_GLEBS_Tromer_sequences.fasta > ../Data/Sequences/Human_domains/Alignments/Human_GLEBS_Tromer_sequences_aln.fasta

# ABBA3 (ABBA_other)
# Both ABBA3 and ABBA4 are labeled as "other", so the way to distinguish them is by their position. 
# ABBA3 starts at position 524 in BUB1 and in position 539 in BUBR1, while ABBA4 starts in position 616 in BUB1
# The /5 allows working only with ABBA3 
grep 'ABBA_other/5' ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt > ../Data/Sequences/Human_domains/Human_ABBA3_Tromer_IDs.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Human_domains/Human_ABBA3_Tromer_IDs.txt -f ../Data/Sequences/Tromer_domain_annotation_plk1.fasta -o ../Data/Sequences/Human_domains/Human_ABBA3_Tromer_sequences.fasta

mafft ../Data/Sequences/Human_domains/Human_ABBA3_Tromer_sequences.fasta > ../Data/Sequences/Human_domains/Alignments/Human_ABBA3_Tromer_sequences_aln.fasta

# PLKBD
grep 'PLKBD' ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt > ../Data/Sequences/Human_domains/Human_PLKBD_Tromer_IDs.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Human_domains/Human_PLKBD_Tromer_IDs.txt -f ../Data/Sequences/Tromer_domain_annotation_plk1.fasta -o ../Data/Sequences/Human_domains/Human_PLKBD_Tromer_sequences.fasta

mafft ../Data/Sequences/Human_domains/Human_PLKBD_Tromer_sequences.fasta > ../Data/Sequences/Human_domains/Alignments/Human_PLKBD_Tromer_sequences_aln.fasta

# KARD
grep 'KARD' ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt > ../Data/Sequences/Human_domains/Human_KARD_Tromer_IDs.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Human_domains/Human_KARD_Tromer_IDs.txt -f ../Data/Sequences/Tromer_domain_annotation_plk1.fasta -o ../Data/Sequences/Human_domains/Human_KARD_Tromer_sequences.fasta

mafft ../Data/Sequences/Human_domains/Human_KARD_Tromer_sequences.fasta > ../Data/Sequences/Human_domains/Alignments/Human_KARD_Tromer_sequences_aln.fasta

# CDII
grep 'CDII' ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt > ../Data/Sequences/Human_domains/Human_CDII_Tromer_IDs.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Human_domains/Human_CDII_Tromer_IDs.txt -f ../Data/Sequences/Tromer_domain_annotation_plk1.fasta -o ../Data/Sequences/Human_domains/Human_CDII_Tromer_sequences.fasta

mafft ../Data/Sequences/Human_domains/Human_CDII_Tromer_sequences.fasta > ../Data/Sequences/Human_domains/Alignments/Human_CDII_Tromer_sequences_aln.fasta

# Kinase
grep 'kinase' ../Data/Sequences/Human_domains/Human_domains_Tromer_IDs.txt > ../Data/Sequences/Human_domains/Human_kinase_Tromer_IDs.txt

python3 002_extract_sequences.py -l ../Data/Sequences/Human_domains/Human_kinase_Tromer_IDs.txt -f ../Data/Sequences/Tromer_domain_annotation_plk1.fasta -o ../Data/Sequences/Human_domains/Human_kinase_Tromer_sequences.fasta

mafft ../Data/Sequences/Human_domains/Human_kinase_Tromer_sequences.fasta > ../Data/Sequences/Human_domains/Alignments/Human_kinase_Tromer_sequences_aln.fasta


