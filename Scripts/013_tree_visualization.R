#### 013_tree_visualization ####

# This script will load the output from MrBayes to visualize the phylogenetic tree

# Load libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggsignif)
library(heatmaply)
library(RColorBrewer)
library(scales)
library(cowplot)
library(ggtree)
library(treeio)

# Set working directory
setwd('<path_to_main_folder>')

# Load tree
tree_madbub <- read.mrbayes(file = 'Results/Species_tree_files/Tree_sequences_aln_mafft_newIDs.nex.con.tre')

# Let's format the tree labels
taxa_labels <- as.data.frame(get.tree(tree_madbub)$tip.label)
colnames(taxa_labels) <- c('Original_label')

taxa_labels2 <- taxa_labels %>%
  mutate(New_label = gsub(pattern = '_MADBUB', replacement = ' (MADBUB)', x = Original_label)) %>%
  mutate(New_label = gsub(pattern = '_MAD', replacement = ' (MAD)', x = New_label)) %>%
  mutate(New_label = gsub(pattern = '_BUB', replacement = ' (BUB)', x = New_label)) %>%
  mutate(New_label = gsub(pattern = '_', replacement = ' ', x = New_label))

# Include the new names  
tree_madbub <-  rename_taxa(tree_madbub, taxa_labels2, Original_label, New_label)

# Plot the tree
p <- tree_madbub %>%
  ggtree() + 
  geom_tiplab(align=T, size = 2.5, linesize=.3, fontface = 'italic') +
  geom_treescale(x=0, y=80, offset=2)
# This is the tree without the highlighted blocks
ggsave('Figures/Tree_species.pdf',
       device = cairo_pdf, width = 14, height = 24, plot = p , dpi = 500)  

# Show the node numbers to allow formatting
p_node_numbers <- p +
  geom_text2(aes(label=node), hjust=-.3, size = 2)
ggsave('Figures/Tree_species_node_numbers.pdf',
       device = cairo_pdf, width = 14, height = 24, plot = p_numbers, dpi = 500) 

# Use the previous figure to change the node numbers in the code below accordingly
p <- p +
  ## Labels for mammals ##
  geom_cladelabel(316, 'Mammals', offset=0.45, offset.text=.05, 
                  barsize=3, color='darkgreen', angle=90, hjust='center', align = TRUE,
                  fontsize = 3) +
  geom_hilight(node=316, fill='darkgreen', alpha=.4) +
  geom_cladelabel(218, 'Mammals', offset=0.45, offset.text=.05, 
                  barsize=3, color='darkgreen', angle=90, hjust='center', align = TRUE,
                  fontsize = 3) +
  geom_hilight(node=218, fill='darkgreen', alpha=.4) +
  ## Labels for fish ##
  geom_cladelabel(407, 'Fish', offset=0.45, offset.text=.05, 
                  barsize=3, color='purple', angle=90, hjust='center', align = TRUE,
                  fontsize = 3) +
  geom_hilight(node=407, fill='purple', alpha=.4) +
  geom_cladelabel(309, 'Fish', offset=0.45, offset.text=.05, 
                  barsize=3, color='purple', angle=90, hjust='center', align = TRUE,
                  fontsize = 3) +
  geom_hilight(node=309, fill='purple', alpha=.4) +
  ## Labels for birds ##
  geom_cladelabel(281, 'Birds', offset=0.45, offset.text=.05, 
                  barsize=3, color='orange', angle=90, hjust='center', align = TRUE,
                  fontsize = 3) +
  geom_hilight(node=281, fill='orange', alpha=.4) +
  geom_cladelabel(378, 'Birds', offset=0.45, offset.text=.05, 
                  barsize=3, color='orange', angle=90, hjust='center', align = TRUE, 
                  fontsize = 3) +
  geom_hilight(node=378, fill='orange', alpha=.4) +
  ## Labels for Reptiles ##
  geom_hilight(node=405, fill='pink', alpha=.4) +
  geom_hilight(node=404, fill='pink', alpha=.4) +
  geom_hilight(node=178, fill='pink', alpha=.4) +
  geom_hilight(node=306, fill='pink', alpha=.4) +
  geom_hilight(node=307, fill='pink', alpha=.4) +
  geom_hilight(node=75, fill='pink', alpha=.4) +
  ## Labels for MAD and BUB ##
  geom_cladelabel(313, 'MAD-like (BUBR1)', offset=0.6, offset.text=.05, 
                  barsize=3, color='steelblue', angle=90, hjust='center', align =TRUE) +
  geom_cladelabel(215, 'BUB-like (BUB1)', offset=0.6, offset.text=.05, 
                  barsize=3, color='red', angle=90, hjust='center', align =TRUE)
p

# Flipping the parent node of all BUB1 sequences and the parent node of all BUBR1 sequences
# allows having the MADBUB sequences at the bottom of the tree
p_flipped <- flip(p, 215, 312) +
  ## Labels for reptiles ##
  geom_strip('Crocodylus porosus (BUB)', 'Salvator merianae (BUB)', color='pink', 
             label="Reptiles",  align = TRUE, offset=0.45, offset.text=.075, 
             barsize=3, angle=90, hjust='center', fontsize = 3) +
  geom_strip('Chrysemys picta bellii (MAD)', 'Salvator merianae (MAD)', color='pink', 
             label="Reptiles",  align = TRUE, offset=0.45, offset.text=.075, 
             barsize=3, angle=90, hjust='center', fontsize = 3)
p_flipped

ggsave('Figures/SupplFigWholeTree.pdf',
       device = cairo_pdf, width = 14, height = 24, plot = p , dpi = 500)  

## Collapse by order ##

# Plot the tree like above with the node numbers and edit accordingly
p_node_numbers <- tree_madbub %>%
  ggtree() + 
  geom_tiplab(align=T, size = 2, linesize=.3) +
  geom_treescale(x=0, y=30, offset=2) +
  geom_text2(aes(label=node), hjust=-.3, size = 2)
p_node_numbers
ggsave('Figures/Tree_species_node_numbers.pdf',
       device = cairo_pdf, width = 14, height = 24, plot = p , dpi = 500) 

# Tree without node numbers and tip labels to be used for the final figure
p <- tree_madbub %>%
  ggtree() + 
  geom_treescale(x=0, y=30, offset=2)
p
ggsave('Figures/Tree_species.pdf',
       device = cairo_pdf, width = 14, height = 24, plot = p , dpi = 500) 

# Add the labels
p <- p +
  ## Labels for mammals ##
  geom_cladelabel(316, 'Mammals', offset=0.45, offset.text=.05, 
                  barsize=3, color='darkgreen', angle=90, hjust='center', align = TRUE) +
  geom_hilight(node=316, fill='darkgreen', alpha=.4) +
  geom_cladelabel(218, 'Mammals', offset=0.45, offset.text=.05, 
                  barsize=3, color='darkgreen', angle=90, hjust='center', align = TRUE) +
  geom_hilight(node=218, fill='darkgreen', alpha=.4) +
  ## Labels for fish ##
  geom_cladelabel(407, 'Fish', offset=0.45, offset.text=.05, 
                  barsize=3, color='purple', angle=90, hjust='center', align = TRUE) +
  geom_hilight(node=407, fill='purple', alpha=.4) +
  geom_cladelabel(309, 'Fish', offset=0.45, offset.text=.05, 
                  barsize=3, color='purple', angle=90, hjust='center', align = TRUE) +
  geom_hilight(node=309, fill='purple', alpha=.4) +
  ## Labels for birds ##
  geom_cladelabel(281, 'Birds', offset=0.45, offset.text=.05, 
                  barsize=3, color='orange', angle=90, hjust='center', align = TRUE) +
  geom_hilight(node=281, fill='orange', alpha=.4) +
  geom_cladelabel(378, 'Birds', offset=0.45, offset.text=.05, 
                  barsize=3, color='orange', angle=90, hjust='center', align = TRUE) +
  geom_hilight(node=378, fill='orange', alpha=.4) +
  ## Labels for Reptiles ##
  geom_hilight(node=405, fill='pink', alpha=.4) +
  geom_hilight(node=404, fill='pink', alpha=.4) +
  geom_hilight(node=178, fill='pink', alpha=.4) +
  geom_hilight(node=306, fill='pink', alpha=.4) +
  geom_hilight(node=307, fill='pink', alpha=.4) +
  geom_hilight(node=75, fill='pink', alpha=.4) +
  # Labels for MAD and BUB
  geom_cladelabel(313, 'MAD', offset=0.6, offset.text=.05, 
                  barsize=3, color='steelblue', angle=90, hjust='center', align =TRUE) +
  geom_cladelabel(215, 'BUB', offset=0.6, offset.text=.05, 
                  barsize=3, color='red', angle=90, hjust='center', align =TRUE)
p

# Flip to have the MADBUBs at the bottom
p_flipped <- flip(p, 215, 312)
p_flipped

p3 <- p_flipped %>%
  # Primates #
  collapse(node=247) %>%
  collapse(node=345) %>%
  # Carnivora #
  collapse(node=325) %>%
  collapse(node=225) %>%
  # Chiroptera #
  collapse(node=343) %>%
  collapse(node=234) %>%
  # Galliformes #
  collapse(node=392) %>%
  collapse(node=289) %>%
  # Marsupialia #
  collapse(node=274) %>%
  collapse(node=371) %>%
  # Afrotheria #
  collapse(node=369) %>%
  collapse(node=272) %>%
  # Anseriformes #
  collapse(node=397) %>%
  collapse(node=294) %>%
  # Paleognathae #
  collapse(node=399) %>%
  collapse(node=296) %>%
  # Passeriformes #
  collapse(node=381) %>%
  collapse(node=282) %>%
  # Accipitriformes #
  collapse(node=390) %>%
  collapse(node=303) %>%
  # Charadriiformes #
  collapse(node=387) %>%
  collapse(node=302) %>%
  # Psittaciformes #
  collapse(node=304) %>%
  collapse(node=388) %>%
  # Osteoglossiformes #
  collapse(node=311) %>%
  collapse(node=409) %>%
  # Xenarthra #
  collapse(node=368) %>%
  collapse(node=271) %>%
  # Cetartiodactyla #
  collapse(node=332) %>%
  collapse(node=236) +
  # Primates #
  geom_point2(aes(subset=(node==247)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 247, label = 'Primates', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==345)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 345, label = 'Primates', fontsize = 3, align = TRUE) +
  # Carnivora #
  geom_point2(aes(subset=(node==325)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 325, label = 'Carnivora', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==225)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 225, label = 'Carnivora', fontsize = 3, align = TRUE) +
  # Marsupialia #
  geom_point2(aes(subset=(node==274)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 274, label = 'Marsupialia', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==369)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 369, label = 'Marsupialia', fontsize = 3, align = TRUE) +
  # Afrotheria #
  geom_point2(aes(subset=(node==272)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 272, label = 'Afrotheria', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==371)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 371, label = 'Afrotheria', fontsize = 3, align = TRUE) +
  # Chiroptera #
  geom_point2(aes(subset=(node==343)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 343, label = 'Chiroptera', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==234)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 234, label = 'Chiroptera', fontsize = 3, align = TRUE) +
  # Galliformes #
  geom_point2(aes(subset=(node==392)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 392, label = 'Galliformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==289)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 289, label = 'Galliformes', fontsize = 3, align = TRUE) +
  # Anseriformes #
  geom_point2(aes(subset=(node==397)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 397, label = 'Anseriformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==294)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 294, label = 'Anseriformes', fontsize = 3, align = TRUE) +
  # Paleognathae #
  geom_point2(aes(subset=(node==399)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 399, label = 'Paleognathae', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==296)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 296, label = 'Paleognathae', fontsize = 3, align = TRUE) +
  # Passeriformes #
  geom_point2(aes(subset=(node==381)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 381, label = 'Passeriformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==282)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 282, label = 'Passeriformes', fontsize = 3, align = TRUE) +
  # Accipitriformes #
  geom_point2(aes(subset=(node==390)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 390, label = 'Accipitriformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==303)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 303, label = 'Accipitriformes', fontsize = 3, align = TRUE) +
  # Charadriiformes #
  geom_point2(aes(subset=(node==387)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 387, label = 'Charadriiformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==302)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 302, label = 'Charadriiformes', fontsize = 3, align = TRUE) +
  # Pscittaciformes #
  geom_point2(aes(subset=(node==304)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 304, label = 'Pscittaciformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==388)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 388, label = 'Pscittaciformes', fontsize = 3, align = TRUE) +
  # Lepisosteiformes #
  geom_point2(aes(subset=(node==108)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 108, label = 'Lepisosteiformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==5)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 5, label = 'Lepisosteiformes', fontsize = 3, align = TRUE) +
  # Polypteriformes #
  geom_point2(aes(subset=(node==85)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 85, label = 'Polypteriformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==188)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 188, label = 'Polypteriformes', fontsize = 3, align = TRUE) +
  # Xenarthra #
  geom_point2(aes(subset=(node==271)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 271, label = 'Xenarthra', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==368)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 368, label = 'Xenarthra', fontsize = 3, align = TRUE) +
  # Monotremata #
  geom_point2(aes(subset=(node==111)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 111, label = 'Monotremata', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==8)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 8, label = 'Monotremata', fontsize = 3, align = TRUE) +
  # Erinacidaceae #
  geom_point2(aes(subset=(node==123)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 123, label = 'Erinacidaceae', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==20)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 20, label = 'Erinacidaceae', fontsize = 3, align = TRUE) +
  # Soricidaceae #
  geom_point2(aes(subset=(node==18)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 18, label = 'Soricidaceae', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==121)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 121, label = 'Soricidaceae', fontsize = 3, align = TRUE) +
  # Perissodactyla #
  geom_point2(aes(subset=(node==181)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 181, label = 'Perissodactyla', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==78)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 78, label = 'Perissodactyla', fontsize = 3, align = TRUE) +
  # Anura #
  geom_point2(aes(subset=(node==185)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 185, label = 'Anura', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==82)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 82, label = 'Anura', fontsize = 3, align = TRUE) +
  # Petromyzontiformes #
  geom_point2(aes(subset=(node==207)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 207, label = 'Petromyzontiformes (MADBUB)', fontsize = 3, align = TRUE) +
  # Araneae #
  geom_point2(aes(subset=(node==211)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 211, label = 'Araneae (MADBUB)', fontsize = 3, align = TRUE) +
  # Enterogona #
  geom_point2(aes(subset=(node==209)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 209, label = 'Enterogona (MADBUB)', fontsize = 3, align = TRUE) +
  # Amphioxiformes #
  geom_point2(aes(subset=(node==208)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 208, label = 'Amphioxiformes (MADBUB)', fontsize = 3, align = TRUE) +
  # Echinoida #
  geom_point2(aes(subset=(node==210)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 210, label = 'Echinoida (MADBUB)', fontsize = 3, align = TRUE) +
  # Cetartiodactyla #
  geom_point2(aes(subset=(node==332)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 332, label = 'Cetartiodactyla', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==236)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 236, label = 'Cetartiodactyla', fontsize = 3, align = TRUE)
p3

p4 <- p3 +
  geom_strip('Crocodylus_porosus_BUB', 'Salvator_merianae_BUB', color='pink', 
             label="Reptiles",  align = TRUE, offset=0.45, offset.text=.05, 
             barsize=3, angle=90, hjust='center') +
  geom_strip('Chrysemys_picta_bellii_MAD', 'Salvator_merianae_MAD', color='pink', 
             label="Reptiles",  align = TRUE, offset=0.45, offset.text=.05, 
             barsize=3, angle=90, hjust='center')
p4

p_order <- p4 %>%
  # Squamata #
  collapse(node=405) %>%
  collapse(node=307) %>%
  # Testudines #
  collapse(node=404) %>%
  collapse(node=306) +
  # Testudines #
  geom_point2(aes(subset=(node==404)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 404, label = 'Testudines', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==306)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 306, label = 'Testudines', fontsize = 3, align = TRUE) +
  # Squamata #
  geom_point2(aes(subset=(node==405)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 405, label = 'Squamata', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==307)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 307, label = 'Squamata', fontsize = 3, align = TRUE) +
  # Osteoglossiformes #
  geom_point2(aes(subset=(node==409)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 409, label = 'Osteoglossiformes', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==311)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 311, label = 'Osteoglossiformes', fontsize = 3, align = TRUE) +
  # Crocodilia #
  geom_point2(aes(subset=(node==75)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 75, label = 'Crocodilia', fontsize = 3, align = TRUE) +
  geom_point2(aes(subset=(node==178)), shape=19, size=2, fill='black') +
  geom_cladelabel(node = 178, label = 'Crocodilia', fontsize = 3, align = TRUE)
p_order

ggsave('Figures/Figure1D.pdf',
       device = cairo_pdf, width = 10, height = 20, plot = p_order, dpi = 500)
