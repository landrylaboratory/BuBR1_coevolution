##################################################
####        009_entropy_analysis              ####
#### This script will look at the entropy     ####
#### values to show them in fugures.          ####
##################################################

# Load libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(magrittr)
library(ggpubr)
library(ggsignif)
library(heatmaply)
library(RColorBrewer)
library(scales)
theme_set(theme_cowplot())

# Set working directory
setwd('<path_to_main_folder>/Data/Evol_results/MAFFT')

#### Load data ####

## Original filter
entropy_bubr1 <- read_delim('BUBR1_entropy_2020_05_08.txt', delim = '\t', col_names = c('entropy')) %>%
  mutate(protein = 'BUBR1', position = row_number())

entropy_bub1 <- read_delim('BUB1_entropy_2020_05_08.txt', delim = '\t', col_names = c('entropy')) %>%
  mutate(protein = 'BUB1', position = row_number())

## Working with the extra duplicates

human_bub1_bubr1_aln <- read_delim('../../Sequences/Human_domains/Alignments/Tables/Human_domains_Tromer_table.txt', delim = '\t')

#### Work on the paired boxplots ####

bub1_bubr1_entropy <- rbind(entropy_bub1, entropy_bubr1)

entropy_common_domains <- inner_join(x = bub1_bubr1_entropy, y = human_bub1_bubr1_aln,
                                     by = c('position' = 'Position', 'protein' = 'Protein'))

#### Add a color for the lines to highlight the critical residues. 

entropy_common_domains %<>%
  mutate(Domain = gsub(x = Domain, pattern = 'kinase', replacement = '(Pseudo)kinase'))
  
entropy_common_domains$Domain <- factor(entropy_common_domains$Domain, 
                                        levels = c('TPR', 'GLEBS', 'ABBA3', 'PLKBD', 'KARD', 'CDII', '(Pseudo)kinase'))

check_paired <- as.data.frame(table(entropy_common_domains$Group) == 2)
check_paired <- cbind(as.numeric(rownames(check_paired)),check_paired)
colnames(check_paired) <- c('Group', 'Paired')

entropy_common_domains_new <- inner_join(x = entropy_common_domains, y = check_paired,
                                         by = c('Group' = 'Group'))


#### Mark the critical residues ####

## Kinase (catalytic triad): ##
# BUB1 K821 - BUBR1 K795
# BUB1 D917 - BUBR1 D882
# BUB1 D946 - BUBR1 D911

## KARD (phosphosites) ##
# BUB1 S655 - BUBR1 S670
# BUB1 S661  - BUBR1 S676
# BUB1 A665  - BUBR1 T680


## MVA mutations
# BUB1 I747 - BUBR1 R727
# BUB1 R840 - BUBR1 R814
# BUB1 A877 - BUBR1 L844
# BUB1 L1047 - BUBR1 L1012

## With MVA mutants ##
# Color code:
# 0: Regular residues (gray)
# 1: Kinase catalytic triad (red)
# 2: KARD phosphosites (blue)
# 3: MVA mutations (green)

entropy_common_domains_new %<>%
    mutate(key_res = ifelse(and(protein == 'BUB1', position %in% c(821, 917, 946)), 1, 
                            ifelse(and(protein == 'BUBR1', position %in% c(795, 882, 911)), 1, 
                                 ifelse(and(protein == 'BUB1', position %in% c(655, 661, 665)), 2, 
                                        ifelse(and(protein == 'BUBR1', position %in% c(670, 676, 680)), 2, 
                                           ifelse(and(protein == 'BUB1', position %in% c(747, 840, 877, 1047)), 3, 
                                                  ifelse(and(protein == 'BUBR1', position %in% c(727, 814, 844, 1012)), 3,0)
                                        )
                                 )
                          )
                    )
  )
)

# Draw the figure with all the paired points
p_boxplots_entropy <- entropy_common_domains_new %>% filter(Paired == TRUE) %>%
  ggplot(aes(x = protein, y = entropy, fill = protein)) +
  scale_fill_manual(values = c('#e6e6e6', '#cccccc')) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  facet_grid(cols = vars(Domain)) +
  geom_line(aes(group = Group, colour = factor(key_res), alpha = factor(key_res))) +
  scale_colour_manual(values = c('black', 'red', 'blue', '#009900')) +
  scale_alpha_manual(values = c(0.1, 0.8, 0.8, 0.8)) +
  geom_point() +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", size=1),
        strip.background = element_rect(fill = 'white'),
        legend.position = 'none'
  ) +
  xlab('Protein') + ylab('Entropy') +
  stat_compare_means(comparisons = list(c('BUB1', 'BUBR1')), paired = T, method = 't.test')

p_final <- p_boxplots_entropy +
  geom_point(data = entropy_common_domains_new %>% filter(Paired == FALSE))      
p_final

ggsave(plot = p_final, device = cairo_pdf, width = 10, height = 7, dpi = 500, 
       filename = '../../../Figures/Figure1E.pdf')
