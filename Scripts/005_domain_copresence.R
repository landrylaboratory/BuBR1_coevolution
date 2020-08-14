##################################################
####    005_domain_copresence                 ####
#### This script looks at table of presence   ####
#### of domains in each sequence and produces ####
#### figures on domain composition.           ####
##################################################

#### Load libraries ####

library(tidyverse)
library(magrittr)
library(ggplot2)
library(Cairo)
library(scales)
library(stringr)

# Set working directory to the scripts folder
setwd('<path to scripts folder>')

#### Define functions ####

# A function for the other definition of correlation (Wu,2003)
get_correlation_wu <- function(col1, col2, data){
  column1 <- data[, col1]
  column2 <- data[, col2]
  
  present_dom1 <- sum(column1)
  present_dom2 <- sum(column2)
  
  N <- nrow(column1)
  z <- sum(as.numeric((column1 + column2) == 2))
  x <- min(present_dom1, present_dom2)
  y <- max(present_dom1, present_dom2)
  
  corr <- (N*z - x*y) / sqrt((N*x - x^2)*(N*y - y^2))
  return(corr)
}

# A function to calculate the new correlation
calc_copresence_corr <- function(df){
  N <- ncol(df)
  vals <- Vectorize(get_correlation_wu, vectorize.args=list("col1","col2"))
  outer(X = 1:N,Y = 1:N,FUN = vals,data=df)
}

#### MADBUB ####

domains_MADBUB <- read.table('../Data/Domain_annotation/MADBUB/domain_table_MADBUB_discrete_new.txt',
                             h = T, sep = '\t')

colnames(domains_MADBUB)
# Save ABBA5 onwards as ABBA_other, CMI2 onwards as CMI_other, Dbox3 onwards as Dbox_other
domains_MADBUB_new <- domains_MADBUB %>% rowwise() %>%
  mutate(ABBA_other = as.numeric(ABBA_5 || ABBA_6 || ABBA_7 || ABBA_8 || ABBA_9 || ABBA_10 || ABBA_11)) %>%
  select(-ABBA_5, -ABBA_6, -ABBA_7, -ABBA_8, -ABBA_9, -ABBA_10, -ABBA_11) %>%
  mutate(CMI_other = as.numeric(CMI_2 || CMI_3 || CMI_4)) %>%
  select(-CMI_2, -CMI_3, -CMI_4) %>%
  mutate(`D-Box_other` = as.numeric(`D.Box_3` || `D.Box_4` || `D.Box_5` || `D.Box_6` || `D.Box_7` || `D.Box_8` || `D.Box_9` || `D.Box_10`)) %>%
  select(-`D.Box_3`, -`D.Box_4`, -`D.Box_5`, -`D.Box_6`, -`D.Box_7`, -`D.Box_8`, -`D.Box_9`, -`D.Box_10`)

# Change the column names for format
colnames(domains_MADBUB_new)
colnames(domains_MADBUB_new) <- c('Sequence', 'ABBA1', 'ABBA2', 'CDII', 'CMI', 'D-Box1', 'GLEBS1',
                                  'KARD', 'KEN1', 'KEN2', 'Kinase', 'TPR', 'ABBA3', 'ABBA4', 
                                  'D-Box2', 'MadaM', 'KARD2', 'GLEBS2', 'PLKBD', 'ABBA_other', 'CMI_other',
                                  'D-Box_other'
  
)
colnames(domains_MADBUB_new)

final_domains_MADBUB <- domains_MADBUB_new %>% select(-Sequence, -`D-Box_other`,
                                                      -ABBA4, -ABBA_other, -KARD2,
                                                      -GLEBS2, -CMI_other, -MadaM)

# Draw the UpSetR figures to see which combinations of features are most common
final_domains_MADBUB_upset <- domains_MADBUB_new %>% select(-`D-Box_other`,
                                                            -ABBA4, -ABBA_other, -KARD2,
                                                            -GLEBS2, -CMI_other, -MadaM) %>%
  ungroup() %>% as.data.frame()

pdf(file = '../Figures/SupplFig1E.pdf', onefile = F, width = 10, height = 7)
upset(final_domains_MADBUB_upset, order.by = "freq", nsets = 16)
grid.text("MADBUB features",x = 0.65, y=0.95, gp=gpar(fontsize=16))
dev.off()


#### Following Tromer et al. Open Biol 2016, look at correlation in domain presence ####
# Calculate correlation (Wu, 2003)
correlation_table_MADBUB <- final_domains_MADBUB  %>% calc_copresence_corr()

colnames(correlation_table_MADBUB) <- colnames(final_domains_MADBUB)
rownames(correlation_table_MADBUB) <- colnames(final_domains_MADBUB)

# Use hierarchical clustering for the heatmap
row_order <- hclust(dist(correlation_table_MADBUB))$order
col_order <- hclust(dist(t(correlation_table_MADBUB)))$order
correlation_table_MADBUB <- correlation_table_MADBUB[row_order, col_order]

correlation_table_MADBUB_new <- as.data.frame(correlation_table_MADBUB)
correlation_table_MADBUB_new$Domain_1 <- rownames(correlation_table_MADBUB_new)
num_cols_MADBUB <- ncol(correlation_table_MADBUB_new)
correlation_table_MADBUB_new <- correlation_table_MADBUB_new[c(num_cols_MADBUB,1:(num_cols_MADBUB - 1))]

tidy_cor_MADBUB <- correlation_table_MADBUB_new %>%
  gather(key = 'Domain_2', value = 'Correlation', -Domain_1) %>%
  mutate(Correlation = as.numeric(Correlation), 
         Domain_2 = as.character(Domain_2),
         Domain_1 = as.character(Domain_1))

right_order <- tidy_cor_MADBUB$Domain_1[1:(num_cols_MADBUB - 1)]
tidy_cor_MADBUB$Domain_1 <- factor(tidy_cor_MADBUB$Domain_1, levels = right_order)
tidy_cor_MADBUB$Domain_2 <- factor(tidy_cor_MADBUB$Domain_2, levels = right_order)

p <- ggplot(tidy_cor_MADBUB, aes(x = Domain_1, y = Domain_2)) +
  geom_tile(aes(fill = Correlation), colour = 'black', size = 1.1) +
  scale_fill_gradient2(name = 'Pearson correlation', low = muted('blue'), mid = 'white', high = muted('red'),
                       limits = c(-0.5, 1)) +
  xlab('') + ylab('') + 
  ggtitle('MADBUB domain co-presence') +
  theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 18), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(), axis.title.y = element_text(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, size = 18, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.justification = 0.5,
        legend.title = element_text(size = 14),
        legend.text = element_text(angle = -90, hjust = 1, vjust = 0.4))
p

ggsave(filename = '../Figures/Figure1C.pdf',
       plot = p, device = cairo_pdf, width = 10, height = 10, dpi = 500)

#### BuBR1 #### 

# Remove the pseudo_kinase and X_pseudo_kinase assignments
domains_BuBR1 <- read.table('../Data/Domain_annotation/BuBR1/domain_table_BuBR1_discrete_new.txt',
                             h = T, sep = '\t')

domain_counts <- domains_BuBR1 %>% select(-Sequence) %>% colSums()

# Replace column names in the original table
colnames(domains_BuBR1)
colnames(domains_BuBR1) <- c('Sequence', 'ABBA1', 'ABBA2', 'KEN1', 'KEN2',
                             'GLEBS', 'D-Box1', 'TPR', 'CDII', 'KARD',
                             '(Pseudo)kinase', 'PLKBD', 'ABBA3', 'D-Box2', 'ABBA4',
                             'CMI'
                             )
colnames(domains_BuBR1)

final_domains_BuBR1 <- domains_BuBR1 %>% select(-Sequence, -ABBA4)

# Draw the UpSetR figures to see which combinations of features are most common
final_domains_BuBR1_upset <- domains_BuBR1 %>% select(-ABBA4) %>%
  ungroup() %>% as.data.frame()

pdf(file = '../Figures/SupplFig1D.pdf', onefile = F, width = 10, height = 7)
upset(final_domains_BuBR1_upset, order.by = "freq", nsets = 16)
grid.text("BUBR1 features",x = 0.65, y=0.95, gp=gpar(fontsize=16))
dev.off()


#### Following Tromer et al. Open Biol 2016, look at correlation in domain presence ####
# Calculate the correlation (Wu, 2003)
correlation_table_BuBR1 <- final_domains_BuBR1 %>% rowwise() %>% calc_copresence_corr()

colnames(correlation_table_BuBR1) <- colnames(final_domains_BuBR1)
rownames(correlation_table_BuBR1) <- colnames(final_domains_BuBR1)

# Use hierarchical clustering for the heatmap
row_order <- hclust(dist(correlation_table_BuBR1))$order
col_order <- hclust(dist(t(correlation_table_BuBR1)))$order
correlation_table_BuBR1 <- correlation_table_BuBR1[row_order, col_order]

correlation_table_BuBR1_new <- as.data.frame(correlation_table_BuBR1)
correlation_table_BuBR1_new$Domain_1 <- rownames(correlation_table_BuBR1_new)
num_cols_BuBR1 <- ncol(correlation_table_BuBR1_new)
correlation_table_BuBR1_new <- correlation_table_BuBR1_new[c(num_cols_BuBR1,1:(num_cols_BuBR1 - 1))]

tidy_cor_BuBR1 <- correlation_table_BuBR1_new %>%
  gather(key = 'Domain_2', value = 'Correlation', -Domain_1) %>%
  mutate(Correlation = as.numeric(Correlation), 
         Domain_2 = as.character(Domain_2),
         Domain_1 = as.character(Domain_1))

right_order <- tidy_cor_BuBR1$Domain_1[1:(num_cols_BuBR1 - 1)]
tidy_cor_BuBR1$Domain_1 <- factor(tidy_cor_BuBR1$Domain_1, levels = right_order)
tidy_cor_BuBR1$Domain_2 <- factor(tidy_cor_BuBR1$Domain_2, levels = right_order)

p <- ggplot(tidy_cor_BuBR1, aes(x = Domain_1, y = Domain_2)) +
  geom_tile(aes(fill = Correlation), colour = 'black', size = 1.1) +
  scale_fill_gradient2(name = 'Pearson correlation', low = muted('blue'), mid = 'white', high = muted('red'),
                       limits = c(-0.5, 1)) +
  xlab('') + ylab('') + 
  ggtitle('BUBR1 domain co-presence') +
  theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 18), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(), axis.title.y = element_text(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, size = 18, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.justification = 0.5,
        legend.title = element_text(size = 14),
        legend.text = element_text(angle = -90, hjust = 1, vjust = 0.4))
p

ggsave(filename = '../Figures/Figure1B.pdf',
       plot = p, device = cairo_pdf, width = 10, height = 10, dpi = 500)

#### BUB1 ####

domains_BUB1 <- read.table('../Data/Domain_annotation/BUB1/domain_table_BUB1_discrete_new.txt',
                           h = T, sep = '\t')

domain_counts_BUB1 <- domains_BUB1 %>% select(-Sequence) %>% colSums()

colnames(domains_BUB1)
colnames(domains_BUB1) <- c('Sequence', 'ABBA3', 'CDII', 'CDII2', 'TPR', 'Kinase', 
                            'GLEBS', 'CMI', 'KARD', 'PLKBD', 'ABBA4'
  
)
colnames(domains_BUB1)

final_domains_BUB1 <- domains_BUB1 %>% select(-Sequence, -ABBA4, -CDII2)

# Draw the UpSetR figures to see which combinations of features are most common
final_domains_BUB1_upset <- domains_BUB1 %>% select(-ABBA4, -CDII2) %>%
  ungroup() %>% as.data.frame()

pdf(file = '../Figures/SupplFig1C.pdf', onefile = F, width = 10, height = 7)
upset(final_domains_BUB1_upset, order.by = "freq", nsets = 16)
grid.text("BUB1 features",x = 0.65, y=0.95, gp=gpar(fontsize=16))
dev.off()

#### Following Tromer et al. Open Biol 2016, look at correlation in domain presence ####
# Calculate correlations (Wu, 2003)
correlation_table_BUB1 <- final_domains_BUB1 %>% rowwise() %>% calc_copresence_corr()

colnames(correlation_table_BUB1) <- colnames(final_domains_BUB1)
rownames(correlation_table_BUB1) <- colnames(final_domains_BUB1)

# Apply hierarchical clustering
row_order <- hclust(dist(correlation_table_BUB1))$order
col_order <- hclust(dist(t(correlation_table_BUB1)))$order
correlation_table_BUB1 <- correlation_table_BUB1[row_order, col_order]

correlation_table_BUB1_new <- as.data.frame(correlation_table_BUB1)
correlation_table_BUB1_new$Domain_1 <- rownames(correlation_table_BUB1_new)
num_cols_BUB1 <- ncol(correlation_table_BUB1_new)
correlation_table_BUB1_new <- correlation_table_BUB1_new[c(num_cols_BUB1,1:(num_cols_BUB1 - 1))]

tidy_cor_BUB1 <- correlation_table_BUB1_new %>%
  gather(key = 'Domain_2', value = 'Correlation', -Domain_1) %>%
  mutate(Correlation = as.numeric(Correlation), 
         Domain_2 = as.character(Domain_2),
         Domain_1 = as.character(Domain_1))

right_order <- tidy_cor_BUB1$Domain_1[1:num_cols_BUB1-1]
tidy_cor_BUB1$Domain_1 <- factor(tidy_cor_BUB1$Domain_1, levels = right_order)
tidy_cor_BUB1$Domain_2 <- factor(tidy_cor_BUB1$Domain_2, levels = right_order)

p <- ggplot(tidy_cor_BUB1, aes(x = Domain_1, y = Domain_2)) +
  geom_tile(aes(fill = Correlation), colour = 'black', size = 1.1) +
  scale_fill_gradient2(name = 'Pearson correlation', low = muted('blue'), mid = 'white', high = muted('red'),
                       limits = c(-0.5, 1)) +
  xlab('') + ylab('') + 
  ggtitle('BUB1 domain co-presence') +
  theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 18), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(), axis.title.y = element_text(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, size = 18, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.justification = 0.5,
        legend.title = element_text(size = 14),
        legend.text = element_text(angle = -90, hjust = 1, vjust = 0.4))
p

ggsave(filename = '../Figures/Figure1A.pdf',
       plot = p, device = cairo_pdf, width = 10, height = 10, dpi = 500)

#### Check if organisms from the same taxonomic order have the same domain composition ####

# Read the list of sequences from the final datasets
bub1_filtered <- read_delim('../Results/MAFFT/BUB1_IDs.txt', delim = '\t', col_names = c('seq_id'))
bubr1_filtered <- read_delim('../Results/MAFFT/BuBR1_IDs.txt', delim = '\t', col_names = c('seq_id'))

domains_bub1_filtered <- domains_BUB1 %>% filter(Sequence %in% bub1_filtered$seq_id) %>%
  rowwise() %>%
  mutate(Prefix = str_match(Sequence, "(.*)P[0-9]+")[2])

domains_bubr1_filtered <- domains_BuBR1 %>% filter(Sequence %in% bubr1_filtered$seq_id) %>%
  rowwise() %>%
  mutate(Prefix = str_match(Sequence, "(.*)P[0-9]+")[2])

# Load the data on the Ensembl species
ensembl_species <- read_delim('../Data/Ensembl_species_names_no_dups.txt', delim = '\t')
ensembl_species %<>% 
  filter(`Species name` != 'Sus scrofa (Pig USMARC)') %>%
  rowwise() %>%
  mutate(Species_scientific = str_match(`Species name`, "(.*) \\(.*\\)")[2])

# Add the species names to the filtered sequences
bub1_domains_filtered_new <- inner_join(x = domains_bub1_filtered, y = ensembl_species,
                                         by = c("Prefix" = "Prefix"))
bubr1_domains_filtered_new <- inner_join(x = domains_bubr1_filtered, y = ensembl_species,
                                         by = c("Prefix" = "Prefix"))

# Load data on taxonomic orders
taxons <- read_delim('../Data/103_Species_orders.txt', delim = '\t', col_names = c('Species_scientific', 'Taxon'))

bub1_domains_filtered_taxons <- inner_join(x = bub1_domains_filtered_new, y = taxons,
                                            by = c("Species_scientific" = "Species_scientific"))
bubr1_domains_filtered_taxons <- inner_join(x = bubr1_domains_filtered_new, y = taxons,
                                            by = c("Species_scientific" = "Species_scientific"))

# Create the function to get the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Use an average to check if all species in each taxon has the same domain composition
# They do not share the same domain composition but are very similar
# The mode can get me a representative of each taxon instead
bub1_consensus_taxon <- bub1_domains_filtered_taxons %>% 
  group_by(Taxon) %>%
  select(-Sequence, -Species_scientific, -`Species name`, -Prefix) %>%
  summarise_all(mean)
write.table(x = bub1_consensus_taxon, file = '../Results/Consensus_BUB1_domains_order.txt', 
            quote = F, sep = '\t', row.names = F, col.names = T)

bubr1_consensus_taxon <- bubr1_domains_filtered_taxons %>% 
  group_by(Taxon) %>%
  select(-Sequence, -Species_scientific, -`Species name`, -Prefix) %>%
  summarise_all(mean)
write.table(x = bubr1_consensus_taxon, file = '../Results/Consensus_BUBR1_domains_order.txt', 
            quote = F, sep = '\t', row.names = F, col.names = T)
