##################################################
####    007_domain_copresence                 ####
#### This script looks at table of presence   ####
#### of domains in each sequence and produces ####
#### figures on domain composition.           ####
##################################################

#### Load libraries ####

library(tidyverse)
library(magrittr)
library(ggplot2)
library(Cairo)
library(UpSetR)
library(scales)

# Set working directory to the scripts folder
setwd(<path_to_scripts>)

#### BuBR1 ####

# Remove the pseudo_kinase and X_pseudo_kinase assignments
domains_BuBR1 <- read.table('../Data/Domain_annotation/BuBR1/domain_table_BuBR1.txt',
                            h = T, sep = '\t') %>%
  select(-X_pseudo_kinase, -pseudo_kinase)

domain_counts <- domains_BuBR1 %>% select(-Sequence) %>% colSums()
# ABBA1      ABBA2 ABBA_other       CDII       KEN1       KEN2        TPR      GLEBS       KARD       Dbox 
#   140        146         85         96        152        146        157        148         92        149 
# CMI     kinase      MadaM  Other_KEN 
#   4         93         11          1 

names(domain_counts)[3] <- "Other ABBA"
names(domain_counts)[10] <- "D-Box"
names(domain_counts)[12] <- "Kinase-like"
names(domain_counts)[14] <- "Other KEN"

new_labels <- c('Sequence')
for(i in 1:length(domain_counts)){
  new_label <- paste(names(domain_counts)[i], ' (', as.character(domain_counts[i]), ')', sep = '')
  new_labels <- c(new_labels, new_label)
}

# Replace column names in the original table
colnames(domains_BuBR1) <- new_labels

# Following Tromer et al. Open Biol 2016, look at correlation in domain presence
correlation_table_BuBR1 <- domains_BuBR1 %>% select(-Sequence) %>% cor()

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
  scale_fill_gradient2(name = 'Pearson correlation', low = muted('blue'), mid = 'white', high = muted('red')) +
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
        legend.text = element_text(angle = -90))
p

# Save the heatmap
ggsave(filename = '../Figures/FigSuppl1A.pdf',
       plot = p, device = cairo_pdf, width = 10, height = 10, dpi = 500)

#### BUB1 ####

domains_BUB1 <- read.table('../Data/Domain_annotation/BUB1/domain_table_BUB1.txt',
                           h = T, sep = '\t')

domain_counts_BUB1 <- domains_BUB1 %>% select(-Sequence) %>% colSums()
# BUB1
# ABBA2 ABBA_other       CDII        TPR     kinase      GLEBS        CMI       Dbox       KARD      MadaM 
#     1        123        163        139        164        125        122          4        104          1 

names(domain_counts_BUB1)[2] <- "Other ABBA"
names(domain_counts_BUB1)[8] <- "D-Box"
names(domain_counts_BUB1)[5] <- "Kinase"

new_labels <- c('Sequence')
for(i in 1:length(domain_counts_BUB1)){
  new_label <- paste(names(domain_counts_BUB1)[i], ' (', as.character(domain_counts_BUB1[i]), ')', sep = '')
  new_labels <- c(new_labels, new_label)
}

# Replace column names in the original table
colnames(domains_BUB1) <- new_labels

# Calculate correlations between BUB1 domains
correlation_table_BUB1 <- domains_BUB1 %>% select(-Sequence) %>% cor()

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
  scale_fill_gradient2(name = 'Pearson correlation', low = muted('blue'), mid = 'white', high = muted('red')) +
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
        legend.text = element_text(angle = -90))
p

# Save the heatmap
ggsave(filename = '../Figures/FigSuppl1F.pdf',
       plot = p, device = cairo_pdf, width = 10, height = 7, dpi = 500)

#### Show the most abundant domain compositions with UpSetR ####

# Save the plot for BuBR1
ggsave(filename = '../Figures/FigSuppl1B.pdf', 
       device = cairo_pdf, width = 10, height = 7, dpi = 500, 
       plot = upset(domains_BuBR1, order.by = "freq", nsets = 16))

# Save the plot for BUB1
ggsave(filename = '../Figures/FigSuppl_BUB1_composition.pdf',
       device = cairo_pdf, width = 10, height = 7, dpi = 500, 
       plot = upset(domains_BUB1, order.by = "freq", nsets = 16))

