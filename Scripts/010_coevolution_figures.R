##################################################
####        010_coevolution_figures           ####
#### This script reads the output from the    ####
#### coevolution analyses and plots the data. ####
##################################################

#### Load libraries, set directory  ####

library(ggplot2)
library(reshape2)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggsignif)
library(cowplot)

setwd(<path_to_Evol_results>)

#### Define functions ####

get_MMI <- function(mut_info, domains){
  
  domain_names <- names(domains)
  domain_mean_mut_info <- as.data.frame(matrix(0, nrow = length(domain_names), ncol = length(domain_names)))
  rownames(domain_mean_mut_info) <- domain_names
  colnames(domain_mean_mut_info) <- domain_names
  
  # Fill the matrix with this loop
  counter_i <- 1
  for(domain_i in domains){
    counter_j <- 1
    for(domain_j in domains){
      
      # Get the mean mutual information between the two current domains
      current_mean <- mean(mut_info[domain_i, domain_j])
      
      domain_mean_mut_info[counter_i, counter_j] <- current_mean
      counter_j <- counter_j + 1
    }
    
    counter_i <- counter_i + 1
  }
  
  return(domain_mean_mut_info)
}

evol_heatmap <- function(mut_info, legend_title, domains, plot_title){
  
  # Get the mean mutual information for the domains
  domain_mean_mut_info <- get_MMI(mut_info, domains)
  
  #### Get the tidy data and the heatmap ####
  row_order <- hclust(dist(domain_mean_mut_info))$order
  col_order <- hclust(dist(t(domain_mean_mut_info)))$order
  new_data <- domain_mean_mut_info[row_order, col_order]
  new_tidy <- melt(new_data)
  
  new_tidy$var1 <- rep(colnames(new_data), length(colnames(new_data)))
  
  # To ensure that we get the same order, factor the variables
  new_tidy$var1 <- factor(new_tidy$var1, levels = unique(new_tidy$var1))
  new_tidy$variable <- factor(new_tidy$variable, levels = unique(new_tidy$variable))
  
  # Apply normalization by the domain vs itself
  norm_terms <- new_tidy %>% filter(variable == var1)
  new_tidy_norm <- left_join(x = new_tidy, y = norm_terms %>% select(var1, value), by = c("var1" = "var1"))
  new_tidy_norm %<>% mutate(norm = value.x / value.y)
  
  p <- ggplot(new_tidy_norm, aes(x = var1, y = variable)) +
    geom_tile(aes(fill = norm), colour = 'black', size = 1.1) +
    scale_fill_gradient2(name = legend_title, low = muted('blue'), mid = 'white', high = muted('red')) +
    xlab('') + ylab('') + 
    ggtitle(plot_title) + 
    theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 18), 
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_text(face = 'bold'), axis.title.y = element_text(face = 'bold'),
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
  
  output <- list()
  output$plot <- p
  output$data <- domain_mean_mut_info
  output$tidy_data <- new_tidy
  output$tidy_data_norm <- new_tidy_norm
  
  return(output)
  
}

#### Load the data ####

mut_info_corr_MDAT_BuBR1 <- as.matrix(read.table('BuBR1_common_seqs_APC.txt', h = F))

# Load info on domains
## Positions for the domains of the human BuBR1 sequence (ENSP00000287598)
domains_human <- list()

#### Annotate the domains ####
domains_human$KEN1 <- c(20:45)
domains_human$TPR <- c(60:201)
domains_human$`D-Box1` <- c(221:233)
domains_human$`D-box2` <- c(552:564)
domains_human$ABBA1 <- c(269:280)
domains_human$KEN2 <- c(292:321)
domains_human$ABBA2 <- c(337:348)
domains_human$GLEBS <- c(393:426)
domains_human$`Other ABBA` <- c(525:536)
domains_human$KARD <- c(665:686)
domains_human$CDII <- c(719:751)
domains_human$`Kinase-like` <- c(770:1016)

#### Plotting ####

# Run the function to get the plots
domain_mean_mut_info_corr_MDAT_BuBR1 <- evol_heatmap(mut_info_corr_MDAT_BuBR1, 'Normalized MIp', domains_human, 'BuBR1 domain co-evolution')

domain_mean_mut_info_corr_MDAT_BuBR1$plot

ggsave(filename = '../../Figures/Fig1B.pdf',
       plot = domain_mean_mut_info_corr_MDAT_BuBR1$plot,
       device = cairo_pdf, width = 10, height = 10, dpi = 500)


#### BUB1 ####

mut_info_corr_MDAT_BUB1 <- as.matrix(read.table('BUB1_common_seqs_APC.txt', h = F))

# Load info on domains

## Positions for the domains of the human BUB1 sequence (ENSP00000302530)
domains_human_BUB1 <- list()

domains_human_BUB1$TPR <- c(9:147)
domains_human_BUB1$GLEBS <- c(233:263)
domains_human_BUB1$CMI <- c(457:478)
domains_human_BUB1$ABBA1 <- c(524:535)
domains_human_BUB1$ABBA2 <- c(616:627) 
domains_human_BUB1$KARD <- c(650:671)
domains_human_BUB1$CDII <- c(740:772)
domains_human_BUB1$Kinase <- c(792:1051)

# Run the function to get the plots
domain_mean_mut_info_corr_MDAT_BUB1 <- evol_heatmap(mut_info_corr_MDAT_BUB1, 'Normalized MIp', domains_human_BUB1, 'BUB1 domain co-evolution')

domain_mean_mut_info_corr_MDAT_BUB1$plot

ggsave(filename = '../../Figures/FigSuppl1E.pdf',
       plot = domain_mean_mut_info_corr_MDAT_BUB1$plot,
       device = cairo_pdf, width = 10, height = 7, dpi = 500)

#### Select only the columns for KARD and kinase from the heatmaps ####

# BuBR1 #
BuBR1_KARD_kinase <- domain_mean_mut_info_corr_MDAT_BuBR1$tidy_data_norm %>%
  filter(var1 %in% c('KARD', 'Kinase-like'))

p_BUBR1 <- ggplot(BuBR1_KARD_kinase, aes(x = variable, y = var1)) +
  geom_tile(aes(fill = norm), colour = 'black', size = 1.1) +
  scale_fill_gradient2(name = 'Normalized MIp', low = muted('blue'), mid = 'white', high = muted('red')) +
  xlab('') + ylab('') + 
  ggtitle('BUBR1 co-evolution, KARD and kinase') + 
  theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 18), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold'), axis.title.y = element_text(face = 'bold'),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, size = 18, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.justification = 0.5,
        legend.title = element_text(size = 14),
        legend.text = element_text(angle = -90))
p_BUBR1


# BUB1 #
BUB1_KARD_kinase <- domain_mean_mut_info_corr_MDAT_BUB1$tidy_data_norm %>%
  filter(var1 %in% c('KARD', 'Kinase'))

p_BUB1 <- ggplot(BUB1_KARD_kinase, aes(x = variable, y = var1)) +
  geom_tile(aes(fill = norm), colour = 'black', size = 1.1) +
  scale_fill_gradient2(name = 'Normalized MIp', low = muted('blue'), mid = 'white', high = muted('red')) +
  xlab('') + ylab('') + 
  ggtitle('BUB1 co-evolution, KARD and kinase') + 
  theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 18), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold'), axis.title.y = element_text(face = 'bold'),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, size = 18, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.justification = 0.5,
        legend.title = element_text(size = 14),
        legend.text = element_text(angle = -90))
p_BUB1

p_KARD_kinase <- plot_grid(p_BUBR1, p_BUB1, nrow = 2)

ggsave(filename = '../../Figures/Fig1B_only_KARD_kinase.pdf',
       plot = p_KARD_kinase,
       device = cairo_pdf, width = 10, height = 10, dpi = 500)


#### Look at the entropy distributions ####

# Load everything for BuBR1
entropy_MDAT_BuBR1 <- read.table('BuBR1_common_seqs_entropy.txt', h = F)
entropy_MDAT_BuBR1$position <- 1:length(entropy_MDAT_BuBR1$V1)
colnames(entropy_MDAT_BuBR1)[1] <- 'entropy'
sorted_entropy_MDAT_BuBR1 <- entropy_MDAT_BuBR1[order(entropy_MDAT_BuBR1$entropy, decreasing = TRUE),]
sorted_entropy_MDAT_BuBR1$domain <- ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$`D-Box1`, 'D-Box1',
                                        ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$`D-Box2`, 'D-Box2',   
                                          ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$KARD, 'KARD',
                                              ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$KEN1, 'KEN1',
                                                     ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$KEN2, 'KEN2',
                                                            ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$TPR, 'TPR',
                                                                   ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$`Kinase-like`, 'kinase-like',
                                                                          ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$GLEBS, 'GLEBS',
                                                                                 ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$ABBA1, 'ABBA1',
                                                                                        ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$`Other ABBA`, 'Other ABBA',
                                                                                               ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$ABBA2, 'ABBA2',
                                                                                                      ifelse(sorted_entropy_MDAT_BuBR1$position %in% domains_human$CDII, 'CDII',
                                                                                                             'Non-domain'))))))))))))

sorted_entropy_MDAT_BuBR1$rank <- 1:length(entropy_MDAT_BuBR1$entropy)
sorted_entropy_MDAT_BuBR1$protein <- 'BuBR1'

# Load everything for BUB1
entropy_MDAT_BUB1 <- read.table('BUB1_common_seqs_entropy.txt', h = F)
entropy_MDAT_BUB1$position <- 1:length(entropy_MDAT_BUB1$V1)
colnames(entropy_MDAT_BUB1)[1] <- 'entropy'
sorted_entropy_MDAT_BUB1 <- entropy_MDAT_BUB1[order(entropy_MDAT_BUB1$entropy, decreasing = TRUE),]
sorted_entropy_MDAT_BUB1$domain <- ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$KARD, 'KARD',
                                          ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$TPR, 'TPR',
                                                 ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$Kinase, 'Kinase',
                                                        ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$GLEBS, 'GLEBS',
                                                               ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$ABBA1, 'ABBA1',
                                                                      ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$ABBA2, 'ABBA2',
                                                                        ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$CMI, 'CMI',
                                                                               ifelse(sorted_entropy_MDAT_BUB1$position %in% domains_human_BUB1$CDII, 'CDII',
                                                                                      'Non-domain'))))))))

sorted_entropy_MDAT_BUB1$rank <- 1:length(entropy_MDAT_BUB1$entropy)
sorted_entropy_MDAT_BUB1$protein <- 'BUB1'

# Put the two datasets together and plot the distributions of entropies side by side
all_entropies <- rbind(sorted_entropy_MDAT_BUB1, sorted_entropy_MDAT_BuBR1)

all_entropies %<>% 
  mutate(key_res = ifelse(and(protein == 'BUB1', position %in% c(655, 661, 665, 821, 917, 946)),
                          0, 
                          ifelse(and(protein == 'BuBR1', position %in% c(670, 676, 680, 795, 882, 911)), 
                                 0, 
                                 ifelse(protein == 'BUB1', 1, 2)))) %>%
  mutate(key_res = as.factor(key_res)) %>%
  mutate(label = ifelse(and(protein == 'BUB1', position %in% c(655, 661, 665, 821, 917, 946)),
                        position, 
                        ifelse(and(protein == 'BuBR1', position %in% c(670, 676, 680, 795, 882, 911)), 
                               position, 
                               NA)))

#### Plot the entropy distributions of the common domains ####

plot_entropies2 <- all_entropies %>% 
  mutate(protein = gsub(pattern = 'BuBR1', replacement = 'BUBR1', x = protein)) %>%
  mutate(domain = gsub(pattern = 'kinase-like', replacement = 'Kinase', x = domain)) %>%
  filter(domain %in% c('ABBA1', 'ABBA2', 'CDII', 'GLEBS', 'KARD', 'Kinase', 'TPR', 'Non-domain')) %>%
  ggplot(aes(x = factor(domain, levels = c('ABBA1', 'ABBA2', 'CDII', 'GLEBS', 'KARD', 'Kinase', 'TPR', 'Non-domain')), 
             y = entropy, fill = protein)) +
  geom_point(aes(colour = key_res, group = protein, shape = key_res, alpha = key_res, size = key_res),
             position = position_jitterdodge(jitter.width = 0.45)) +
  scale_colour_manual(values = c('black', '#ff4d4d', '#3399ff'), guide = 'none') +
  scale_shape_manual(values = c(15, 16, 16), guide = 'none') +
  scale_alpha_manual(values = c(1, 0.5, 0.5), guide = 'none') +
  scale_size_manual(values = c(2, 1.5, 1.5), guide = 'none') +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  scale_fill_manual(values = c('#ff4d4d', '#3399ff')) +
  xlab('Domain') + ylab('Entropy') + 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold'), axis.title.y = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
        axis.line = element_line()) +
  stat_compare_means(aes(group = protein), method = 't.test', paired = FALSE, label = 'p.format') +
  labs(fill = '') + ylim(0,2.25)

plot_entropies2

ggsave(filename = '../../Figures/Fig1C.pdf',
       plot = plot_entropies2,
       device = cairo_pdf, width = 10, height = 7, dpi = 500)

#### Work with the merged sequences ####

mut_info_corr_MDAT_merged <- as.matrix(read.table('merged_common_seqs_APC.txt', h = F))

domains_human <- list()

#### Annotate the domains ####

domains_human$BUB1.TPR <- c(9:147)
domains_human$BUB1.GLEBS <- c(233:263)
domains_human$BUB1.CMI <- c(457:478)
domains_human$BUB1.ABBA1 <- c(524:535)
domains_human$BUB1.ABBA2 <- c(616:627) 
domains_human$BUB1.KARD <- c(650:671)
domains_human$BUB1.CDII <- c(740:772)
domains_human$BUB1.Kinase <- c(792:1051)

domains_human$BuBR1.KEN1 <- c(20:45) + 1085
domains_human$BuBR1.TPR <- c(60:201) + 1085
domains_human$`BuBR1.D-Box1` <- c(221:233) + 1085
domains_human$`BuBR1.D-Box2` <- c(552:564) + 1085
domains_human$BuBR1.ABBA1 <- c(269:280) + 1085
domains_human$BuBR1.KEN2 <- c(292:321) + 1085
domains_human$BuBR1.ABBA2 <- c(337:348) + 1085
domains_human$BuBR1.GLEBS <- c(393:426) + 1085
domains_human$`BuBR1.Other ABBA` <- c(525:536) + 1085
domains_human$BuBR1.KARD <- c(665:686) + 1085
domains_human$BuBR1.CDII <- c(719:751) + 1085
domains_human$`BuBR1.Kinase-like` <- c(770:1016) + 1085

domain_mean_mut_info_corr_MDAT_merged <- evol_heatmap(mut_info_corr_MDAT_merged, 'Normalized MIp', domains_human, 'Merged BUB1-BuBR1')

# Prepare the axes labels and their colors
protein_colors <- sapply(domain_mean_mut_info_corr_MDAT_merged$tidy_data$var1,
                         function(x) ifelse(substr(x, 1, 5) == 'BuBR1', return('#3399ff'), return('#ff4d4d'))) 

ordered_labels <- as.character(levels(domain_mean_mut_info_corr_MDAT_merged$tidy_data_norm$var1))
domain_split <- strsplit(ordered_labels, split = '[.]')
new_names <- unlist(lapply(domain_split, function(x) return(x[2])))

final_plot <- domain_mean_mut_info_corr_MDAT_merged$plot +
  theme(axis.text.x = element_text(colour = protein_colors),
        axis.text.y = element_text(colour = protein_colors)) + 
  scale_x_discrete(labels = new_names) +
  scale_y_discrete(labels = new_names)

ggsave(filename = '../../Figures/FigSuppl6D.pdf',
       plot = final_plot,
       device = cairo_pdf, width = 10, height = 10, dpi = 500)

#### Select only the columns for KARD and kinase from the heatmap ####

KARD_kinase <- domain_mean_mut_info_corr_MDAT_merged$tidy_data_norm %>%
  filter(var1 %in% c('BUB1.KARD', 'BuBR1.KARD', 'BUB1.Kinase', 'BuBR1.Kinase-like'))

KARD_kinase$var1 <- factor(KARD_kinase$var1, levels = c('BUB1.KARD', 'BuBR1.KARD', 'BUB1.Kinase', 'BuBR1.Kinase-like'))

# Prepare labels and their colors
ordered_labels_y <- as.character(levels(KARD_kinase$var1))
domain_split_y <- strsplit(ordered_labels_y, split = '[.]')
new_names_y <- unlist(lapply(domain_split_y, function(x) return(x[2])))
new_colors <- c('#ff4d4d', '#3399ff', '#ff4d4d', '#3399ff')

p_merged <- ggplot(KARD_kinase, aes(x = variable, y = var1)) +
  geom_tile(aes(fill = norm), colour = 'black', size = 1.1) +
  scale_fill_gradient2(name = 'Normalized MIp', low = muted('blue'), mid = 'white', high = muted('red')) +
  xlab('') + ylab('') + 
  ggtitle('') + 
  theme(plot.title = element_text(hjust = 0.5, face = 'plain', size = 18), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold'), axis.title.y = element_text(face = 'bold'),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = -90, size = 18, hjust = 0, vjust = 0.5, colour = protein_colors),
        axis.text.y = element_text(size = 18, colour = new_colors),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.justification = 0.5,
        legend.title = element_text(size = 14),
        legend.text = element_text(angle = -90)) + 
  scale_x_discrete(labels = new_names) +
  scale_y_discrete(labels = new_names_y)
p_merged

ggsave(filename = '../../Figures/FigSuppl6D_only_KARD_kinase.pdf',
       plot = p_merged,
       device = cairo_pdf, width = 10, height = 5, dpi = 500)
