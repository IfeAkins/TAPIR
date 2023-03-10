knitr::opts_chunk$set(echo = TRUE)


library(tidyverse)
library(magrittr)
library(ggthemes)
library(scales)
library(patchwork)



kraken_enriched_extra_G230302_df <- kraken_G230229_modified_df %>% mutate(week = 'week8E') %>% 
  bind_rows(kraken_G230302_modified_df %>% mutate(week = 'week8'))


plot_bar_kraken_enriched_extra <- function(kraken_combined_report_final_df, species_colors){
  library(tidyverse)
  library(magrittr)
  library(ggthemes)
  library(scales)
  library(viridis)
  
  # minor rearrangement
  
  # rearrange sample and species name
  # sample names
  sample_name_elements <- kraken_combined_report_final_df %>% filter(!Sample %in% c('Swab', 'swab', 'Water', 'water')) %>% 
    distinct(Sample) %>% pull() %>% unique() %>% sort()
  
  sample_name_elements <- c('Swab', 'Water') %>% append(sample_name_elements)
  
  # species
  kraken_combined_report_final_species <- kraken_combined_report_final_df %>% filter(!species %in% c('unclassified', 'unassigned', 'Others')) %>% 
    distinct(species) %>% pull() %>% unique() %>% sort()
  
  kraken_combined_report_final_species %<>% c('Others', 'unclassified', 'unassigned')
  
  # reorder species
  kraken_combined_report_final_df$species <- factor(kraken_combined_report_final_df$species, levels = kraken_combined_report_final_species)
  # reorder week
  # kraken_combined_report_final_df$week <- factor(kraken_combined_report_final_df$week, levels = c('week8', 'week8E'))
  
  # plot
  kraken_combined_report_final_plot <- kraken_combined_report_final_df %>% 
    ggplot(aes(x = factor(Sample, levels = sample_name_elements), y = percentage, fill = species, label = Read_category)) + 
    geom_col(position = position_stack(reverse = TRUE)) +
    xlab('Sample') + ylab("Proportion (Species)") + labs(fill = 'Species') + 
    facet_grid(week ~ factor(Site, levels = c("Neg_Ctrl", "Anus", "Nose")), scales = 'free_x', space='free_x') +
    theme_bw() +
    geom_text(position = position_stack(vjust = 0.5, reverse = TRUE), color = "white") +
    scale_fill_manual(breaks = kraken_combined_report_final_species, values = species_colors) + 
    theme(
      axis.title.y = element_text(size = 20, face = 'bold'),
      legend.position = "bottom",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      axis.title.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      strip.text.x = element_text(size = 16, face = 'bold'),
      strip.text.y = element_text(size = 16, face = 'bold'),
      strip.background = element_rect(colour="black",
                                      fill="white"),
      axis.text.y = element_text(size = 25, face = 'bold'),
      axis.text.x = element_text(size = 14, face = 'bold', angle = 45, hjust = 1))
  
  return(kraken_combined_report_final_plot)
}


kraken_enriched_extra_G230307_plot <- plot_bar_kraken_enriched_extra(kraken_enriched_extra_G230302_df, kraken_G230229_species_name_and_colour)
kraken_enriched_extra_G230307_plot


jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230302_kraken_enriched_extra_plot.jpeg", width = 12000, height = 7000, units = 'px', res = 600)
print(kraken_enriched_extra_G230307_plot)
dev.off()

kraken_enriched_extra_G230302_df %>% 
ggplot(aes(x = Sample, y = percentage, fill = species, label = Read_category)) + 
  geom_col(position = position_stack(reverse = TRUE)) +
  xlab('Sample') + ylab("Proportion (Species)") + labs(fill = 'Species') + 
  facet_grid(week ~ factor(Site, levels = c("Neg_Ctrl", "Anus", "Nose")), scales = 'free_x', space='free_x') +
  theme_bw() +
  geom_text(position = position_stack(vjust = 0.5, reverse = TRUE), color = "white")
