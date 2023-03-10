knitr::opts_chunk$set(echo = TRUE)


library(tidyverse)
library(magrittr)
library(ggthemes)
library(scales)
library(patchwork)



kraken_bracken_G230302_df <- kraken_G230302_modified_df %>% mutate(type = 'KRAKEN') %>% 
  bind_rows(bracken_G230302_modified_df %>% 
              mutate(type = 'BRACKEN',
                     percentage = percentage_total_reads,
                     cladeReads = kraken_assigned_reads,
                     taxonReads = new_est_reads,
                     samplename = Filename,
                     taxRank = taxonomy_lvl,
                     taxID = taxonomy_id,
                     species = name)) %>% 
  select(-c(percentage_total_reads, kraken_assigned_reads, new_est_reads, Filename, 
            taxonomy_lvl, taxonomy_id, name, added_reads,
            fraction_total_reads))


plot_bar_kraken_bracken <- function(kraken_combined_report_final_df, species_colors){
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
  # reorder type
  kraken_combined_report_final_df$type <- factor(kraken_combined_report_final_df$type, levels = c('KRAKEN', 'BRACKEN'))
  
  # plot
  kraken_combined_report_final_plot <- kraken_combined_report_final_df %>% 
    ggplot(aes(x = factor(Sample, levels = sample_name_elements), y = percentage, fill = species, label = Read_category)) + 
    geom_col(position = position_stack(reverse = TRUE)) +
    xlab('Sample') + ylab("Proportion (Species)") + labs(fill = 'Species') + 
    facet_grid(type ~ factor(Site, levels = c("Neg_Ctrl", "Anus", "Nose")), scales = 'free_x', space='free_x') +
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


kraken_bracken_G230302_plot <- plot_bar_kraken_bracken(kraken_bracken_G230302_df, kraken_G230302_species_name_and_colour)
kraken_bracken_G230302_plot


jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230302_kraken_bracken_plot.jpeg", width = 12000, height = 7000, units = 'px', res = 600)
print(kraken_bracken_G230302_plot)
dev.off()

