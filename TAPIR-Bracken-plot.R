library(tidyverse)
library(magrittr)
library(ggthemes)
library(scales)
library(kableExtra)
library(viridis)


bracken_G230302 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/combined_bracken_S_G230302.txt')

head(bracken_G230302) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")


current_species_name_and_colours <- readxl::read_xlsx('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/colors_by_gram_230303.xlsx')

head(current_species_name_and_colours) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")


bracken_parse_df <- function(bracken_combined_report){
  library(tidyverse)
  library(magrittr)
  library(ggthemes)
  library(scales)
  
  # parse kraken combined report
  bracken_combined_report_df <- bracken_combined_report %>% 
    mutate(percentage_total_reads = fraction_total_reads * 100) %>% 
    mutate(Site = case_when(str_detect(Filename, '[Ss]wab|[Ww]ater') ~ 'Neg_Ctrl',
                            str_detect(Filename, 'A') ~ 'Anus',
                            str_detect(Filename, 'B|N') ~ 'Nose')) %>% 
    mutate(Sample = str_extract(Filename, "^[a-zA-Z0-9]+")) %>% 
    mutate(Read_category = case_when(new_est_reads < 1000 ~ "",
                                     new_est_reads >= 1000 & new_est_reads <= 10000 ~ "*",
                                     new_est_reads > 10000 ~ "**")) %>% 
    mutate(name = if_else(percentage_total_reads >= 1, name, 'Others')) %>% 
    mutate(name = if_else(!str_detect(name, 'virus|phage'), name, 'Others'))
  
  # calculate Others
  capture_others_df <- bracken_combined_report_df %>% 
    filter(name == 'Others') %>%
    group_by(Filename) %>%
    mutate(sumpercentage = sum(percentage_total_reads),
           sumkrakenassignedreads = sum(kraken_assigned_reads),
           sumaddedreads = sum(added_reads),
           sumnewestreads = sum(new_est_reads),
           sumfractiontotalreads = sum(fraction_total_reads),
           sumpercentagetotalreads = sum(percentage_total_reads)) %>%
    ungroup() %>% 
    select(-c(percentage_total_reads, kraken_assigned_reads, added_reads, 
              new_est_reads, taxonomy_id, fraction_total_reads, 
              percentage_total_reads, Read_category)) %>% 
    distinct() %>% 
    rename(kraken_assigned_reads = sumkrakenassignedreads,
           added_reads = sumaddedreads,
           new_est_reads = sumnewestreads,
           fraction_total_reads = sumfractiontotalreads,
           percentage_total_reads = sumpercentagetotalreads) %>% 
    mutate(taxonomy_id = 0) %>% 
    mutate(Read_category = case_when(new_est_reads < 1000 ~ "",
                                     new_est_reads >= 1000 & new_est_reads <= 10000 ~ "*",
                                     new_est_reads > 10000 ~ "**")) %>% 
    distinct() %>% 
    select(Filename, name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, 
           fraction_total_reads, percentage_total_reads, Site, Sample, Read_category)
  
  # merge both dataframes
  bracken_combined_report_semifinal_df <- bracken_combined_report_df %>%
    filter(name != 'Others') %>% 
    dplyr::bind_rows(capture_others_df)
  
  # capture unassigned
  capture_unassigned_df <- bracken_combined_report_semifinal_df %>% group_by(Filename) %>%
    mutate(sumkrakenassignedreads = sum(kraken_assigned_reads),
           sumaddedreads = sum(added_reads),
           sumnewestreads = sum(new_est_reads),
           sumfractiontotalreads = sum(fraction_total_reads),
           sumpercentagetotalreads = sum(percentage_total_reads)) %>% 
    select(-c(name, percentage_total_reads, kraken_assigned_reads, added_reads, 
              new_est_reads, taxonomy_id, fraction_total_reads, 
              percentage_total_reads, Read_category)) %>% 
    distinct() %>% 
    mutate(kraken_assigned_reads = round(sumkrakenassignedreads/sumpercentagetotalreads * 100 - sumkrakenassignedreads),
           added_reads = round(sumaddedreads/sumpercentagetotalreads * 100 - sumaddedreads),
           new_est_reads = round(sumnewestreads/sumpercentagetotalreads * 100 - sumnewestreads),
           fraction_total_reads = round(sumfractiontotalreads/sumpercentagetotalreads * 100 - sumfractiontotalreads),
           percentage_total_reads = 100 - sumpercentagetotalreads,
           taxonomy_id = 0,
           taxonomy_lvl = 'U',
           name = 'unassigned') %>% 
    distinct() %>% 
    mutate(Read_category = case_when(new_est_reads < 1000 ~ "",
                                     new_est_reads >= 1000 & new_est_reads <= 10000 ~ "*",
                                     new_est_reads > 10000 ~ "**")) %>% 
    distinct() %>% ungroup() %>%
    select(Filename, name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, 
           fraction_total_reads, percentage_total_reads, Site, Sample, Read_category)
  
  # merge both dataframes
  bracken_combined_report_final_df <- bracken_combined_report_semifinal_df %>%
    dplyr::bind_rows(capture_unassigned_df)
  
  return(bracken_combined_report_final_df)
}



bracken_G230302_modified_df <- bracken_parse_df(bracken_G230302)

head(bracken_G230302_modified_df) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")



bracken_G230302_species <- bracken_G230302_modified_df %>% 
  filter(!name %in%  c('Others', 'unassigned')) %>% distinct(name) %>% 
  pull() %>% unique() %>% sort() %>% c('Others', 'unassigned')

print(bracken_G230302_species)


setdiff(bracken_G230302_species, current_species_name_and_colours$species_name)



scales::viridis_pal()(length(unique(bracken_G230302_modified_df$name)))
# I chose the first one




# filter only species of interest for the plot and convert dataframe to vector
bracken_G230302_species_name_and_colour <- current_species_name_and_colours %>% 
  filter(species_name %in% c(bracken_G230302_species)) %>% dplyr::pull(species_colour, species_name)



plot_bar_bracken <- function(bracken_combined_report_final_df, species_colors){
  library(tidyverse)
  library(magrittr)
  library(ggthemes)
  library(scales)
  library(viridis)
  
  # minor rearrangement
  
  # rearrange sample and species name
  # sample names
  sample_name_elements <- bracken_combined_report_final_df %>% 
    filter(!Sample %in% c('Swab', 'swab', 'Water', 'water')) %>% 
    distinct(Sample) %>% pull() %>% unique() %>% sort()
  
  sample_name_elements <- c('Swab', 'Water') %>% append(sample_name_elements)
  
  # species
  bracken_combined_report_final_species <- bracken_combined_report_final_df %>% 
    filter(!name %in% c('Others', 'unassigned')) %>% 
    distinct(name) %>% pull() %>% unique() %>% 
    sort() %>% c('Others', 'unassigned')
  
  
  # reorder species
  bracken_combined_report_final_df$name <- factor(bracken_combined_report_final_df$name, levels = bracken_combined_report_final_species)
  
  # plot
  bracken_combined_report_final_plot <- bracken_combined_report_final_df %>% 
    ggplot(aes(x = factor(Sample, levels = sample_name_elements), y = percentage_total_reads, fill = name, label = Read_category)) + 
    geom_col(position = position_stack(reverse = TRUE)) +
    xlab('Sample') + ylab("Proportion (Species)") + labs(fill = 'Species') + 
    facet_grid(. ~ factor(Site, levels = c("Neg_Ctrl", "Anus", "Nose")), scales = 'free_x', space='free_x') +
    theme_bw() +
    geom_text(position = position_stack(vjust = 0.5, reverse = TRUE), color = "white") +
    scale_fill_manual(breaks = bracken_combined_report_final_species, values = species_colors) + 
    theme(
      axis.title.y = element_text(size = 20, face = 'bold'),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
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
  
  return(bracken_combined_report_final_plot)
}


bracken_G230302_plot <- plot_bar_bracken(bracken_G230302_modified_df, bracken_G230302_species_name_and_colour)
bracken_G230302_plot


write.table(bracken_G230302_modified_df, 'Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_bracken_per_sequence_date/G230302_combined_parsed_bracken.tsv', row.names = FALSE, sep = '\t')

jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230302_bracken_plot.jpeg", width = 10000, height = 5000, units = 'px', res = 600)
print(bracken_G230302_plot)
dev.off()

