knitr::opts_chunk$set(echo = TRUE)

#First load required libraries. If you don't have these yet, kindly install them, so as to avoid errors downstream.
library(tidyverse)
library(magrittr)
library(ggthemes)
library(scales)
library(kableExtra)
library(viridis)

#Read in the data (combined kraken report file) and combined kraken unclassified file

kraken_G230317 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_kraken_report_per_sequence_date/combined_kraken_report_S_G230317.txt')
kraken_G230317_unclassified <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_kraken_report_per_sequence_date/combined_kraken_report_U_G230317.txt')

# merge
kraken_G230317_merge_df <- kraken_G230317 %>% bind_rows(kraken_G230317_unclassified)

head(kraken_G230317_merge_df) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")


#Read in data containing existing species and colours from Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/MicrobialSpeciesNameAndColour.tsv.

current_species_name_and_colours <- readxl::read_xlsx('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/colors_by_gram_230303.xlsx')

head(current_species_name_and_colours) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")


#Load the customized function below. The function helps to parse the combined kraken report file.

kraken3_parse_df <- function(kraken_combined_report){
  library(tidyverse)
  library(magrittr)
  library(ggthemes)
  library(scales)
  
  # parse kraken combined report
  kraken_combined_report_df <- kraken_combined_report %>% 
    mutate(Site = case_when(str_detect(samplename, '[Ss]wab|[Ww]ater') ~ 'Neg_Ctrl',
                            str_detect(samplename, 'A') ~ 'Anus',
                            str_detect(samplename, 'B|N') ~ 'Nose')) %>% 
    mutate(Sample = str_extract(samplename, "^[a-zA-Z0-9]+")) %>% 
    mutate(Read_category = case_when(taxonReads < 1000 ~ "",
                                     taxonReads >= 1000 & taxonReads <= 10000 ~ "*",
                                     taxonReads > 10000 ~ "**")) %>% 
    mutate(species = if_else(percentage > 1 | species == 'unclassified', species, 'Others')) %>% 
    mutate(species = if_else(!str_detect(species, 'virus|phage'), species, 'Others'))
  
  # calculate 'Others'
  capture_others_df <- kraken_combined_report_df %>% 
    filter(species == 'Others') %>%
    group_by(samplename) %>%
    mutate(sumpercentage = sum(percentage),
           sumcladereads = sum(cladeReads),
           sumtaxonReads = sum(taxonReads)) %>%
    ungroup() %>% 
    select(-c(percentage, cladeReads, taxonReads, taxRank, taxID, Read_category)) %>% 
    distinct() %>% 
    rename(percentage = sumpercentage,
           cladeReads = sumcladereads,
           taxonReads = sumtaxonReads) %>% 
    mutate(taxID = 0,
           taxRank = 'None') %>% 
    mutate(Read_category = case_when(taxonReads < 1000 ~ "",
                                     taxonReads >= 1000 & taxonReads <= 10000 ~ "*",
                                     taxonReads > 10000 ~ "**")) %>% 
    distinct() %>% 
    select(samplename, percentage, cladeReads, taxonReads, taxRank, taxID, 
           species, Site, Sample, Read_category)
  
  # merge both dataframes
  kraken_combined_report_semifinal_df <- kraken_combined_report_df %>%
    filter(species != 'Others') %>% 
    dplyr::bind_rows(capture_others_df)
  
  # capture unassigned
  capture_unassigned_df <- kraken_combined_report_semifinal_df %>% group_by(samplename) %>%
    mutate(sumpercentage = sum(percentage),
           sumcladereads = sum(cladeReads),
           sumtaxonReads = sum(taxonReads)) %>% 
    select(-c(species, percentage, cladeReads, taxonReads, taxRank, taxID, Read_category)) %>% 
    distinct() %>% 
    mutate(taxx = round(sumtaxonReads/sumpercentage * 100 - sumtaxonReads),
           clad = round(sumcladereads/sumpercentage * 100 - sumcladereads),
           unassigned = 100 - sumpercentage) %>% 
    distinct() %>% 
    rename(percentage = unassigned,
           cladeReads = clad,
           taxonReads = taxx) %>% 
    mutate(species = 'unassigned',
           taxID = 0,
           taxRank = 'None') %>% 
    mutate(Read_category = case_when(taxonReads < 1000 ~ "",
                                     taxonReads >= 1000 & taxonReads <= 10000 ~ "*",
                                     taxonReads > 10000 ~ "**")) %>% 
    distinct() %>% ungroup() %>% 
    select(samplename, percentage, cladeReads, taxonReads, taxRank, taxID, 
           species, Site, Sample, Read_category)
  
  # merge both dataframes
  kraken_combined_report_final_df <- kraken_combined_report_semifinal_df %>%
    dplyr::bind_rows(capture_unassigned_df)
  
  return(kraken_combined_report_final_df)
}


#Call the function in order to generate a modified kraken report file containing relevant metadata (e.g Site)
kraken_G230317_modified_df <- kraken3_parse_df(kraken_G230317_merge_df)

head(kraken_G230317_modified_df) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")

#Show species contained in modified kraken report
kraken_G230317_species <- kraken_G230317_modified_df %>% 
  filter(!species %in% c('unclassified', 'unassigned', 'Others')) %>% distinct(species) %>% 
  pull() %>% unique() %>% sort()

kraken_G230317_species %<>% c('unclassified', 'unassigned', 'Others')

print(kraken_G230317_species)

print(current_species_name_and_colours$species_name)

#Identify species in the current df not present in the list of species
setdiff(kraken_G230317_species, current_species_name_and_colours$species_name)

#Assign colors to the new species. You can search for suitable colors online,
#bearing in mind that color-blind-friendly colours are preferable. 
#Ensure that colours are unique to each species. Check whether or not colour chosen is already assigned to a microbial species.

#new_species_and_colours <- c("unassigned" = "#778899")

#Create a dataframe of the new species and colours, merge these with existing dataframe of species and colours
# change to dataframe and maintain column header names like the existing df
#new_species_and_colours_df <- data.frame(species_name = new_species_and_colours %>% names(),
                                        #species_colour = new_species_and_colours %>% unique())

# merge with existing current_species_names_and_colours_df or the data you read from Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/MicrobialSpeciesNameAndColour.tsv
# exclude 'unclassified' and 'Others' and then merge again
#current_species_names_and_colours_minus_unclassifiedandothers_df <- current_species_name_and_colours %>%
  #bind_rows(new_species_and_colours_df) %>%
  #arrange(species_name) %>% 
  #filter(!species_name %in% c('unclassified', 'unassigned', 'Others'))

#current_species_names_and_colours_only_unclassifiedandothers_df <- current_species_name_and_colours %>% 
  #bind_rows(new_species_and_colours_df) %>%
  #arrange(species_name) %>% 
  #filter(species_name %in% c('unclassified', 'unassigned', 'Others'))

#current_species_names_and_colours_update1_df <- current_species_names_and_colours_minus_unclassifiedandothers_df %>%
  #bind_rows(current_species_names_and_colours_only_unclassifiedandothers_df) %>% distinct()

current_species_name_and_colours_df <- as.data.frame(current_species_name_and_colours)

#Filter for only the species present in the combined kraken output
# filter only species of interest for the plot and convert dataframe to vector
kraken_G230317_species_name_and_colour <- current_species_name_and_colours %>% 
  filter(species_name %in% c(kraken_G230317_species)) %>% dplyr::pull(species_colour, species_name)

#This function helps in plotting bar charts.
plot_bar_kraken <- function(kraken, species_colors){
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
  
  # plot
  kraken_combined_report_final_plot <- kraken_combined_report_final_df %>% 
    ggplot(aes(x = factor(Sample, levels = sample_name_elements), y = percentage, fill = species, label = Read_category)) + 
    geom_col(position = position_stack(reverse = TRUE)) +
    xlab('Sample') + ylab("Proportion (Species)") + labs(fill = 'Species') + 
    facet_grid(. ~ factor(Site, levels = c("Neg_Ctrl", "Anus", "Nose")), scales = 'free_x', space='free_x') +
    theme_bw() +
    geom_text(position = position_stack(vjust = 0.5, reverse = TRUE), color = "white") +
    scale_fill_manual(breaks = kraken_combined_report_final_species, values = species_colors) + 
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
  
  return(kraken_combined_report_final_plot)
}


kraken_G230317_modified_df <- kraken3_parse_df(kraken_G230317_merge_df)

head(kraken_G230317_modified_df) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")


#Show species contained in modified kraken report
kraken_G230317_species <- kraken_G230317_modified_df %>% 
  filter(!species %in% c('unclassified', 'unassigned', 'Others')) %>% distinct(species) %>% 
  pull() %>% unique() %>% sort()

kraken_G230317_species %<>% c('unclassified', 'unassigned', 'Others')

print(kraken_G230317_species)

#Identify species in the current df not present in the list of species
setdiff(kraken_G230317_species, current_species_name_and_colours$species_name)
                      

new_species_and_colours <- c("unassigned" = "#778899")

# change to dataframe and maintain column header names like the existing df
new_species_and_colours_df <- data.frame(species_name = new_species_and_colours %>% names(),
                                         species_colour = new_species_and_colours %>% unique())

# merge with existing current_species_names_and_colours_df or the data you read from Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/MicrobialSpeciesNameAndColour.tsv
# exclude 'unclassified' and 'Others' and then merge again
current_species_names_and_colours_minus_unclassifiedandothers_df <- current_species_name_and_colours %>%
  bind_rows(new_species_and_colours_df) %>%
  arrange(species_name) %>% 
  filter(!species_name %in% c('unclassified', 'unassigned', 'Others'))

current_species_names_and_colours_only_unclassifiedandothers_df <- current_species_name_and_colours %>% 
  bind_rows(new_species_and_colours_df) %>%
  arrange(species_name) %>% 
  filter(species_name %in% c('unclassified', 'unassigned', 'Others'))

current_species_names_and_colours_update1_df <- current_species_names_and_colours_minus_unclassifiedandothers_df %>%
  bind_rows(current_species_names_and_colours_only_unclassifiedandothers_df) %>% distinct()


# filter only species of interest for the plot and convert dataframe to vector
kraken_G230317_species_name_and_colour <- current_species_names_and_colours_update1_df %>% 
  filter(species_name %in% c(kraken_G230317_species)) %>% dplyr::pull(species_colour, species_name)


plot_bar_kraken <- function(kraken_combined_report_final_df, species_colors){
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
  
  # plot
  kraken_combined_report_final_plot <- kraken_combined_report_final_df %>% 
    ggplot(aes(x = factor(Sample, levels = sample_name_elements), y = percentage, fill = species, label = Read_category)) + 
    geom_col(position = position_stack(reverse = TRUE)) +
    xlab('Sample') + ylab("Proportion (Species)") + labs(fill = 'Species') + 
    facet_grid(. ~ factor(Site, levels = c("Neg_Ctrl", "Anus", "Nose")), scales = 'free_x', space='free_x') +
    theme_bw() +
    geom_text(position = position_stack(vjust = 0.5, reverse = TRUE), color = "white") +
    scale_fill_manual(breaks = kraken_combined_report_final_species, values = species_colors) + 
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
  
  return(kraken_combined_report_final_plot)
}


kraken_G230317_plot <- plot_bar_kraken(kraken_G230317_modified_df, kraken_G230317_species_name_and_colour)
kraken_G230317_plot


write.table(kraken_G230317_modified_df, 'Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/combined_kraken_report_per_sequence_date/G230317_combined_parsed_kraken_report.tsv', row.names = FALSE, sep = '\t')

jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230317_kraken_plot.jpeg", width = 10000, height = 5000, units = 'px', res = 600)
print(kraken_G230317_plot)
dev.off()


