knitr::opts_chunk$set(echo = TRUE)

library(tidyverse)
library(magrittr)
library(ggthemes)
library(scales)
library(kableExtra)
library(dplyr)


hDNA_microbialDNA_G230414 <- readr::read_tsv('Q:/IUK-A-MIGE/PROJECTS/TAPIR/files/hDNA_microbialDNA_proportion/mapped_stats_hDNA_G230414.txt')
head(hDNA_microbialDNA_G230414) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")


hDNA_microbialDNA_G230414_long_df <- hDNA_microbialDNA_G230414 %>% select(-c(Total, mapped_reads, unmapped_reads)) %>% 
  pivot_longer(-Filename, names_to = 'Read_proportion', values_to = 'Proportion') %>% 
  bind_cols(hDNA_microbialDNA_G230414 %>% select(-c(Total, proportion_mapped, proportion_unmapped)) %>% 
              pivot_longer(-Filename, names_to = 'Read', values_to = 'Counts') %>% select(-Filename))

head(hDNA_microbialDNA_G230414_long_df) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")


hDNA_microbialDNA_G230414_df <- hDNA_microbialDNA_G230414_long_df %>% 
  mutate(Site = case_when(str_detect(Filename, '[Ss]wab|[Ww]ater') ~ 'Neg_Ctrl',
                          str_detect(Filename, 'A') ~ 'Anus',
                          str_detect(Filename, 'B|N') ~ 'Nose')) %>% 
  mutate(Sample = str_extract(Filename, "^[a-zA-Z0-9]+")) %>% 
  mutate(Read_category = case_when(Counts < 1000 ~ "",
                                   Counts >= 1000 & Counts <= 10000 ~ "*",
                                   Counts > 10000 ~ "**")) 

head(hDNA_microbialDNA_G230414_df) %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")



human_contamination_plot2 <- function(bam_metadata_merged_df, ordered_sample_names){
  library(tidyverse)
  library(ggthemes)
  
  data_df_plot <- bam_metadata_merged_df %>% 
    mutate(Read = if_else(Read == 'mapped_reads', 'human DNA', 'microbial DNA')) %>% 
    ggplot(aes(x = factor(Sample, levels = ordered_sample_names), y = Proportion, fill = Read, label = Read_category)) + 
    geom_col(position = position_stack(reverse = TRUE)) +
    xlab('Sample') + ylab("Proportion (Reads)") +
    facet_grid(. ~ factor(Site, levels = c("Neg_Ctrl", "Anus", "Nose")), scales = 'free_x', space='free_x') +
    theme_bw() +
    geom_text(position = position_stack(vjust = 0.1, reverse = TRUE), color = "white", size = 10) +
    scale_fill_manual(values = c('Black', 'Grey')) +
    theme(
      axis.title.y = element_text(size = 20, face = 'bold'),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 30, face = 'bold'),
      legend.position="right",
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
      # plot.margin=unit(c(0,1,0,3.1),"cm"),
      axis.text.x = element_text(size = 20, face = 'bold', angle = 45, hjust = 1))
  
  return(data_df_plot)
}




sample_name_elements_G230414 <- hDNA_microbialDNA_G230414_df %>% 
  filter(!Sample %in% c('Swab', 'swab', 'Water', 'water')) %>% 
  distinct(Sample) %>% pull() %>% unique() %>% sort()

sample_name_elements_G230414 <- c('Swab', 'Water') %>% append(sample_name_elements_G230414)

sample_name_elements_G230414



hDNA_microbialDNA_G230414_plot <- human_contamination_plot2(hDNA_microbialDNA_G230414_df, sample_name_elements_G230414)

hDNA_microbialDNA_G230414_plot
  
jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230414_hDNA_microbialDNA_plot.jpeg", width = 10000, height = 5000, units = 'px', res = 600)
print(hDNA_microbialDNA_G230414_plot)
dev.off()