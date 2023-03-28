# dataframe picked from bracken_plot_230119.Rmd, bracken_plot_230126.Rmd, and bracken_plot_230202.Rmd

sample_list_G230119 <- bracken_G230119_modified_df %>% filter(!Sample %in% c('Swab', 'Water')) %>% pull(Sample) %>% unique()
sample_list_G230126 <- bracken_G230126_modified_df %>% filter(!Sample %in% c('Swab', 'Water')) %>% pull(Sample) %>% unique()
sample_list_G230202 <- bracken_G230202_modified_df %>% filter(!Sample %in% c('Swab', 'Water')) %>% pull(Sample) %>% unique()

common_samplenames_G230119toG230202 <- Reduce(intersect, list(sample_list_G230119, sample_list_G230126, sample_list_G230202))

combine_G230119toG230208_selection_df <- bracken_G230119_modified_df %>% bind_rows(bracken_G230126_modified_df) %>% 
  bind_rows(bracken_G230202_modified_df) %>%  bind_rows(bracken_G230208_modified_df) %>% 
  filter(Sample %in% c(common_samplenames_G230119toG230202)) %>% 
  mutate(sequencedate = case_when(str_detect(Filename, '-2-') ~ 'week2',
                                  str_detect(Filename, '-3-') ~ 'week3',
                                  str_detect(Filename, '-4') ~ 'week4',
                                  str_detect(Filename, '-5') ~ 'week5'))

combine_G230119toG230202_selection_species <- combine_G230119toG230202_selection_df %>% 
  filter(!name %in%  c('Others', 'unassigned')) %>% distinct(name) %>% 
  pull() %>% unique() %>% sort() %>% c('Others', 'unassigned')

# setdiff(combine_G230119toG230202_selection_species, current_species_name_and_colours$species_name)
combine_G230119toG230202_selection_species_name_and_colour <- current_species_names_and_colours_update1_df %>% 
  filter(species_name %in% c(combine_G230119toG230202_selection_species)) %>% dplyr::pull(species_colour, species_name)

# plot
bracken_combine_G230119toG230208_selection_plot <- plot_bar_bracken_select(combine_G230119toG230208_selection_df, combine_G230119toG230202_selection_species_name_and_colour)

ggsave("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230119toG230202_common_sample_selection.jpeg", bracken_combine_G230119toG230202_selection_plot, device = "jpeg", width = 6000, height = 8000, units = 'px', dpi = 600)

jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230119toG230208_common_sample_selection_bracken.jpeg", width = 10000, height = 8000, units = 'px', res = 600)
print(bracken_combine_G230119toG230208_selection_plot)
dev.off()


################### kraken together #############################

combine_G230119toG230208_selection_kraken_df <- kraken_G230119_modified_df %>% bind_rows(kraken_G230126_modified_df) %>% 
  bind_rows(kraken_G230202_modified_df) %>% bind_rows(kraken_G230208_modified_df) %>% filter(Sample %in% c(common_samplenames_G230119toG230202)) %>% 
  mutate(sequencedate = case_when(str_detect(samplename, '-2-') ~ 'week2',
                                  str_detect(samplename, '-3-') ~ 'week3',
                                  str_detect(samplename, '-4') ~ 'week4',
                                  str_detect(samplename, '-5') ~ 'week5'))

combine_G230119toG230208_selection_kraken_species <- combine_G230119toG230208_selection_kraken_df %>% 
  filter(!species %in%  c('Others', 'unassigned')) %>% distinct(species) %>% 
  pull() %>% unique() %>% sort() %>% c('Others', 'unassigned')

combine_G230119toG230208_selection_kraken_species_name_and_colour <- current_species_names_and_colours_update1_df %>% 
  filter(species_name %in% c(combine_G230119toG230208_selection_kraken_species)) %>% dplyr::pull(species_colour, species_name)


# plot
kraken_combine_G230119toG230208_selection_plot <- plot_bar_kraken_select(combine_G230119toG230208_selection_kraken_df, combine_G230119toG230208_selection_kraken_species_name_and_colour)

jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230119toG230208_common_sample_selection_kraken.jpeg", width = 10000, height = 8000, units = 'px', res = 600)
print(kraken_combine_G230119toG230208_selection_plot)
dev.off()

################### kraken week5 and week5extra together #############################
sample_list_G230206 <- kraken_G230206_modified_df %>% filter(!Sample %in% c('Swab', 'Water')) %>% pull(Sample) %>% unique()
sample_list_G230208 <- kraken_G230208_modified_df %>% filter(!Sample %in% c('Swab', 'Water')) %>% pull(Sample) %>% unique()

common_samplenames_G230206toG230208 <- intersect(sample_list_G230206, sample_list_G230208)

combine_G230206toG230208_selection_kraken_df <- kraken_G230206_modified_df %>% 
  mutate(sequencedate = 'week5Extra') %>% bind_rows(kraken_G230208_modified_df %>% mutate(sequencedate = 'week5')) %>% 
  filter(Sample %in% c(common_samplenames_G230206toG230208)) %>% filter(Site != 'Nose') 


combine_G230206toG230208_selection_kraken_species <- combine_G230206toG230208_selection_kraken_df %>% 
  filter(!species %in%  c('Others', 'unassigned')) %>% distinct(species) %>% 
  pull() %>% unique() %>% sort() %>% c('Others', 'unassigned')

combine_G230206toG230208_selection_kraken_species_name_and_colour <- current_species_names_and_colours_update1_df %>% 
  filter(species_name %in% c(combine_G230206toG230208_selection_kraken_species)) %>% dplyr::pull(species_colour, species_name)

# plot
kraken_combine_G230206toG230208_selection_plot <- plot_bar_kraken_select(combine_G230206toG230208_selection_kraken_df, combine_G230206toG230208_selection_kraken_species_name_and_colour)

jpeg("Q:/IUK-A-MIGE/PROJECTS/TAPIR/figures/G230206toG230208_common_sample_selection_kraken.jpeg", width = 12000, height = 7000, units = 'px', res = 600)
print(kraken_combine_G230206toG230208_selection_plot)
dev.off()
