# LOAD LIBRARIES
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)
library(dplyr)
library(readr)
library(tidyr)

# VARIABLES

project <- 'Jamy_2022'
cell <- 'cell'
marker <- 'rDNA'
sim <- 'sim97'

# READ & PREPARE DATA

# read OTU abundance table and create a long table with an OTU id and its abundance in a given sample
otu_data <- read_tsv(file = sprintf('%s_OTU_summary_table_%s_%s_%s.tsv', project, marker, cell, sim), 
                     na = 'NA') %>%
  select(otu_id, sample, abundance) 



# read taxonomic assignment table
taxonomy <- read_tsv(file = sprintf('%s_otu_summary_table_%s_%s_%s.tsv', project, marker, cell, sim),
                     na = 'NA') %>%
  select(otu_id, kingdom, domain, phyllum, class, order, family, genus, species)

# join OTU abundance table and taxonomy tables together into a long table and
# create a relative abundance column (abundance divided by sum of the abundances
# for a given sample)

# pick a taxon of interest if you want:
spec_taxonomy <- taxonomy %>% 
  filter(order == 'Vampyrellida')

# use taxonomy or spec taxonomy if you want a specific taxonomic group 
otu_rel_abund <-inner_join(otu_data, spec_taxonomy, by='otu_id') %>%
  rename(abund = abundance) %>%
  rename(sample_id = sample) %>%
  group_by(sample_id) %>% 
  mutate(rel_abund = abund / sum(abund)) %>% 
  ungroup() %>% 
  pivot_longer(cols=c('kingdom', 'domain', 'phyllum', 'class', 'order',
                      'family', 'genus', 'species'),
               names_to = 'level',
               values_to = 'taxon')
  #%>%
  #mutate(sample_id = factor(sample_id,
  #                          levels = c('Total', 'Mock', 'NH1', 'NH4', 'Sim17',
  #                                     'Sim22', 'Th16', 'Th38', 'Th40', 'A3', 
  #                                     'X17007')))

# adjust some taxon names in the relative abundance table
otu_rel_abund <- otu_rel_abund %>% 
  mutate(taxon = str_replace_all(taxon, c('Bacteria_X' = 'Bacteria',
                                          'Archaea_X' = 'Archaea',
                                          'Thalassomyxa_clade' = '*Thalassomyxa* clade',
                                          'Vampyrellida_clade_B1' = 'clade B1', 
                                          'Vampyrellida_clade_B2' = 'clade B2',
                                          'Vampyrellida_clade_A' = 'clade A',
                                          'Vampyrellida_clade_PurpleVamp' = 'clade Purple-Vamp')),
         taxon = str_replace_na(taxon, replacement = 'Unclassified')) 

# create a relative abundance table for a taxon level of interest (taxon_level variable)
taxon_level <- 'species'


taxon_rel_abund <- otu_rel_abund %>% 
  filter(level == taxon_level) %>% 
  group_by(sample_id, taxon) %>% 
  summarise(rel_abund = 100*(sum(rel_abund)), .groups = 'drop') %>% 
  replace(is.na(.), 0)


# print sample names
unique(otu_data$sample)

# renaming the sample names in Suthaus 2022 data
taxon_rel_abund <- taxon_rel_abund %>% 
  mutate(sample_place = case_when(sample_id == 'Mock' ~ 'Mock',
                                  sample_id == 'NH1' ~ 'Neuen-<br>hähnen',
                                  sample_id == 'NH4' ~ 'Neuen-<br>hähnen',
                                  sample_id == 'Sim17' ~ 'Simmelried',
                                  sample_id == 'Sim22' ~ 'Simmelried',
                                  sample_id == 'Th16' ~ 'Thielenbruch',
                                  sample_id == 'Th38' ~ 'Thielenbruch',
                                  sample_id == 'Th40' ~ 'Thielenbruch',
                                  sample_id == 'A3' ~ 'Deep sea',
                                  sample_id == 'X17007' ~ 'Deep sea',
                                  sample_id == 'X579' ~ 'Deep sea', 
                                  sample_id == 'FS1' ~ 'FS samples', 
                                  sample_id == 'FS2' ~ 'FS samples', 
                                  sample_id == 'FS7' ~ 'FS samples', 
                                  sample_id == 'FS9' ~ 'FS samples', 
                                  sample_id == 'TS_3C_D' ~ 'TS samples', 
                                  sample_id == 'TS_3C' ~ 'TS samples', 
                                  sample_id == 'TS_4b2_D' ~ 'TS samples', 
                                  sample_id == 'TS_4b2' ~ 'TS samples', 
                                  sample_id == 'TS_9C' ~ 'TS samples', 
                                  sample_id == 'TS_HM22' ~ 'TS samples'))

# renaming the sample names in Jamy 2022 data
taxon_rel_abund <- taxon_rel_abund %>% 
  mutate(sample_id = case_when(sample_id == 'ERR6454469' ~ 'Jädraås SE (soil)',
                                  sample_id == 'ERR6454466' ~ 'Skogaryd SE (soil)',
                                  sample_id == 'ERR6454465' ~ 'Erken SE (sediment)',
                                  sample_id == 'ERR6454467' ~ 'Kallkäls SE (soil)',
                                  sample_id == 'ERR6454464' ~ 'Whapmagoostui<br>CA (permafrost)',
                                  sample_id == 'ERR6454468' ~ 'El Yunque PR (soil)',
                                  sample_id == 'ERR6454462' ~ 'Skogaryd SE (water)',
                                  sample_id == 'ERR6454463' ~ 'Svartberget SE (water)',
                                  sample_id == 'ERR6454461' ~ 'Erken SE (water)'))



# pool the taxonomic groups under certain percentage abundance (perc_treshold) 
# to reduce number of taxa
perc_treshold <- 15

taxon_pool <- taxon_rel_abund %>%
  group_by(sample_id, taxon) %>%
  summarise(mean = mean(rel_abund), .groups = 'drop') %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean) < perc_treshold,
            mean = mean(mean),
            .groups = 'drop')


# GROUPED STACKED BARCHART

inner_join(taxon_rel_abund, taxon_pool, by = 'taxon')  %>% 
  mutate(taxon = if_else(pool, sprintf('Other (< %s%%)', perc_treshold), taxon)) %>% 
  group_by(sample_id, taxon) %>% # sample_place 
  summarise(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups = 'drop') %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE)) %>% 
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0)) +
  #facet_grid(~sample_place, scale = 'free_x', space = 'free', switch = 'y') +
  labs(x=NULL,
       y='Mean Relative Abundance (%)') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.key.size = unit(6.5, 'pt'),
        axis.text.x = ggtext::element_markdown(angle=45, hjust=1, size=7), #element_text(size=8, angle=45, hjust=1),
        legend.text = ggtext::element_markdown(),
        strip.background = element_blank(),
        strip.text.x = ggtext::element_markdown(size=6.5)) +
  scale_fill_brewer(palette = "Paired")

ggsave(sprintf('%s_%s_%s_%s_%s_stacked_bar.png', project, marker, cell, sim, taxon_level), 
       device='png',
       width=7, 
       height=5,
       dpi=300)


# GROUPED UNSTACKED BARCHART
# how taxa vary across the samples
inner_join(taxon_rel_abund, taxon_pool, by = 'taxon')  %>% 
  mutate(taxon = if_else(pool, sprintf('Other (< %s%%)', perc_treshold), taxon)) %>% 
  group_by(sample_id, sample_place, taxon) %>% 
  summarise(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups = 'drop') %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = FALSE)) %>% 
  ggplot(aes(fill=sample_place, y=taxon, x=rel_abund)) +
  geom_col(position = position_dodge(), width = 0.7) +
  scale_x_continuous(expand = c(0, 0)) +
  # facet_grid(~sample_place, scale = 'free_x', space = 'free', switch = 'y') +
  labs(y=NULL,
       x='Mean Relative Abundance (%)') +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.key.size = unit(12, 'pt'),
        strip.background = element_blank(),
        legend.position = c(0.85, 0.8),
        panel.grid.major.x = element_line(color = 'lightgray', size = 0.2)) +
  scale_fill_brewer(palette = "Paired")

ggsave(sprintf('unstacked_bar_%s_%s_%s.tiff', marker, cell, taxon_level), width=5, height=7)







