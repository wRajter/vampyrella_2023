# LOAD LIBRARIES
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)

# VARIABLES
cell <- 'cellCombined'
marker <- 'Full18S'
sim <- 'sim99'

# READ & PREPARE DATA

# read OTU abundance table and create a long table with an OTU id and its abundance in a given sample
otu_data <- read_tsv(file = sprintf('otu_summary_table_%s_%s_%s.tsv', marker, cell, sim), 
                     na = 'NA') %>%
  select(ID, A3, Mock, NH1, NH4, Sim17, Sim22, 
         Th16, Th38, Th40, X17007) %>%
  rename(otu_id = ID) %>%
  pivot_longer(-otu_id, names_to = 'sample_id', values_to = 'abund')


# read taxonomic assignment table
taxonomy <- read_tsv(file = sprintf('../../raw_data/OTU_summary_tables/otu_summary_table_%s_%s_%s.tsv', marker, cell, sim),
                     na = 'NA') %>%
  select(ID, Kingdom, Domain, Phyllum, Class, Order, Family, Genus, Species) %>% 
  rename(otu_id = ID) %>% 
  rename_all(tolower)

# join OTU abundance table and taxonomy tables together into a long table and
# create a relative abundance column (abundance divided by sum of the abundances
# for a given sample)

# pick a taxon of interest if you want:
spec_taxonomy <- taxonomy %>% 
  filter(order == 'Vampyrellida')

# use taxonomy or spec taxonomy if you want a specific taxonomic group 
otu_rel_abund <-inner_join(otu_data, taxonomy, by='otu_id') %>% 
  group_by(sample_id) %>% 
  mutate(rel_abund = abund / sum(abund)) %>% 
  ungroup() %>% 
  pivot_longer(cols=c('kingdom', 'domain', 'phyllum', 'class', 'order',
                      'family', 'genus', 'species'),
               names_to = 'level',
               values_to = 'taxon') %>% 
  mutate(sample_id = factor(sample_id,
                            levels = c('Total', 'Mock', 'NH1', 'NH4', 'Sim17',
                                       'Sim22', 'Th16', 'Th38', 'Th40', 'A3', 
                                       'X17007')))

# adjust some taxon names in the relative abundance table
otu_rel_abund <- otu_rel_abund %>% 
  mutate(taxon = str_replace_all(taxon, c('Bacteria_X' = 'Bacteria',
                                          'Archaea_X' = 'Archaea',
                                          'Thalassomyxa_clade' = '*Thalassomyxa* clade')),
         taxon = str_replace_na(taxon, replacement = 'Unclassified')) 

# create a relative abundance table for a taxon level of interest (taxon_level variable)
taxon_level <- 'phyllum'


taxon_rel_abund <- otu_rel_abund %>% 
  filter(level == taxon_level) %>% 
  group_by(sample_id, taxon) %>% 
  summarise(rel_abund = 100*(sum(rel_abund)), .groups = 'drop') %>% 
  replace(is.na(.), 0)



taxon_rel_abund <- taxon_rel_abund %>% 
  mutate(sample_place = case_when(sample_id == 'Mock' ~ 'Mock<br>sample',
                                  sample_id == 'NH1' ~ 'Neuenhähnen',
                                  sample_id == 'NH4' ~ 'Neuenhähnen',
                                  sample_id == 'Sim17' ~ 'Simmelried',
                                  sample_id == 'Sim22' ~ 'Simmelried',
                                  sample_id == 'Th16' ~ 'Thielenbruch',
                                  sample_id == 'Th38' ~ 'Thielenbruch',
                                  sample_id == 'Th40' ~ 'Thielenbruch',
                                  sample_id == 'A3' ~ 'Deep sea',
                                  sample_id == 'X17007' ~ 'Deep sea'))

# pool the taxonomic groups under certain percentage abundance (perc_treshold) 
# to reduce number of taxa
perc_treshold <- 10

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
  group_by(sample_id, sample_place, taxon) %>% 
  summarise(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups = 'drop') %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE)) %>% 
  ggplot(aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(~sample_place, scale = 'free_x', space = 'free', switch = 'y') +
  labs(x=NULL,
       y='Mean Relative Abundance (%)') +
  theme_classic() +
  theme(legend.title = element_blank(),
        # legend.key.size = unit(12, 'pt'),
        legend.text = ggtext::element_markdown(),
        strip.background = element_blank(),
        strip.text.x = ggtext::element_markdown()) +
  scale_fill_brewer(palette = "Paired")

ggsave(sprintf('stacked_bar_%s_%s_%s_%s.png', marker, cell, sim, taxon_level), 
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







