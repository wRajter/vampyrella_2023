# LOAD LIBRARIES
library(tidyverse)
library(readxl)
library(ggtext)
library(RColorBrewer)


# READ & PREPARE DATA
asv_data <- read_tsv(file = 'asv_summary_table_Full18S_cell2.tsv', na = 'NA') %>%
  select(ID, Total, A3, Mock, NH1, NH4, Sim17, Sim22, 
         Th16, Th38, Th40, X17007) %>%
  rename(asv_id = ID) %>%
  pivot_longer(-asv_id, names_to = 'sample_id', values_to = 'abundance')

taxonomy <- read_tsv(file = 'asv_summary_table_Full18S_cell2.tsv', na = 'NA') %>%
  select(ID, Kingdom, Domain, Phyllum, Class, Order, Family, Genus, Species) %>% 
  rename(asv_id = ID) %>% 
  rename_all(tolower)


asv_rel_abund <-inner_join(asv_data, taxonomy, by='asv_id') %>% 
  group_by(sample_id) %>% 
  mutate(rel_abun = abundance / sum(abundance)) %>% 
  ungroup() %>% 
  select(-abundance) %>% 
  pivot_longer(cols=c('kingdom', 'domain', 'phyllum', 'class', 'order',
                      'family', 'genus', 'species'),
               names_to = 'level',
               values_to = 'taxon') %>% 
  mutate(sample_id = factor(sample_id,
                            levels = c('Total', 'Mock', 'NH1', 'NH4', 'Sim17',
                                       'Sim22', 'Th16', 'Th38', 'Th40', 'A3', 
                                       'X17007')))

taxon_rel_abund <- asv_rel_abund %>% 
  filter(level=='domain') %>% 
  group_by(sample_id, taxon) %>% 
  summarise(rel_abun = sum(rel_abun), .groups = 'drop') %>% 
  group_by(sample_id, taxon) %>% 
  summarise(mean_rel_abund = 100*mean(rel_abun), .groups = 'drop') %>% 
  mutate(taxon = str_replace_all(taxon,
                             c('Bacteria_X' = 'Bacteria', 'Archaea_X' = 'Archaea')),
         taxon = str_replace_na(taxon, replacement = 'Unclassified')) 
  
taxon_pool <- taxon_rel_abund %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean_rel_abund) < 6, .groups = 'drop')



# STACKED BARCHART

taxon_pool <- taxon_rel_abund %>% 
  group_by(taxon) %>% 
  summarise(pool = max(mean_rel_abund) < 6, .groups = 'drop')

inner_join(taxon_rel_abund, taxon_pool, by = 'taxon')  %>% 
  mutate(taxon = if_else(pool, 'Other (< 6%)', taxon)) %>% 
  group_by(sample_id, taxon) %>% 
  summarise(mean_rel_abund = sum(mean_rel_abund), .groups = 'drop') %>% 
  ggplot(aes(x=sample_id, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x=NULL,
       y='Mean Relative Abundance (%)') +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.key.size = unit(12, 'pt')) +
  scale_fill_brewer(palette = "Dark2")

ggsave('stacked_bar_domain.tiff', width=7, height=5)


# HEATMAP

inner_join(taxon_rel_abund, taxon_pool, by = 'taxon')  %>% 
  mutate(taxon = if_else(pool, 'Other (< 6%)', taxon)) %>% 
  group_by(sample_id, taxon) %>% 
  summarise(mean_rel_abund = sum(mean_rel_abund), .groups = 'drop') %>%
  ggplot(aes(x=sample_id, fill=mean_rel_abund, y= reorder(taxon, desc(taxon)))) +
  geom_tile() +
  geom_text(aes(label=format(round(mean_rel_abund, 1)), nsmall=1), size=2.6) +
  scale_fill_gradient("Mean<br>Relative<br>Abundance<br>(%)",
                      low = '#FFFFFF', high = '#FF0000',
                      expand = c(0, 0)) +
  labs(x=NULL,
       y=NULL) +
  theme_classic() +
  theme(legend.title = element_markdown(size = 9, 'pt'),
        legend.title.align = 0.5,
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed(ratio = 0.7)
  
ggsave('heatmap_domain.tiff', width=7, height=5)



# grouped box plot and range plot 

taxon_pool <- taxon_rel_abund %>% 
  group_by(taxon) %>% 
  summarize(mean=mean(mean_rel_abund), .groups="drop") %>% 
  group_by(taxon) %>%
  summarize(pool = max(mean) < 3,
            mean = mean(mean),
            .groups="drop")

inner_join(taxon_rel_abund, taxon_pool, by = 'taxon')  %>% 
  mutate(taxon = if_else(pool, 'Other (< 6%)', taxon)) %>% 
  group_by(sample_id, taxon) %>% 
  summarise(mean_rel_abund = sum(mean_rel_abund), .groups = 'drop') %>%
ggplot(aes(x=taxon, y=mean_rel_abund, fill=sample_id)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4,
                                              dodge.width =0.8),
              pch=21, stroke=0, size=1.8) +
  stat_summary(fun=median, geom = "crossbar",
               position = position_dodge(width=0.8),
               width=0.8,
               size=0.25,
               show.legend=FALSE) +
  scale_y_continuous() +
  labs(x=NULL,
       y="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        legend.position = c(0.8, 0.9),
        legend.background = element_rect(color="black", fill = NA),
        legend.margin = margin(t=-5, r=3, b=3)
  )


