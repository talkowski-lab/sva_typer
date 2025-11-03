library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggthemes)


here::i_am(".gitignore")
setwd(here::here())


total_df <- read_delim('test/hyper_param_test/total_results.txt')
metaparam_df <- read_delim('test/hyper_param_test/hyperparam_table.txt')


iter_df <- total_df %>%
  filter(region %in% c('hexamer_region', 'VNTR_region')) %>%
  select(ID, region, iteration, start, end) %>%
  mutate(length = end - start) %>%
  mutate(region = gsub("_region", "", region)) %>%
  pivot_wider(names_from = region, values_from = c(start, end, length), names_glue="{region}_{.value}") #%>%

df <- iter_df %>%
  left_join(metaparam_df, by = c('iteration'='run')) %>%
  pivot_longer(cols = colnames(metaparam_df)[-1], names_to='hyperparam', values_to='value') %>%
  filter(!(ID %in% c("CH10_SVA_F_2", "SVA"))) %>%
  mutate(ID = factor(ID, levels=c("SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F", "CH10_SVA_F_1")))
  


# Hexamer Graph
df %>%
  ggplot(aes(x = value, y = hexamer_length)) +
  geom_jitter(color='grey', alpha=0.2) +
  geom_smooth(method='loess') +
  facet_grid(ID ~ hyperparam, scales='free_x') + 
  coord_cartesian(ylim = c(0, 250)) + 
  theme_clean() + 
  theme(
    strip.text.x = element_text(size=12, face='bold'),
    strip.clip = "off"
  )
ggsave('test/hyper_param_test/hexamer_graph.png', width=13, height=6, dpi=300)
# VNTR Graph
df %>%
  ggplot(aes(x = value, y = VNTR_length)) +
  geom_jitter(color='grey', alpha=0.2) +
  geom_smooth(method='loess') +
  facet_grid(ID ~ hyperparam, scales='free_x') + 
  theme_clean() + 
  theme(
    strip.text.x = element_text(size=12, face='bold'),
    strip.clip = "off"
  )
ggsave('test/hyper_param_test/VNTR_graph.png', width=13, height=6, dpi=300)
