library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggthemes)


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
  pivot_longer(cols = colnames(metaparam_df)[-1], names_to='hyperparam', values_to='value') 
  


# Hexamer Graph
df %>%
  ggplot(aes(x = value, y = hexamer_region)) +
  geom_jitter(color='grey', alpha=0.2) +
  geom_smooth(method='loess') +
  facet_grid(ID ~ hyperparam, scales='free_x') + 
  coord_cartesian(ylim = c(0, 250)) + 
  theme_clean() + 
  theme(
    strip.text.x = element_text(size=12, face='bold'),
  )
ggsave('test/hyper_param_test/hexamer_graph.png', width=12, height=8, dpi=300)
# VNTR Graph
df %>%
  ggplot(aes(x = value, y = VNTR_region)) +
  geom_jitter(color='grey', alpha=0.2) +
  geom_smooth(method='loess') +
  facet_grid(ID ~ hyperparam, scales='free_x') + 
  theme_clean() + 
  theme(
    strip.text.x = element_text(size=12, face='bold'),
  )
ggsave('test/hyper_param_test/VNTR_graph.png', width=12, height=8, dpi=300)
