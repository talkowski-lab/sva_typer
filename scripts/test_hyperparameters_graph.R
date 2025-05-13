library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggthemes)


total_df <- read_delim('test/hyper_param_test/total_results.txt')
metaparam_df <- read_delim('test/hyper_param_test/hyperparam_table.txt')


df <- total_df %>%
  filter(region %in% c('hexamer_region', 'VNTR_region')) %>%
  mutate(length = end - start) %>%
  select(ID, region, iteration, length) %>%
  pivot_wider(names_from = region, values_from = length) %>%
  left_join(metaparam_df, by = c('iteration'='run')) %>%
  pivot_longer(cols = colnames(metaparam_df)[-1], names_to='hyperparam', values_to='value')
  

df %>%
  ggplot(aes(x = value, y = VNTR_region)) +
  geom_jitter(color='grey', alpha=0.2) +
  geom_smooth(method='loess') +
  facet_grid(ID ~ hyperparam, scales='free_x') + 
  theme_clean() + 
  theme(
    strip.text.x = element_text(size=12, face='bold'),
  )
