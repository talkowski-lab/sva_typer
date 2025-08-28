library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggthemes)
library(glue)

total_df <- read_delim('../test/benchmarking_seqs_hmmered.txt')


iter_df <- total_df %>%
  filter(region %in% c('hexamer_region', 'VNTR_region')) %>%
  select(ID, region, start, end) %>%
  mutate(length = end - start) %>%
  mutate(region = gsub("_region", "", region)) %>%
  pivot_wider(names_from = region, values_from = c(start, end, length), names_glue="{region}_{.value}") #%>%


matchup_df <- read_delim("../test/benchmarking_seqs.txt")
matchup_df <- matchup_df %>%
  left_join(iter_df %>%
    select(child_hap=ID, child_hexamer_length=hexamer_length, child_VNTR_length=VNTR_length)
  ) %>%
  left_join(iter_df %>%
    select(parent_hap=ID, parent_hexamer_length=hexamer_length, parent_VNTR_length=VNTR_length)
  )

matchup_df <- matchup_df %>%
  mutate(SVA_ID = sapply(strsplit(child_hap, "_"), function(s) glue("{s[3]}_{s[4]}"))) %>%
  mutate(
    family = sapply(strsplit(child_hap, "_"), function(s) s[1]),
    child_hap = sapply(strsplit(child_hap, "_"), function(s) glue("{s[1]}_{s[2]}")),
    parent_hap = sapply(strsplit(parent_hap, "_"), function(s) glue("{s[1]}_{s[2]}")),
  ) %>%
  select(SVA_ID, family, child_hap, parent_hap, everything()) %>%
  mutate(hex_diff = child_hexamer_length - parent_hexamer_length) %>%
  mutate(VNTR_diff = child_VNTR_length - parent_VNTR_length) 

matchup_df %>%
  write_delim("~/test/all_matchup_HGSVC2.txt", delim="\t" )

xxx <- matchup_df %>%
  group_by(SVA_ID, family) %>%
  filter(abs(sum(hex_diff)) > 5 | abs(sum(VNTR_diff)) > 5)

xxx %>% 
  write_delim("~/test/filtered_matchup_HGSVC2.txt", delim="\t")
