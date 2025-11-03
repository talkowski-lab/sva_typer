library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggthemes)
library(ggh4x)
library(ggforce)
library(glue)

here::i_am(".gitignore")
setwd(here::here())


total_df <- read_delim('test/benchmarking_seqs_hmmered.txt')


iter_df <- total_df %>%
  filter(region %in% c('hexamer_region', 'VNTR_region')) %>%
  select(ID, region, start, end) %>%
  mutate(length = end - start) %>%
  mutate(region = gsub("_region", "", region)) %>%
  pivot_wider(names_from = region, values_from = c(start, end, length), names_glue="{region}_{.value}") #%>%


matchup_df <- read_delim("test/benchmarking_seqs.txt")
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
  write_delim("~/test/all_matchup_HGSVC3.txt", delim="\t" )


stats_matchup <- matchup_df %>%
  group_by(family) %>%
  summarize(
    cor_hex = cor(child_hexamer_length, parent_hexamer_length, use = "pairwise.complete.obs"),
    cor_VNTR = cor(child_VNTR_length, parent_VNTR_length, use = "pairwise.complete.obs"),
    diff_hex = mean(abs(hex_diff), na.rm=TRUE),
    diff_VNTR = mean(abs(VNTR_diff), na.rm=TRUE)
  ) %>%
  pivot_longer(-family) %>%
  separate_wider_delim(name, "_", names=c("stat", "tr"))

stats_matchup %>%
  filter(stat == "cor") %>%
  ggplot(aes(x=family, y=value * value, fill=family)) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
  scale_fill_tableau() +
  facet_row(tr ~ .) +
  labs(y="R^2 between parent and child") +
  theme_clean() + 
  theme(
    legend.position = "none"
  )

ggsave("test/figures/trio_repeat_correlations.pdf", width=5, height=4)

stats_matchup %>%
  filter(stat == "diff") %>%
  ggplot(aes(x=family, y=value, fill=family)) +
  geom_bar(stat="identity") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
  scale_fill_tableau() +
  facet_row(tr ~ .) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(y="Mean diff between parent and child") +
  theme_clean() + 
  theme(
    legend.position = "none"
  )
ggsave("test/figures/trio_repeat_diff.pdf", width=5, height=4)

