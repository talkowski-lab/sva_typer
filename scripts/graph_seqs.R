library(ggplot2)
library(readr)
library(stringr)
library(ggthemes)

df <- read_delim('../test/seqs_hmmered.txt')

df %>%
  filter(!(region %in% c("hexamer_region", "VNTR_region"))) %>%
  mutate(region = if_else(grepl("skip", region, fixed=TRUE), "skip", region)) %>%
  group_by(ID) %>%
  mutate(total_length = max(end)) %>%
  mutate(mid = (start+end)/2, length = (end - start)) %>%
    ggplot(aes(y=reorder(ID, total_length))) +
    geom_tile(aes(x=mid, fill=region, width=length), height=0.95) +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_manual(
      values=c("hex"='red', "VNTR_1"="#92cbdf", "VNTR_2"="#2596be", "VNTR_3"="#134b5f", "skip"='grey')
    ) + 
    labs(x = 'Position', y = 'Haplotype', fill = 'Repeat type') +
    theme_clean() +
    theme(
      panel.grid.major.y = element_blank()
    )
  
