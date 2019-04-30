library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

data <- read_tsv(args[1], col_names=c("graph", "strain", "identity"))
samples <- read_tsv("strains.tsv", col_names = c("strain", "clade"))

data %>%
  spread(key = graph, value = identity) %>%
  inner_join(samples) %>%
  ggplot(aes(x=construct, cactus, color=strain, pch=clade)) +
  geom_point() +
  coord_cartesian(xlim=c(0.5,1), ylim=c(0.5,1)) +
  geom_abline(intercept=0) +
  theme_bw() +
  labs(x="percent identity on construct graph", y="percent identity on cactus graph")
  #geom_vline(aes(xintercept=construct, color=strain)) +
  #geom_hline(aes(yintercept=cactus, color=strain))

ggsave(args[2], device = "png", width = 6, height=5)
