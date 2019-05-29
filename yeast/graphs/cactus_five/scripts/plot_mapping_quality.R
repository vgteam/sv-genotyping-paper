library(tidyverse)
library(ggrepel)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

construct <- read_tsv(args[1], col_names = c("mapq_gt0", "mapq_ge10", "mapq_ge20", "mapq_ge30", "mapq_ge40", "mapq_ge50", "mapq_ge60", "id_ge100", "id_ge90", "id_ge50", "all", "graph", "sample"))
cactus <- read_tsv(args[2], col_names = c("mapq_gt0", "mapq_ge10", "mapq_ge20", "mapq_ge30", "mapq_ge40", "mapq_ge50", "mapq_ge60", "id_ge100", "id_ge90", "id_ge50", "all", "graph", "sample"))

construct_new <- construct %>%
  gather("mapq_gt0", "mapq_ge10", "mapq_ge20", "mapq_ge30", "mapq_ge40", "mapq_ge50", "mapq_ge60", "id_ge100", "id_ge90", "id_ge50", key="filter", value="construct_number") %>%
  mutate(construct_fraction = construct_number / all) %>%
  select(sample, filter, construct_fraction)

cactus_new <- cactus %>%
  gather("mapq_gt0", "mapq_ge10", "mapq_ge20", "mapq_ge30", "mapq_ge40", "mapq_ge50", "mapq_ge60", "id_ge100", "id_ge90", "id_ge50", key="filter", value="cactus_number") %>%
  mutate(cactus_fraction = cactus_number / all) %>%
  select(sample, filter, cactus_fraction)

construct_new %>%
  inner_join(cactus_new, by=c("sample", "filter")) %>%
  filter(filter %in% c("mapq_gt0", "mapq_ge10", "mapq_ge20", "mapq_ge30", "mapq_ge40", "mapq_ge50", "mapq_ge60")) %>%
  mutate(filter = factor(filter, levels = c("mapq_gt0", "mapq_ge10", "mapq_ge20", "mapq_ge30", "mapq_ge40", "mapq_ge50", "mapq_ge60"), labels = c("0", "10", "20", "30", "40", "50", "60"))) %>%
  mutate(clade=ifelse(sample %in% c("UWOPS91-917.1","UFRJ50816","YPS138","N44","CBS432"), 'paradoxus', 'cerevisiae')) %>%
  mutate(ingraph=ifelse(sample %in% c("UFRJ50816", "YPS128", "CBS432", "SK1", "S288c"), 'included', 'excluded')) %>%
  ggplot(aes(construct_fraction, cactus_fraction, color=sample, alpha=ingraph, pch=clade)) +
  geom_point(aes(size=filter)) +
  geom_line() +
  # geom_label_repel(aes(label = sample),
  #                  box.padding   = 0.15, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50') +
  labs(color="Strain", size="Mapping quality threshold", x="Mapped read fraction on construct graph", y="Mapped read fraction on cactus graph", alpha="Included in graph", pch="Clade") +
  coord_cartesian(xlim=c(0.6,1), ylim=c(0.6,1)) +
  geom_abline(intercept=0) +
  theme_bw()

ggsave(args[3], device = "png", width = 8, height=7)
