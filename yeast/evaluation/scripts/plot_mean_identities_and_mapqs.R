library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

five <- read_tsv(args[1], col_names=c("strain", "cactus_id", "vcf_id", "linear_id", "cactus_mapq", "vcf_mapq", "linear_mapq")) %>% 
    mutate(clade=ifelse(strain %in% c("UWOPS919171","UFRJ50816","YPS138","N44","CBS432"), 'paradoxus', 'cerevisiae')) %>%
    mutate(ingraph=ifelse(strain %in% c("UFRJ50816", "YPS128", "CBS432", "SK1", "S288c"), 'included', 'excluded')) %>%
    mutate(graph = "Five strains")
all <- read_tsv(args[2], col_names=c("strain", "cactus_id", "vcf_id", "linear_id", "cactus_mapq", "vcf_mapq", "linear_mapq")) %>%
    mutate(clade=ifelse(strain %in% c("UWOPS919171","UFRJ50816","YPS138","N44","CBS432"), 'paradoxus', 'cerevisiae')) %>%
    mutate(ingraph='included') %>%
    mutate(graph = "All strains")

total <- rbind(five, all) %>%
    mutate(diff_id = cactus_id - vcf_id) %>%
    mutate(diff_mapq = cactus_mapq - vcf_mapq)

total %>%
  ggplot(aes(x=strain, y=diff_id, fill=graph, alpha=ingraph)) +
  geom_col(position=position_dodge()) +
  scale_alpha_discrete(range=c(.4, 1)) +
  facet_grid(~clade, scales='free') +
  labs(fill="Graph type", alpha="During graph\nconstruction",
       x="Yeast strain",
       y="Average delta in mapping identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text=element_text(size=14))

ggsave(args[3], device = "png", width = 6, height=5)

total %>%
  ggplot(aes(x=strain, y=diff_mapq, fill=graph, alpha=ingraph)) +
  geom_col(position=position_dodge()) +
  scale_alpha_discrete(range=c(.4, 1)) +
  facet_grid(~clade, scales='free') +
  labs(fill="Graph type", alpha="During graph\nconstruction",
       x="Yeast strain",
       y="Average delta in mapping quality") +
  theme_bw() +
  guides(fill=FALSE, alpha=FALSE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text=element_text(size=14))

ggsave(args[4], device = "png", width = 6, height=5)
