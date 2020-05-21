library(readr)
library(dplyr)

sets <- c("GSE45827", "GSE86374", "tcga")

l <- lapply(sets, function(set){
  deg <- read_tsv(paste0("data/", set, "_deg.tsv"))
  write_tsv(deg %>% select(ensembl), paste0("data/", set, "_universe.txt"), 
            col_names = F)
  deg <- deg %>% filter(adj_p < 0.05)
  up <- deg %>% filter(logfc > 2) %>% select(ensembl)
  down <- deg %>% filter(logfc < -2) %>% select(ensembl)
  write_tsv(up, paste0("data/", set, "_up.txt"), col_names = F)
  write_tsv(down, paste0("data/", set, "_down.txt"), col_names = F)
})