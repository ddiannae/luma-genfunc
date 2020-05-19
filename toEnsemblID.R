library(org.Hs.eg.db)
library(readr)
library(dplyr)

sets <- c("GSE45827", "GSE86374")

l <- lapply(sets, function(set){
  deg <- read_tsv(paste0("data/", set, "_luma_vs_healthy.txt"), 
                  col_names = c("id", "adj_p", "p", "t", "B", "logfc", "entrez_id"), 
                  col_types = "cddddd--c",
                  skip = 1)
  deg <- deg %>% filter(!is.na(entrez_id))
  deg <- deg %>% mutate(entrez_id = strsplit(entrez_id , "///"))
  deg <- deg %>% tidyr::unnest(cols = c(entrez_id))
  deg$ensembl <- mapIds(org.Hs.eg.db,
                        keys = deg$entrez_id,
                        column="ENSEMBL",
                        keytype="ENTREZID",
                        multiVals="list")
  deg <- deg %>% tidyr::unnest(cols = c(ensembl))
  deg <- deg %>% filter(!is.na(ensembl))
  deg <- deg %>% filter(!duplicated(ensembl))   
  deg <- deg %>% dplyr::select(ensembl, everything()) %>% dplyr::select(-id)
  write_tsv(deg, paste0("data/", set, "_deg.tsv"))
})

tcga <- read_tsv("data/luma-deg-ebayes.tsv")
colnames(tcga) <- c("ensembl", "logfc", "p", "adj_p", "B")
write_tsv(tcga, "data/tcga_deg.tsv")
