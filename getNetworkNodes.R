library(readr)
library(dplyr)

sets <- c("GSE45827", "GSE86374", "tcga")

all_deg <- lapply(sets, function(set){
  deg <- read_tsv(paste0("data/", set, "_deg.tsv"))
  deg <- deg %>% filter(adj_p < 0.05)
  deg$set <- toupper(set)
  return(deg %>% dplyr::select(ensembl, logfc, set)) 
})

all_deg <- bind_rows(all_deg)
all_deg %>% filter(abs(logfc) > 2) %>% group_by(set) %>% tally()

all_deg <- tidyr::spread(all_deg, key = "set", value = "logfc", fill = 0)  
colnames(all_deg)[1] <- "ensemblID"
vertices <- read_tsv(file = paste0("data/luma-20127-vertices.tsv")) 
vertices <- vertices %>% select(ensemblID)
vertices <- vertices %>% left_join(all_deg, by = "ensemblID")
write_tsv(vertices, "data/luma-20127-vertices-deg.tsv")
