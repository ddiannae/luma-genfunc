library(readr)
library(dplyr)
library(ComplexHeatmap)

sets <- c("GSE45827", "GSE86374", "tcga")

all_deg <- lapply(sets, function(set){
  deg <- read_tsv(paste0("data/", set, "_deg.tsv"))
  deg <- deg %>% filter(adj_p < 0.05)
  deg$set <- toupper(set)
  return(deg %>% dplyr::select(ensembl, logfc, set)) 
})

all_deg <- bind_rows(all_deg)
all_deg <- tidyr::spread(all_deg, key = "set", value = "logfc", fill = 0)  
all_deg <- all_deg %>% filter(tcga != 0)
dm <- as.matrix(all_deg %>% select(-ensembl))
rownames(dm) <- all_deg$ensembl

Heatmap(dm, name = "logFC", show_row_names = FALSE,
        column_names_rot = 0, column_names_gp = gpar(fontsize = 20),
        column_names_centered = TRUE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 15)))