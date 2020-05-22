library(stringr)
library(readr)
library(dplyr)
library(VennDiagram)
library(scales)

## SET 1
upr_s1 <- read_tsv("webgestalt/GSE45827_up/enrichment_results_wg_result1589950101.txt")
upr_s1 <- upr_s1 %>% select(geneSet, description, enrichmentRatio, FDR)
downr_s1 <- read_tsv("webgestalt/GSE45827_down/enrichment_results_wg_result1589949772.txt")
downr_s1 <- downr_s1 %>% select(geneSet, description, enrichmentRatio, FDR)

# Chart
v <- venn.diagram(
  main = "GSE45827",
  x = list(upr_s1$geneSet, downr_s1$geneSet),
  category.names = c("GO genes up" , "GO genes down"),
  filename = "figures/GSE45827_webgestalt.png",
  output = TRUE,
  main.fontfamily = "sans",
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  col = c("#a85763", '#57A89C'),
  fill = c(alpha("#a85763",0.8), alpha('#57A89C',0.8)),
  # Numbers
  cex = 0.8,
  fontfamily = "sans",
  # Set names
  cat.fontfamily = "sans",
  cat.pos = c(0, 0),
  cat.dist = c(0.03, 0.03),
  cat.col = c("#a85763", '#57A89C'),
  cat.cex = 0.6
)

### Set TCGA
upr_s3 <- read_tsv("webgestalt/TCGA_up/enrichment_results_wg_result1589951738.txt")
upr_s3 <- upr_s3 %>% select(geneSet, description, enrichmentRatio, FDR)
downr_s3 <- read_tsv("webgestalt/TCGA_down/enrichment_results_wg_result1589951617.txt")
downr_s3 <- downr_s3 %>% select(geneSet, description, enrichmentRatio, FDR)

# Chart
v <- venn.diagram(
  main = "TCGA",
  x = list(upr_s3$geneSet, downr_s3$geneSet),
  category.names = c("GO genes up" , "GO genes down"),
  filename = "figures/TCGA_webgestalt.png",
  output = TRUE,
  main.fontfamily = "sans",
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  margin = 0.02,
  # Circles
  lwd = 2,
  col = c("#E78618", '#E0E718'),
  fill = c(alpha("#E78618", 0.8), alpha('#E0E718',0.3)),
  # Numbers
  cex = 0.6,
  fontfamily = "sans",
  ext.text = FALSE,
  # Set names
  cat.fontfamily = "sans",
  cat.pos = c(30, 0),
  cat.dist = c(0.03, -0.06),
  cat.col = c("#E78618", '#E0E718'),
  cat.cex = 0.6
)

# Chart
v <- venn.diagram(
  x = list(upr_s1$geneSet, downr_s1$geneSet, upr_s3$geneSet, downr_s3$geneSet),
  category.names = c("GSE45827 up" , "GSE45827 down", 
                     "TCGA up" , "TCGA down"),
  filename = "figures/all_venn.png",
  output = TRUE,
  main.fontfamily = "sans",
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  margin = 0.01,
  #  # Circles
  lwd = 2,
  lty = 'blank',
  col = c("#a85763", '#57A89C', "#E78618", '#E0E718'),
  fill = c(alpha("#a85763",0.8), alpha('#57A89C',0.8),
           alpha("#E78618", 0.8), alpha('#E0E718',0.3)),
  #  # Numbers
  cex = 0.4,
  fontfamily = "sans",
  #  
  #  # Set names
  cat.col = c("#a85763", '#57A89C', "#E78618", '#E0E718'),
  cat.cex = 0.5,
  # Set names
  cat.fontfamily = "sans",
  cat.pos = c(165, 200, 0, 0),
  cat.dist = c(0.45, 0.45, 0.1,0.1)
)

gsea_results <- read_tsv("gsea/TCGA/gsea_report_for_1_1590041610752.xls")
gsea_results <- gsea_results %>% select(NAME, `FDR q-val`) %>% 
  rename("name" = "NAME",  "fdr" = "FDR q-val") %>% filter(fdr < 0.25)
gsea_results <- gsea_results %>% 
  mutate(name = tolower(str_remove(str_replace_all(name, "_", " "), "GO ")))
write_tsv(gsea_results, "gsea/tcga_simple_results.tsv")
                        