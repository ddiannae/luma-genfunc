library(readr)
library(tidyr)
library(stringr)
library(dplyr)

sets <- c("GSE45827", "GSE86374", "tcga")

##### SET 1
set <- sets[1]
series <- read_tsv(paste0("data/", set, "_series_matrix.txt"), 
                skip = 31, n_max = 35, col_names = F)
series <- t(series)[, c(1, 2)]
colnames(series) <- c("title", "geo")
series <- series[-1, ]
series <- as_tibble(series)

luma <- series %>% filter(str_detect(title, "Luminal A")) %>% select(geo) %>% unlist
names(luma) <- NULL
healthy <- series %>% filter(str_detect(title, "Normal")) %>% select(geo) %>% unlist
names(healthy) <- NULL

matrix <- read_tsv(paste0("data/", set, "_series_matrix.txt"), 
                   skip = 68, col_names = F)
colnames(matrix) <- matrix[1, ]
colnames(matrix)[1] <- "Name" 
matrix  <- matrix[c(-1, -nrow(matrix)), ]  
matrix <- matrix %>% mutate_at(vars(matches("GSM")), as.numeric) %>% 
  select(Name, all_of(healthy), all_of(luma))

deg <- read_tsv(paste0("data/", set, "_luma_vs_healthy.txt"))
deg <- deg %>% select(ID, Gene.ID)
colnames(deg) <- c("Name", "Description")

matrix <- matrix %>% left_join(deg, by = "Name") %>% 
  select(Name, Description, everything())

write_tsv(matrix, paste0("data/", set, "_exprGSEA.txt"))

cls <- paste(paste(length(c(healthy, luma)), "2", "1", sep = "\t"), 
        paste("#", "healthy", "luma", sep = "\t"), 
        paste(c(rep.int(0, times = length(healthy)), rep.int(1, times = length(luma))), collapse = "\t"),
       sep = "\n")
write(cls, file = paste0("data/", set, "_GSEA.cls"))

#### SET 2
set <- sets[2]
series <- read_tsv(paste0("data/", set, "_series_matrix.txt"), 
                   skip = 31, n_max = 32, col_names = F)
series <- t(series)[, c(9, 1)]
colnames(series) <- c("title", "geo")
series <- series[-1, ]
series <- as_tibble(series)

luma <- series %>% filter(str_detect(title, "Luminal A")) %>% select(geo) %>% unlist
names(luma) <- NULL
healthy <- series %>% filter(str_detect(title, "Normal")) %>% select(geo) %>% unlist
names(healthy) <- NULL

matrix <- read_tsv(paste0("data/", set, "_series_matrix.txt"), 
                   skip = 64, col_names = F)
colnames(matrix) <- matrix[1, ]
colnames(matrix)[1] <- "Name" 
matrix  <- matrix[c(-1, -nrow(matrix)), ]  
matrix <- matrix %>% mutate_at(vars(matches("GSM")), as.numeric) %>% 
  select(Name, all_of(healthy), all_of(luma))

deg <- read_tsv(paste0("data/", set, "_luma_vs_healthy.txt"))
deg <- deg %>% select(ID, Gene.ID)
deg <- deg %>% mutate(ID = as.character(ID))
colnames(deg) <- c("Name", "Description")

matrix <- matrix %>% inner_join(deg, by = "Name") %>% 
  select(Name, Description, everything())

write_tsv(matrix, paste0("data/", set, "_exprGSEA.txt"))

cls <- paste(paste(length(c(healthy, luma)), "2", "1", sep = "\t"), 
             paste("#", "healthy", "luma", sep = "\t"), 
             paste(c(rep.int(0, times = length(healthy)), rep.int(1, times = length(luma))), collapse = "\t"),
             sep = "\n")
write(cls, file = paste0("data/", set, "_GSEA.cls"))


### SET TCGA
set <- sets[3]
luma_matrix <- read_tsv("data/luma_cpm10_arsyn.tsv")
healthy_matrix <- read_tsv("data/healthy_cpm10_arsyn.tsv") 
healthy <- colnames(healthy_matrix)[-1]
luma <- colnames(luma_matrix)[-1]
matrix <- luma_matrix %>% inner_join(healthy_matrix, by = "gene")
colnames(matrix)[1] <- "Name"
matrix <- matrix %>% mutate(Description = "NA") %>% 
  select(Name, Description, all_of(healthy), all_of(luma))

cls <- paste(paste(length(c(healthy, luma)), "2", "1", sep = "\t"), 
             paste("#", "healthy", "luma", sep = "\t"), 
             paste(c(rep.int(0, times = length(healthy)), rep.int(1, times = length(luma))), collapse = "\t"),
             sep = "\n")
