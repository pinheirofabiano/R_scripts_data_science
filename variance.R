
# Set workspace for access the input files

setwd("/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results")

# install tidyverse and matrixStats

install.packages("tidyverse")
install.packages("matrixStats")
library(tidyverse)
library(matrixStats)

# Loading normalized expression file
log2cpm <- read.table(file = "log2cpm_expression.tsv", header = T, sep = "\t", row.names = 1)

groups <- c("control", "sepsis")
tissues <- c("cortex", "hippocampus", "heart", "lung", "kidney", "colon")

results <- NULL

for(group in groups){
for(tissue in tissues){
  log2cpm_per_tissue <- log2cpm %>%
    select(starts_with(group) & ends_with(tissue))

  variance_per_tissue<- rowVars(as.matrix(log2cpm_per_tissue[sapply(log2cpm_per_tissue, is.numeric)]))
  
  results[[paste0(group, ".", tissue)]] <- variance_per_tissue
  
  }
}

results<- as.data.frame(results)

rownames(results) <- rownames(log2cpm)


write.table(x = results, file = "log2cpm_with_variance_per_tissue.tsv", quote = F, sep = "\t")

  
