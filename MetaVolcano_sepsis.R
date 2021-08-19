
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MetaVolcanoR")

install.packages("data.table")
install.packages("tidyverse")

library(MetaVolcanoR)
library(data.table)
library(tidyverse)

setwd("/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/projeto_AUTOPSIAS/data_analysis/results/")

kidney.sepsis1 <- fread(file = "DEGs_IndividualvsAllControls_per_sex_FILTERED/ALL_COUNTS/edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis1.kidney.tsv", header = T)
kidney.sepsis2 <- fread(file = "DEGs_IndividualvsAllControls_per_sex_FILTERED/ALL_COUNTS/edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis2.kidney.tsv", header = T)
kidney.sepsis3 <- fread(file = "DEGs_IndividualvsAllControls_per_sex_FILTERED/ALL_COUNTS/edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis3.kidney.tsv", header = T)
kidney.sepsis4 <- fread(file = "DEGs_IndividualvsAllControls_per_sex_FILTERED/ALL_COUNTS/edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis4.kidney.tsv", header = T)
kidney.sepsis5 <- fread(file = "DEGs_IndividualvsAllControls_per_sex_FILTERED/ALL_COUNTS/edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis5.kidney.tsv", header = T)
kidney.sepsis6 <- fread(file = "DEGs_IndividualvsAllControls_per_sex_FILTERED/ALL_COUNTS/edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis6.kidney.tsv", header = T)
kidney.sepsis7 <- fread(file = "DEGs_IndividualvsAllControls_per_sex_FILTERED/ALL_COUNTS/edgeR_all_FDR_0.1_log2FC_1_infected_VS_control_individual_sepsis7.kidney.tsv", header = T)

kidney.sepsis <- list(kidney_sepsis1 = kidney.sepsis1,
                     kidney_sepsis2 = kidney.sepsis2,
                     kidney_sepsis3 = kidney.sepsis3,
                     kidney_sepsis4 = kidney.sepsis4,
                     kidney_sepsis5 = kidney.sepsis5,
                     #kidney_sepsis6 = kidney.sepsis6,
                     kidney_sepsis7 = kidney.sepsis7)

# DEGs_combined <- fread("DEGs_IndividualvsAllControls_per_sex_FILTERED/All_edgeR_combined.txt", header = T)
# DEGs_combined_FDR0.1 <- DEGs_combined[, paste0(organ, "|", "sepsis", [0-9], "|" logFRD)


dirOut <- "Meta_Volcano/Metacombs"
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T)
}


# Combining-approach

degs_comb_kidney <- combining_mv(diffexp= kidney.sepsis,
                               pcriteria='PValue', 
                               foldchangecol='logFC',
                               genenamecol='gene_ID',
                               geneidcol=NULL,
                               metafc='Median',
                               metathr=0.01, 
                               collaps=TRUE,
                               jobname="kidney",
                               outputfolder="Meta_Volcano/Metacombs",
                               draw='HTML')

###############
DEGs_kidney <- degs_comb_kidney@metaresult
kidney_input <- degs_comb_kidney@input


inputname_comb <- degs_comb_kidney@inputnames
inputname_comb <- rep(inputname_comb, each = 2)

for (i in 2:length(colnames(kidney_input))) {
  colnames(kidney_input)[i] <- paste(colnames(kidney_input)[i], inputname_comb[i-1], sep = "_")
}

DEGs_kidney$fdr <- p.adjust(DEGs_kidney$metap, method = "fdr")

################

# Upregulated DEGs

# upregulated_DEGs_kidney <- degs_comb_kidney@metaresult %>%
#   filter(metafc > 0)
# 
# # Downregulated DEGs
# 
# downregulated_DEGs_kidney <- degs_comb_kidney@metaresult %>%
#   filter(metafc < 0)


write.table(DEGs_kidney, file = ("Meta_Volcano/Metacombs/DEGs_kidney.tsv"), quote = F, sep = "\t", row.names = F)

write.table(kidney_input, file = ("Meta_Volcano/Metacombs/kidney_input.tsv"), quote = F, sep = "\t", row.names = F)


# Plot MetaVolcano

degs_comb_kidney@MetaVolcano

# save Plot as TIFF!!


