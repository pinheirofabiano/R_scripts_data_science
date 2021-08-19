
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("CEMiTool")
library(CEMiTool)

setwd("/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results")

dirOut <- "CEMiToll"
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T)
}


log2cpm_expression <- read.table(file = "log2cpm_expression.tsv", header = T, sep = "\t", row.names = 1)

expr_heart <- log2cpm_expression[, grep("heart", colnames(log2cpm_expression), fixed = TRUE, value = F)]


phenodata <- read.table(file = "phenodata.tsv", header = T, sep = "\t", row.names = 1)
phenodata <- cbind(rownames(phenodata), phenodata)
rownames(phenodata) <- NULL
colnames(phenodata) <- c("Sample", "Class", "Tissue")

phenodata_heart <- phenodata[which(phenodata$Tissue == "heart"),]

    cem_heart <- cemitool(expr_heart, annot = phenodata_heart, 
                  #gmt, 
                  #interactions, 
                  filter = TRUE, 
                  filter_pval = 0.1,
                  apply_vst = TRUE, eps = 0.1, cor_method = c("pearson"), cor_function = "cor", network_type = "unsigned", tom_type = "signed", set_beta = NULL, force_beta = TRUE, sample_name_column = "Sample", class_column = "Class", merge_similar = TRUE, rank_method = "median",
                  #ora_pval = 0.05, gsea_scale = TRUE, gsea_min_size = 15, gsea_max_size = 1000, min_ngen = 30, diss_thresh = 0.8, 
                  plot = TRUE, plot_diagnostics = TRUE, order_by_class = TRUE, center_func = "median", directed = FALSE, verbose = FALSE)
    
    
summary_heart <- mod_summary(cem_heart)
generate_report(cem_heart, directory = paste0(dirOut, "/Report_heart"))
write_files(cem_heart, directory = paste0(dirOut, "/Tables_heart"))
save_plots(cem_heart, "all", directory = paste0(dirOut, "/Plots_heart"))

#cem <- plot_gsea(cem)
#show_plot(cem, "gsea")
#cem <- plot_profile(cem)
#plots <- show_plot(cem, "profile")


# write.table(x = expr_heart, file = "expr_heart.tsv", quote = F, sep = "\t", row.names = T)


