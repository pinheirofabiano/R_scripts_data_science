
browseVignettes("CEMiTool")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("CEMiTool")
library(CEMiTool)

BiocManager::install("edgeR")
library(edgeR)

setwd("/Users/fabianopinheiro/Desktop")

dirOut <- "CEMiToll"
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T)
}

###########################################################
# count matrix preparation
# open dataframe
herwanto <- read.table("~/Desktop/herwanto_with_symbols.tsv", header = T, sep = "\t")

# drop samples with NA in the first column
herwanto <- filter(herwanto, !is.na(herwanto$symbols))

# assign 0 to NAs
herwanto <- replace(herwanto,is.na(herwanto),0)

# remove duplicated rows
herwanto <- subset(herwanto, !duplicated(herwanto$symbols))
duplicated(herwanto)

# assign row names to column 1
rownames(herwanto) <- herwanto$symbols
herwanto <- herwanto[,-1]
write.table(herwanto, file = "herwanto_count_matrix_with_rownames.tsv", quote = F, sep = "\t", row.names = T)

#######################################################

#normalization and define gender

herwanto <- read.table("~/Desktop/herwanto_count_matrix_with_rownames.tsv", header = T, sep = "\t")

herwanto_logCPM <- cpm(herwanto, log = T)
herwanto_logCPM <- as.data.frame(herwanto_logCPM)

#phenodata
phenodata <- read.table(file = "~/Desktop/phenodata_herwanto.txt", header = T, sep = "\t")
phenodata <- t(phenodata)
colnames(phenodata)  <- c("Sample", "Class")
rownames(phenodata) <- NULL
phenodata <- data.frame(phenodata)
class(phenodata)

    cem <- cemitool(expr = herwanto_logCPM, annot = phenodata,
                  #gmt, 
                  #interactions, 
                  filter = TRUE, 
                  filter_pval = 0.1,
                  apply_vst = F, eps = 0.1, cor_method = c("pearson"), cor_function = "cor", network_type = "unsigned", tom_type = "signed", set_beta = NULL, force_beta = T, sample_name_column = "Sample", class_column = "Class", merge_similar = TRUE, rank_method = "mean",
                  #ora_pval = 0.05, gsea_scale = TRUE, gsea_min_size = 15, gsea_max_size = 1000, min_ngen = 30, diss_thresh = 0.8, 
                  plot = TRUE, plot_diagnostics = TRUE, order_by_class = TRUE, center_func = "mean", directed = FALSE, verbose = FALSE)
    
summary_herwanto_logCPM <- mod_summary(cem)
generate_report(cem, directory = paste0(dirOut, "/Report_herwanto_logCPM"))
write_files(cem, directory = paste0(dirOut, "/Tables_herwanto_logCPM"))
save_plots(cem, "all", directory = paste0(dirOut, "/Plots_herwanto_logCPM"))

#cem <- plot_gsea(cem)
#show_plot(cem, "gsea")
#cem <- plot_profile(cem)
#plots <- show_plot(cem, "profile")


write.table(herwanto, file = "herwanto_count_matrix_final.csv", quote = F, sep = ",", row.names = T)


