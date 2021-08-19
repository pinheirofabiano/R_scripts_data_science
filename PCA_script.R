rm(list = ls())
#!/usr/bin/env Rscript
# library("optparse")
# 
# option_list = list(
# 		   make_option(c("-e", "--expression-file"), type="character", default=NULL, help="Expression file name", metavar="character"),
# 		   make_option(c("-o", "--out"), type="character", default="out.txt", help="output file name [default= %default]", metavar="character"))
# 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)

library("jsonlite")

params <- fromJSON("~/sandbox/working/2020_Sandra/v4/analysis/scripts/config.json")

if(!dir.exists(params$workdir)) {
  dir.create(params$workdir, recursive = T)
}

setwd(params$workdir)

norm.exp  <- read.delim(file = paste0("analysis/geneCounts/", stringi::stri_trans_tolower(params$normalizeMethod), "_expression.tsv"), header = T, stringsAsFactors = F)
phenodata <- read.delim(file = "analysis/geneCounts/phenodata.tsv", header = T)

phenodata <- phenodata[ phenodata$Sample %in% colnames( norm.exp ), ]

norm.exp <- norm.exp[ ,c( "Symbol", as.character( phenodata$Sample ) ) ]
phenodata <- phenodata[ phenodata$Sample %in% colnames( norm.exp ), ]

rownames( norm.exp ) <- as.character( norm.exp$Symbol )
norm.exp$Symbol <- NULL

# transpose the data to have variables (genes) as columns
data_for_PCA <- t(norm.exp)

## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned  

# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)

suppressMessages(library(ggplot2))
pdf(file = paste0(params$PCA, "PCA_PropExplainedVariance_", stringi::stri_trans_tolower(params$normalizeMethod), ".pdf"), width = 2, height = 3)
eig_pc <- eig_pc[1:10]
ggplot(data = data.frame(eigen = eig_pc, PC = paste0("Dim ", 1:length(eig_pc))), aes(x = reorder(PC, -eigen), eigen)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Dimensions") +
  ylab("Explained var.")
dev.off()

## calculate MDS
mds <- as.data.frame(cmdscale(dist(data_for_PCA), k=3)) # Performs MDS analysis
colnames(mds) <- c("Dim1", "Dim2", "Dim3")

mds$Condition <- phenodata$Class

# mds_uninfected <- mds[mds$Condition == "Uninfected",]
# mds_La <- mds[mds$Condition == "La",]
# mds_Lb <- mds[mds$Condition == "Lb",]
# mds_Li <- mds[mds$Condition == "Li",]

suppressMessages(library(ggalt))

pdf(file=paste0(params$PCA, "/PCA_PC1_vs_PC2_", stringi::stri_trans_tolower(params$normalizeMethod), ".pdf"), width = 4, height = 4)
ggplot(mds, aes(Dim1, Dim2, color=Condition, shape=Condition)) +
  geom_point(size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8)) +
  xlab(paste0("Dim 1 (", format(eig_pc[1], digits = 4), "% explained var.)")) +
  ylab(paste0("Dim 2 (", format(eig_pc[2], digits = 4), "% explained var.)")) +
  coord_cartesian(xlim = 1.2 * c(min(mds$Dim1), max(mds$Dim1)),
                  ylim = 1.2 * c(min(mds$Dim2), max(mds$Dim2))) +   # change axis limits
  # geom_encircle(data = mds_uninfected, aes(x=Dim1, y=Dim2)) +   # draw circles
  # geom_encircle(data = mds_La, aes(x=Dim1, y=Dim2)) +
  # geom_encircle(data = mds_Lb, aes(x=Dim1, y=Dim2)) +
  # geom_encircle(data = mds_Li, aes(x=Dim1, y=Dim2)) +
  coord_fixed() +
  xlim(c(ceiling(min(mds$Dim1))-1,ceiling(max(mds$Dim1))+1)) +
  ylim(c(ceiling(min(mds$Dim2))-1,ceiling(max(mds$Dim2))+1))
dev.off()
