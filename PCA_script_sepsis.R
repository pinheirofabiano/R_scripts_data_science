

counts  <- read.delim(file = '/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results/counts.tsv', header = T, sep = "\t", stringsAsFactors = F)
phenodata <- read.delim(file = "/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results/phenodata.tsv", header = T, sep = "\t", stringsAsFactors = F)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)
counts$meanG <- apply(counts[ , 2:ncol(counts) ], 1, mean )
counts <- counts[ order(counts[ , 'Symbol' ], counts[ , 'meanG' ] ), ]
counts <- counts[ !duplicated(counts$Symbol), ]
rownames(counts) <- counts$Symbol
counts <- counts[ , 2:( ncol(counts) - 1 ) ]


cpm.exp <- normalizeQuantiles(log2(cpm(counts)+1))
write.table(x = data.frame(Symbol = rownames(cpm.exp), cpm.exp), file = "/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results/log2cpm_expression.tsv", quote = F, sep = "\t", row.names = F)




phenodata <- phenodata[ phenodata$Sample %in% colnames(cpm.exp), ]


# transpose the data to have variables (genes) as columns
data_for_PCA <- t(cpm.exp)

## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k = 3, eig = TRUE)
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned  

# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)

suppressMessages(library(ggplot2))




pdf(file = paste0("/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results/log2cpm_expression", ".pdf"), width = 2, height = 3)
# dev.off()

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
mds$Tissue <- phenodata$Tissue

mds_control <- mds[mds$Condition == "control",]
mds_sepsis <- mds[mds$Condition == "sepsis",]
mds_colon <- mds[mds$Tissue == "colon", ]
mds_cortex <- mds[mds$Tissue == "cortex",]
mds_heart <- mds[mds$Tissue == "heart",]
mds_hippocampus <- mds[mds$Tissue == "hippocampus",]
mds_kidney <- mds[mds$Tissue == "kidney",]
mds_lung <- mds[mds$Tissue == "lung",]

suppressMessages(library(ggalt))

pdf(file=paste0("/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results/log2cpm_expression", ".pdf"), width = 4, height = 4)

ggplot(mds, aes(Dim1, Dim2, color=Tissue, shape=Condition)) +
  geom_point(size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8)) +
  xlab(paste0("Dim 1 (", format(eig_pc[1], digits = 4), "% explained var.)")) +
  ylab(paste0("Dim 2 (", format(eig_pc[2], digits = 4), "% explained var.)")) +
  coord_cartesian(xlim = 1.2 * c(min(mds$Dim1), max(mds$Dim1)),
                  ylim = 1.2 * c(min(mds$Dim2), max(mds$Dim2))) +   # change axis limits
  geom_encircle(data = mds_control, aes(x=Dim1, y=Dim2)) +
  geom_encircle(data = mds_sepsis, aes(x=Dim1, y=Dim2)) +           # draw circles
  geom_encircle(data = mds_colon, aes(x=Dim1, y=Dim2)) +
  geom_encircle(data = mds_cortex, aes(x=Dim1, y=Dim2)) +
  geom_encircle(data = mds_heart, aes(x=Dim1, y=Dim2)) +
  geom_encircle(data = mds_hippocampus, aes(x=Dim1, y=Dim2)) +
  geom_encircle(data = mds_kidney, aes(x=Dim1, y=Dim2)) +
  geom_encircle(data = mds_lung, aes(x=Dim1, y=Dim2)) +
  coord_fixed() +
  xlim(c(ceiling(min(mds$Dim1))-1,ceiling(max(mds$Dim1))+1)) +
  ylim(c(ceiling(min(mds$Dim2))-1,ceiling(max(mds$Dim2))+1))

ggplot(mds_lung, aes(Dim1, Dim2, color=Condition, shape=Condition)) +
  geom_point(size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 8)) +
  xlab(paste0("Dim 1 (", format(eig_pc[1], digits = 4), "% explained var.)")) +
  ylab(paste0("Dim 2 (", format(eig_pc[2], digits = 4), "% explained var.)")) +
  coord_cartesian(xlim = 1.2 * c(min(mds_lung$Dim1), max(mds_lung$Dim1)),
                  ylim = 1.2 * c(min(mds_lung$Dim2), max(mds_lung$Dim2))) +   # change axis limits
  coord_fixed() +
  xlim(c(ceiling(min(mds_lung$Dim1))-1,ceiling(max(mds_lung$Dim1))+1)) +
  ylim(c(ceiling(min(mds_lung$Dim2))-1,ceiling(max(mds_lung$Dim2))+1)) + ggtitle('LUNG')

# tissues = colon, cortex, heart, hippocampus, kidney, lung
  
dev.off()
