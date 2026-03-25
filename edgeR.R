
#Setting up your environment####

#Clear your workspace
rm(list=ls())

##Installing required libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install('biomaRt')
BiocManager::install("ComplexHeatmap")
install.packages('circlize')
install.packages('openxlsx')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('FactoMineR')
install.packages('plotly')

#Load required libraries
library(edgeR)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(FactoMineR)
library(circlize)
library(ComplexHeatmap)
library(plotly)
library(biomaRt)

##reading raw expression
expr0 <- read.table(file = 'GSE205748_read_counts_PsA.csv', 
                    header = TRUE,
                    row.names = 1)

##reading metadata (sample information)
meta1 <- read.table('GSE205748_series_matrix_edit.txt',
                    header = TRUE)

##selecting rows and columns
#dataframe[row,column], i.e.:
meta1[1,] #first row
meta1[,1] #first column


####Module 2 - Data cleaning and preparation####
#checking if column order of expr table is the same as row order of meta1
identical(colnames(expr0), meta1$Sample_code)

#rearranging expression table column order to match sample information order
expr1 <- expr0[,meta1$Sample_code]
rm(expr0) #removing expr0

#checking again if column order of expr table is the same as row order of meta1
identical(colnames(expr1), meta1$Sample_code)

#check for missing values
table(is.na(expr1))

####Submodule - Normalizing counts####
#factors are categorical variables
factor(meta1$Tissue_type)

#determine any factors to include in your analysis design
group1 <- as.factor(meta1$Tissue_type)
group1

#creating the model matrix - With an intercept term or without an intercept term
design1 <- model.matrix(~group1)
design1

design2 <- model.matrix(~group1+0)
design2

#compare matrices
head(design1,3)
head(design2,3)

#matrix choose (this only changes the way you define your comparisons later at the contrasts)
design <- design2
head(design,3)
rm(design1,design2)

#create differential gene expression object
d <- DGEList(counts=expr1, group=group1 )

#Calculate normalization factors. Essential step for normalization!
d <- calcNormFactors(d)

#check current expression, las = 2 turns x axis labels vertically
boxplot(cpm(d), las=2)

#check with log values to eliminate the effect of outliers
boxplot(cpm(d, log=TRUE), las=2)

##Filtering lowly expressed genes - keep genes that have more than one count per million (cpm) for at least 5% of samples

# 1) cpm function normalizes counts. 
# 2) cpm(d) >1 checks if the gene expression is more than 1 cpm
# 3) rowSums adds the number of TRUE per Row (aka the number of samples with more than 1 cpm)
# 4) dim(d) return the dimensions of the dataset; dim(d)[1] is the number of genes, dim(d)[2] is the number of samples
# 5) ceiling rounds the given number to the closest higher integer
keep <- rowSums( cpm( d ) > 1 ) >= ceiling(0.05*dim(d)[2])
head(keep)

#keep only rows with TRUE
d1 <- d[keep, ]
#Check how many genes were kept
dim(d1)[1] #Number
(dim(d1)[1]/dim(d)[1])*100 #Percentage

#Create a new DGElist object with the new counts
d2 <- DGEList(counts = d1, group = group1)
#recalculate normalization factors
d2 <- calcNormFactors(d2)

#get the normalized counts
norm_d <- cpm(d2)

## get the log values of the data
logCPM <- cpm(d2, log=TRUE)

#export your normalized values
write.table(norm_d, file = "GSE205748_cpm.txt", sep='\t')
write.table(logCPM, file = "GSE205748_logcpm.txt", sep = '\t')



####Module 3 - Exploratory data analysis####


#create a directory to save your plots
dir.create("Plots")

##Export boxplots to pdf file
#open file
pdf(file = "Plots/GSE205748_boxplots.pdf", width = 16, height = 8)
#convert graphics window to 2 columns
par(mfrow=c(1,2))
#Check sample count distribution again with boxplot
boxplot(norm_d, las=2)
#with log values 
boxplot(logCPM, las=2)
#close graphics device to save pdf
dev.off()
#convert graphics window back to 1 column
par(mfrow=c(1,1))

####Submodule - Principal Component Analysis####

#To perform the PCA we need to use the logCPM values and transpose the table
tlogcpm <- t(logCPM)

#performing the PCA
pcacpm <- PCA(tlogcpm, 
              scale.unit = T, #We scale for better PCA results - controversy
              graph = F) #False as we will create our own graphs

#isolate the data for plotting 
data_pca <- as.data.frame(pcacpm$ind$coord)

#isolate variance explained by each PC
head(pcacpm$eig) #we need the second column
class(pcacpm$eig) # this is a matrix, not a data.frame
pca_perc <- pcacpm$eig[,2] #choosing the second column


ggplot(data_pca, 
       aes(data_pca[, 1], data_pca[,2],
           color = meta1$Tissue_type)) + 
  geom_point(key_glyph = "point") + 
  labs(title = 'PCA - GSE205748', 
       x = paste0("PC1(", round(pca_perc[1], 2), ")"), 
       y = paste0("PC2(", round(pca_perc[2], 2), ")"),
       color = 'Tissue type') + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))


#generalize the plotting function
PCx=1 #Change accordingly
PCy=2 #Change accordingly

#plot
g <- ggplot(data_pca, 
            aes(data_pca[, PCx], data_pca[,PCy],
            color = meta1$Tissue_type)) + 
  geom_point(key_glyph = "point") + 
  labs(title = 'PCA - GSE205748', 
       x = paste0("PC",PCx, "(", round(pca_perc[PCx], 2), ")"), 
       y = paste0("PC", PCy, "(", round(pca_perc[PCy], 2), ")"),
       color = 'Tissue type') + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))
print(g)

##export to pdf
#open file
pdf(file = 'Plots/PCA_GSE205748.pdf', width = 8, height = 6)
#print the plot
print(g)
#close the file
dev.off()


####Module 4 - Differential expression analysis####

####Submodule 1 - Dispersion and model fit####
# Using a common estimate across all genes.
#d3 <- estimateGLMCommonDisp( d2, design, verbose=TRUE) 

# Fitting an estimate based on the mean-variance trend across the dataset, such 
#that genes similar abundances have similar variance estimates (trended dispersion)
#d3 <- estimateGLMTrendedDisp(d2, design) 

# Computing a genewise dispersion (tagwise dispersion). Needs one of 
#the former two as prerequisite. Best for multifactorial analysis
#d3 <- estimateGLMTagwiseDisp(d2, design) 

#Performs all 3 of the above dispersion estimations
d3 <- estimateDisp( d2, design, verbose=TRUE)

#Fit your model
fit <- glmQLFit(d3, design)

####Submodule 2 - Contrast define and DEA####

#Check design matrix categories
head(design,3)

##Create Contrast parameter.
#PsA Lesional vs Healthy would be PsA_les - Healthy, therefore based on the column order of the design matrix:
contr <- c(-1,1,0) # this corresponds to group1Healthy*-1 + group1PsA_les*1 + group1PsA_uninv*0, therefore PsA_les - Healthy

#You can also make the contrast parameter using makeContrasts if that is convenient.
contr2 <- makeContrasts(group1PsA_les-group1Healthy, levels = design) 

##Run DE test
lrt <- glmQLFTest(fit, contrast = contr )

#Check the resulting object
str(lrt)

#check the gene results
head(lrt$table)

##Multiple comparisons adjustment
#False Discovery Rate Correction (Benjamini-Hochberg)
lrt$table$fdr <- p.adjust(lrt$table$PValue, method="BH")

##Results Filtering

#filtering for fdr <= 0.01 (strict) and absolute logFC >=1 (strict)
top <- lrt$table[c(lrt$table$fdr <= 0.01 & abs(lrt$table$logFC)>=1),]

#filtering for fdr <= 0.05 (normal) and absolute logFC >=0.58 (normal)
top2 <- lrt$table[c(lrt$table$fdr <= 0.05 & abs(lrt$table$logFC)>=0.58),]

#export results
write.xlsx(top, file = "DE_Results_GSE205748_FDR_0_01_logFC1.xlsx", rowNames=TRUE, colNames=TRUE)
write.xlsx(top2, file = "DE_Results_GSE205748_FDR_0_05_logFC0_58.xlsx", rowNames=TRUE, colNames=TRUE)

#save workspace
save.image("GSE205748_DE_results.RData")

####Module 5 - Result annotation####
#choose biomart version
mart <- useEnsembl(biomart = 'ensembl', 
                   dataset = 'hsapiens_gene_ensembl', 
                   version = 105) 

#search for annotation based on ensembl_gene_id identifiers
annot <- getBM(filters= "ensembl_gene_id", #which identifier you are using
               attributes= c("ensembl_gene_id", 
                             "description", 
                             "start_position", 
                             "end_position", 
                             "strand", 
                             "hgnc_symbol"), #which attributes you want to collect
               values=rownames(lrt$table), #the names of your genes
               mart= mart) #the mart you defined previously

#creating final dataframe
lrt$table$ensembl_gene_id <- rownames(lrt$table)

Final_version <-  lrt$table %>%
  left_join(annot, by = 'ensembl_gene_id')
head(Final_version)

write.xlsx(Final_version, file="GSE205748_DE_results_PsA_les_vs_healthy.xlsx", colNames=TRUE, rowNames=TRUE)


####Module 6 - Results Visualization####
####Submodule 1 - Volcano plot####

ggplot(data=Final_version, 
       aes(x=logFC, y=-log10(fdr))) + 
  geom_point() + 
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") + 
  theme_minimal()



##adding color to our genes according to their significance level
DE_res <- Final_version

DE_res <- DE_res %>%
  mutate(significance = case_when(
    logFC >=1 & fdr <= 0.01 ~ 'Upregulated',
    logFC <=-1 & fdr <= 0.01 ~ 'Downregulated',
    abs(logFC) < 1 | fdr > 0.01 ~ 'Not significant'
  ))

##plotting again, this time with color
#opening file
pdf("Plots/VolcanoPlot_GSE205748.pdf", height = 8, width = 8)
#preparing the plot
g <- ggplot(DE_res,
       aes(x = logFC,
           y = -log10(fdr),
           color = significance,
           label = hgnc_symbol) # we will use rownames with plotly below, it is not required here for the volcano
       ) + 
  geom_point() + 
  geom_vline(xintercept = c(-1,1), color = 'red') + 
  geom_hline(yintercept = -log10(0.01), color = 'red') + 
  labs(x = 'LogFC', y = '-Log10(FDR)',color = 'Significance') + 
  scale_color_manual(values = c('blue','grey90','red'))+
  theme_minimal() 

#printing
print(g)

#closing the file
dev.off()

#Interactive volcano plot for better identification of DE significant genes
library(plotly)
print(g)
ggplotly()

####Submodule 2 - Heatmap####
#Select the top 50 upregulated genes
Up50 <- DE_res %>%
  slice_max(order_by = logFC, n =50)

#Select the top 50 downregulated genes
Down50 <- DE_res %>%
  slice_min(order_by = logFC, n =50)

#bind dataframes
Top100 <- bind_rows(Up50, Down50)

##Isolate the expression of top 100
logCPM_100 <- logCPM[Top100$ensembl_gene_id,] 


#save heatmap to a file
pdf('Plots/GSE205748_heatmap_top100_DE.pdf', width = 8, height =8)

#Create the heatmap
Heatmap(logCPM_100,
        row_labels = Top100$hgnc_symbol, 
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7), 
        heatmap_legend_param = list(title = "LogCPM\nexpression"), 
        top_annotation = HeatmapAnnotation(Condition = meta1$Tissue_type, 
                                           which = 'column', 
                                           col = list(Condition = c(PsA_les = 'turquoise4',
                                                                    PsA_uninv = 'red3',
                                                                    Healthy = 'green2')
                                                      )
                                           )
        )


#close the file
dev.off()

##save final workspace
date <- Sys.Date()
save.image(file = paste0(date,"_GSE205748_DE_analysis_complete.RData"))

#save session info to a text file for reproducibility purposes
sink(file = 'session_info_GSE205748_DE_analysis_complete.txt')
sessionInfo()
sink()


