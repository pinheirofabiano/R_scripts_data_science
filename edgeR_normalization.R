
#Clear your workspace
rm(list=ls())

##Installing required libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages('dplyr')

#Load required libraries
library(edgeR)
library(dplyr)


##reading raw expression

getwd()
AMPs <- read.table(file = 'AMPs_counts_herwanto.tsv', 
                    header = TRUE,
                    row.names = 1)
AMPs <- AMPs[,-1]

##reading metadata (sample information, same names as in the raw expression table)
meta <- read.table('metadata_herwanto.tsv',
                    header = TRUE)


####Module 2 - Data cleaning and preparation####
#checking if column order of expr table is the same as row order of meta1
#identical(colnames(expr0), meta1$Sample_code)

#rearranging expression table column order to match sample information order
#expr1 <- expr0[,meta1$Sample_code]
#rm(expr0) #removing expr0

#checking again if column order of expr table is the same as row order of meta1
identical(colnames(AMPs[,2:29]), meta$sample)

#check for missing values
table(is.na(AMPs))

####Submodule - Normalizing counts####
#factors are categorical variables
factor(meta$condition)

#determine any factors to include in your analysis design
group1 <- as.factor(meta$condition)
group1

#creating the model matrix - With an intercept term or without an intercept term
#design1 <- model.matrix(~group1)
#design1

design2 <- model.matrix(~group1+0)
design2


#matrix choose (this only changes the way you define your comparisons later at the contrasts)
design <- design2
head(design,3)
rm(design2)

#create differential gene expression object
dge <- DGEList(counts=AMPs, group=group1 )

#Calculate normalization factors. Essential step for normalization!
norm <- calcNormFactors(dge, method = "TMM")

#create TMM normalized table

tmm.exp <- normalizeQuantiles(cpm(norm, log = F))
class(tmm.exp)
tmm.exp <- as.data.frame(tmm.exp)
tmm.exp$symbols <- AMPs$symbols
herwanto_tmm_final <- tmm.exp[,c(29, 1:28)]

write.table(herwanto_tmm_final, file = "herwanto_tmm_counts.tsv", sep='\t')

