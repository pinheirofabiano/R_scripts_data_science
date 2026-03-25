
colnames(tsalik_raw_counts_final)[1] <- "symbols"


df <- merge(herwanto_raw_gene_symbols, tsalik_raw_counts_final, by = "symbols")

unique_df <- df[!duplicated(df$symbols), ]

df2 <- unique_df

colnames(df2)

df2 <- df2[,-c(31:49)]
df2 <- df2[,-2]

write.table(df2, "~/Desktop/metanalysis_sepsis.tsv", sep="\t")

rownames(df2) = make.names(df2$symbols, unique=TRUE)
df2 <- df2[,-1]

colnames(df2) <- metanalysis_metadata$sample_ID

# download from Bioconductor

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")

# download from Github

remove.packages("sva")

install.packages("devtools")

devtools::install_github("zhangyuqing/sva-devel")

packageVersion("sva")

#calculate principal components for the NOT adjusted data

vsd <- vst(deg, blind = F)

plotPCA(vsd, intgroup = "condition")

# run combat-seq

count_matrix <- metanalysis_shock_healthy_only_t1_final
batches <- c(rep(1, 44), rep(2, 6), rep(3,4), rep(4, 13))
groups <- c(rep(1, 44), rep(2, 23))

metanalysis_shock_healthy_only_t1_final_adjusted <- sva::ComBat_seq(count_matrix, batch=batches)

metanalysis_shock_healthy_only_t1_final_adjusted <- as.data.frame(metanalysis_shock_healthy_only_t1_final_adjusted)


write.table(metanalysis_shock_healthy_only_t1_final_adjusted, "~/Desktop/metanalysis_braga_healthy_only_t1_final_adjusted.tsv", sep="\t")

#calculate principal components for the adjusted data

vsd <- vst(deg, blind = F)

plotPCA(vsd, intgroup = "condition")




