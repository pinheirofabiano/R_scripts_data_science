
install.packages("BiocManager")

BiocManager::install(c("tidyverse", "BioNERO"))

library(tidyverse)
library(BioNERO)

df_final <- zhang_pneumonia_TPMs

#replace NCBI sample ID for project sample IDs
final_df <- df |>
  select(pancreatitis_metadata$GEO_ID)

colnames(final_df)

colnames(df_final) <- lupus_metadata$sample

final_df <- cbind(df$GeneID, final_df)
colnames(final_df)[1] <- "GeneID"

#order columns alphabetically
#final_df3 <- final_df2[,order(colnames(final_df2))]

# check male:female ratio in the healthy and disease groups

which(rownames(df) == "XIST")
df[40056,]

females <- c("healthy01", "healthy02", "healthy04","healthy05","healthy06","healthy07",
                       "healthy08","healthy10","healthy11","healthy12","healthy13","healthy14",
                       "healthy17","healthy19","healthy20","healthy21","healthy22",
                       "healthy23","healthy24","healthy26","healthy27","healthy29")
colnames(final_df)

df_final2 <- df_final[,-c(5, 59, 72, 76, 89, 98)]

df_final <- df_final2
final_df <- final_df[,-c(1,2)]

lupus_metadata$sample <- rownames(lupus_metadata)

colnames(df_final) <- lupus_metadata$sample

#save final 
write.table(df_final, "~/Desktop/lupus_TPM_final.tsv", sep="\t")

df_final <- zhang_pneumonia_TPMs

# remove rows with NAs in the column "symbols"
final_df <- df |>
 drop_na(symbols)

# rename row names with Symbols (first column)
rownames(df_final) = make.names(df_final$symbols, unique=TRUE)
df_final <-  df_final[,-c(1,2)]

df_final <- shock_TPM



# remove non-expressed genes
df_final <- remove_nonexp(df_final, method = "median", min_exp = 2) 

# filter low expressed genes
df_final <- filter_by_variance(df_final, n = 4000)
?filter_by_variance

#PC correction
df_final <- PC_correction(df_final)
df_final <- as.data.frame(df_final)

shock_metadata <- GSE310929_shock_metadata[,c(1,8)]
rownames(shock_metadata) = make.names(shock_metadata$symbols, unique=TRUE)
shock_metadata2 <- shock_metadata[,-1]
shock_metadata2 <- as.data.frame((shock_metadata2))
rownames(shock_metadata2) = make.names(shock_metadata$Sample.ID, unique=TRUE)

pancreatitis_metadata_final <- pancreatitis_metadata_final[,-c(1,2)]
nomes <- pancreatitis_metadata_final$patient_ID
pancreatitis_metadata_final <- pancreatitis_metadata_final[,-1]
pancreatitis_metadata_final <- as.data.frame(pancreatitis_metadata_final)
rownames(pancreatitis_metadata_final) <- nomes
colnames(pancreatitis_metadata_final) <- "class"
class(pancreatitis_metadata_final)


# Plot PCA
?plot_PCA
p_pca <- plot_PCA(
  df_final,
  metadata = metadata1,
  metadata_cols = NULL,
)

p_pca

# Find optimal beta power to which correlation coefficients will be raised
# Herwanto (5000 nodes, remove non exp = 2, beta power = 18), covid (1000 nodes, 
#remove non exp = 5, beta power = 20),
#pancreatitis (800 nodes, remove non exp = 5, PC correction, beta power = 13)

sft <- SFT_fit(
  df_final, 
    net_type = "signed",
    cor_method = "pearson"
)

sft$power
sft$plot

# Infer a GCN
gcn <- exp2gcn(
  df_final, 
  net_type = "signed",
  SFTpower = 12,
  cor_method = "pearson"
)


names(gcn)
head(gcn$genes_and_modules)
amps <- c("CAMP", "HEPCIDIN", "HAMP", "S100A8", "S100A9", "DEFA1", "DEFA1B", "DEFA3",
          "DEFA4", "DEFB1", "DEFB4A", "DEFB103A", "DEFB104A")
amps_genes <- which(gcn$genes_and_modules$Genes %in% amps)

amps_genes_final <- gcn$genes_and_modules[amps_genes,]

pneumonia_module_saddlebrown <- which(gcn$genes_and_modules$Modules == "saddlebrown")

pneumonia_module_saddlebrown_final <- gcn$genes_and_modules[pneumonia_module_saddlebrown,]

length(pancreatitis_module_darkgreen_final$Genes)

shock_module_black <- which(gcn$genes_and_modules$Modules == "black")

shock_module_black_final <- gcn$genes_and_modules[shock_module_black,]

#save AMPs vs modules list
write.table(pneumonia_module_saddlebrown_final, "~/Desktop/pneumonia_module_saddlebrown_final.tsv", sep="\t")

# module grey is a bin!!
plot_ngenes_per_module(gcn)

?module_trait_cor
class(metadata_herwanto)
colnames(shock_metadata2) <- "class"
metadata1$metadata1 <- factor(metadata1$metadata1, levels=c("healthy", "pneumonia"))


# Calculating module-trait correlations
me_trait <- module_trait_cor(
  exp = df_final,
  MEs = gcn$MEs,
  metadata = metadata1,
  cor_method = "pearson"
)

# Inspecting the results
head(me_trait)

plot_module_trait_cor(me_trait)

?plot_expression_profile

?plot_expression_profile
plot_expression_profile(
  exp = df_final, 
  net = gcn,
  modulename = "black",
  metadata = shock_metadata2
)


# Functional analysis

install.packages("enrichR")
library(enrichR)

setEnrichrSite("Enrichr")
websiteLive <- TRUE

dbs <- listEnrichrDbs()

dbs <- c("Reactome_Pathways_2024", "KEGG_2021_Human", "GO_Biological_Process_2025",
         "GO_Molecular_Function_2025")

if (websiteLive) {
  enriched <- enrichr(pneumonia_module_saddlebrown_final$Genes, dbs)
}

#Reactome database 2024

if (websiteLive) enriched[["Reactome_Pathways_2024"]]

if (websiteLive) plotEnrich(enriched[[1]], showTerms = 20, numChar = 80, y = "Count", orderBy = "P.value")

#KEGG database 2021

if (websiteLive) enriched[["KEGG_2021_Human"]]

if (websiteLive) plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

# Gene Ontology Biological Process 2025

if (websiteLive) enriched[["GO_Biological_Process_2025"]]

if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 80, y = "Count", orderBy = "P.value")

# Gene Ontology Molecular Function 2025

if (websiteLive) enriched[["GO_Molecular_Function_2025"]]

if (websiteLive) plotEnrich(enriched[[4]], showTerms = 20, numChar = 80, y = "Count", orderBy = "P.value")


# identify hubs

hubs <- get_hubs_gcn(exp = df_final, net = gcn)

head(hubs)

# network visualization

edge_list <- get_edge_list(
  net = gcn,
  module = "bisque4",
  filter = TRUE
)

head(edges)

plot_gcn(
  edgelist_gcn = edge_list,
  net = gcn,
  color_by = "module",
  hubs = hubs
)

#network statistics

install.packages("igraph")
library(igraph)

## Graph from adjacency matrix undirected and weighted

?graph_from_adjacency_matrix
gg <- graph_from_adjacency_matrix(gcn$adjacency_matrix, "undirected", diag = F, weighted = T)

V(gg)
E(gg)
mean_distance(gg)
transitivity(gg)

gg2 <- graph_from_data_frame(edge_list, directed = F)
plot(gg2)


# degree
c1 <- strength(gg)
c1 <- as.data.frame(c1)
c1[,2] <- rownames(c1)
amps_c1 <- which(c1$V2 %in% amps)
c1[amps_c1,]

# closeness
c2 <- closeness(gg)
c2 <- as.data.frame(c2)
c2[,2] <- rownames(c2)
amps_c2 <- which(c2$V2 %in% amps)
c2[amps_c2,]



# betweenness
c3 <- betweenness(gg)
c3 <- as.data.frame(c3)
c3[,2] <- rownames(c3)
amps_c3 <- which(c3$V2 %in% amps)
c3[amps_c3,]





