DEG <- read.delim(file = "/Users/fabiano/Desktop/tsalik_count_matrix_final.tsv")
which(tsalik_count_matrix_final$Gene.symbols=="GAPDH")
AMPs <- tsalik_count_matrix_final[c(2315, 2317, 7115, 34158, 16459, 16461, 16463, 16455, 16453, 16499, 16498, 16476, 16496, 16478, 23386),]
write.table(AMPs, "/Users/fabiano/Desktop/AMPs_counts_tsalik.tsv")

AMPs <- AMPs_counts_tsalik

DEFA4_sirs_raw <- unlist(AMPs[8, c(3:20)], use.names = F)
DEFA4_sepsis_raw <- unlist(AMPs[8, c(21:50)], use.names = F)

GAPDH_sirs <- unlist(AMPs[15, c(3:20)], use.names = F)
GAPDH_sepsis <- unlist(AMPs[15, c(21:50)], use.names = F)

DEFA4_sirs <- DEFA4_sirs_raw/GAPDH_sirs
DEFA4_sepsis <- DEFA4_sepsis_raw/GAPDH_sepsis

#calcular a média dos controles
mean_sirs <- mean(DEFA4_sirs)

# calcular fold change

DEFA4_sirs_FC <- DEFA4_sirs/mean_sirs
DEFA4_sepsis_FC <- DEFA4_sepsis/mean_sirs

#calculate max length of vectors
max_length <- max(length(healthy), length(sepsis))

#set length of each vector equal to max length
length(healthy) <- max_length                      
length(sepsis) <- max_length 

#wilcox test


healthy <- TAE_sepsis_meta_braga$delta[1:44]
sepsis <- TAE_sepsis_meta_braga$delta[45:67]

total <- cbind(healthy, sepsis)
total <- as.data.frame(total)

wilcox.test(TAE_sepsis_meta_braga$delta[1:44], TAE_sepsis_meta_braga$delta[45:67], data = total)


TAE_sepsis_meta_braga$chrono_Age <- GSE310929_shock_metadata_only_t1$Gender

TAE_sepsis_meta_braga$delta <- TAE_sepsis_meta_braga$chrono_Age - TAE_sepsis_meta_braga$RNAAge
  
  
# Load ggplot2
install.packages("ggplot2")
library(ggplot2)

# boxplot

ggplot(total) +
  geom_boxplot(aes(x="HEALTHY", y=healthy, fill="blue", alpha=0.2, outlier.color = "black")) + 
  geom_boxplot(aes(x= "SEPSIS", y=sepsis, fill="blue", alpha=0.2, outlier.color = "black")) + 
  ylab("delta") + xlab("") +
  coord_flip() + theme(panel.grid.major = element_blank(), legend.position = "none", panel.grid.minor = element_blank())


