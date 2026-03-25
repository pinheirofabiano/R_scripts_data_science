if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RNAAgeCalc")

library("RNAAgeCalc")


df_final <- df

rownames(df_final) = make.names(df_final$symbols, unique=TRUE)
df_final <-  df_final[,-1]

fpkm_file <- count2FPKM(df_final, genelength = NULL, idtype = "SYMBOL")

TAE_sepsis_hancock <- predict_age(
  df_final,
  tissue = "blood",
  exprtype = "counts",
  idtype = "ENTREZID",
  stype = "all",
  signature = NULL,
  genelength = NULL,
  chronage = NULL,
  maxp = NULL
)


TAE_sepsis_hancock$chrono_Age <- metadata_50sepsis$Age

TAE_sepsis_hancock$delta <- TAE_sepsis_hancock$RNAAge - TAE_sepsis_hancock$chrono_Age

library(tidyr)
library(dplyr)

# Sample Data
df <- data.frame(combined = c("S1_T1", "S2_T1", "S3_T2"), value = c(10, 15, 20))

# Split column "combined" into "SampleID" and "Timepoint" using "_"
df <- TAE_sepsis %>%
  separate(col = IDs, into = c("sampleID", "timepoint"), sep = "_", convert = T)



# Create the spaghetti plot
ggplot(df, aes(x = df$timepoint, y = df$RNAAge, group = df$sampleID)) +
  geom_line() + geom_smooth(aes(group = 1), method = "loess", color = "red", se = FALSE) +
  labs(title = "Biological Age Over Time in Septic Patients - Whole Blood",
       x = "Time Points (16h, 48h and Day 7)",
       y = "RNA Age (years)") +
  theme_minimal()


df2 <- df[, c("sampleID", "timepoint", "RNAAge")]

# one-way ANOVA

df2$timepoint <- as.factor(df2$timepoint)

anova_sepsis <- aov(df2$RNAAge ~ df2$timepoint + Error(df2$sampleID/df2$timepoint), data = df2)

summary(anova_sepsis)

#save final count matrix
write.table(TAE_sepsis_hancock, "~/Desktop/TAE_sepsis_hancock.tsv", sep="\t")







  