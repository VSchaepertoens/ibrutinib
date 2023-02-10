#Analysis of metabolomics data
#Pre-processing of data https://github.com/biosustain/snakemake_UmetaFlow
#Data-cleaning of data https://github.com/Functional-Metabolomics-Lab/Statistical-analysis-of-non-targeted-LC-MSMS-data/blob/main/Combined_Notebooks/Statistical_Analysis_Summerschool_22_python.ipynb

#loading libraries
library(tidyverse)

#Loading dataset----------------------------------------------------------------
data_matrix_scaled <- read.table("data/Scaled-feature-table-for-statistical-analysis(3).tsv", 
                                 header = TRUE, 
                                 check.names = FALSE)

#Data restructuring-------------------------------------------------------------

data_matrix_scaled_transposed <- t(data_matrix_scaled)

data_matrix <- matrix(as.numeric(data_matrix_scaled_transposed[-c(1:3),]),
                      ncol = 18)
rownames(data_matrix) <- rownames(data_matrix_scaled_transposed[4:676,])

#Load meta data
pheno_data <- as.data.frame(data_matrix_scaled[,c(1:3)])
pData <- read.csv("data/pData_removedBlanks_removedPool.csv")
pheno_data$replicate <- pData$repetition
pheno_data$experiment <- pData$experiment
pheno_data$sample_name <- paste(pheno_data$ATTRIBUTE_Sample,
                                pheno_data$replicate, 
                                pheno_data$experiment, 
                                sep = "_")

colnames(data_matrix) <-  rownames(pheno_data)
rownames(pheno_data) == colnames(data_matrix) 

#Calculate row means of technical replicates------------------------------------------
data_matrix_scaled$experiment <- pheno_data$experiment
data_matrix_scaled <- data_matrix_scaled[,-c(1:3)]

data_matrix_scaled_means <- data_matrix_scaled %>%
  group_by(experiment) %>%
  summarise_all(mean)

data_matrix_scaled_means$treatment <- c('control','control','control','ibrutinib','ibrutinib','ibrutinib')

#Final datamatrix and pheno data
data_matrix_scaled_transposed <- t(data_matrix_scaled_means)
data_matrix <- matrix(as.numeric(data_matrix_scaled_transposed[2:674,]), ncol = 6)


meta_data <- data.frame(experiment = data_matrix_scaled_means$experiment,
                        treatment = data_matrix_scaled_means$treatment)

meta_data$sample_name <- paste(meta_data$treatment,
                               meta_data$experiment, 
                               sep = "_")
rownames(data_matrix) <- rownames(data_matrix_scaled_transposed[2:674,])
colnames(data_matrix) <- meta_data$sample_name
rownames(meta_data) <- meta_data$sample_name
rownames(meta_data) == colnames(data_matrix) 

pheno_data <- meta_data



