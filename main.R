#Analysis of metabolomics data
#Pre-processing of data https://github.com/biosustain/snakemake_UmetaFlow
#Cleaning of data https://github.com/Functional-Metabolomics-Lab/Statistical-analysis-of-non-targeted-LC-MSMS-data/blob/main/Combined_Notebooks/Statistical_Analysis_Summerschool_22_python.ipynb
#data is not log2 transformed

#loading libraries
library(tidyverse)
library(pheatmap)
library(limma)
library(data.table)


#Loading dataset----------------------------------------------------------------
data_matrix_scaled <- read.table("data/RP_pos/RP_pos_10ppm_Scaled-feature-table-for-statistical-analysis.tsv", 
                                 header = TRUE, 
                                 check.names = FALSE)

#Data restructuring-------------------------------------------------------------

data_matrix_scaled_transposed <- t(data_matrix_scaled)

data_matrix <- matrix(as.numeric(data_matrix_scaled_transposed[-c(1:2),]),
                      ncol = 25)
rownames(data_matrix) <- rownames(data_matrix_scaled_transposed[3:2797,])

#Load meta data
pheno_data <- as.data.frame(data_matrix_scaled[,c(1:2)])
#pData <- read.table("data/RP_pos/RP_pos_10ppm_MetaData-28022023-11-36.tsv")
#pheno_data$replicate <- pData$repetition
#pheno_data$experiment <- pData$experiment
pheno_data$replicate <- c(rep(c(1,2,3),6), rep(c(1),7))
pheno_data$vial <- c(rep(c(1),3),rep(c(2),3),rep(c(3),3),rep(c(4),3),rep(c(5),3),rep(c(6),3), c(7, 8, 9, 10, 11, 12, 13))
pheno_data$sample_name <- paste(pheno_data$ATTRIBUTE_Sample,
                                pheno_data$replicate, 
                                pheno_data$vial, 
                                sep = "_")

colnames(data_matrix) <-  rownames(pheno_data)
rownames(pheno_data) == colnames(data_matrix) 


#Data exploration----------------------------------------------------------------
#Plot density and boxplot  
plot(density(data_matrix))

boxplot(data_matrix,
        las = 2,
        names = pheno_data$sample_name,
        ylim = c(-4,6))

#Plot heatmap of correlations
pheatmap(cor(data_matrix,method = "spearman"), 
         labels_col = pheno_data$sample_name, 
         labels_row = pheno_data$sample_name)

#Plot PCA
plotMDS(data_matrix, 
        labels = pheno_data$sample_name, 
        gene.selection = "common", 
        var.explained = TRUE)

mds <- plotMDS(data_matrix, 
               labels = pheno_data$sample_name, 
               gene.selection = "common", 
               var.explained = TRUE)

var_explained <- as.data.frame(mds$var.explained[1:7]*100)
colnames(var_explained) <- c("variance")

ggplot(var_explained, aes(x = rownames(var_explained), y = variance)) +
  geom_col() +
  geom_text(aes(label = round(variance, digits = 2)),vjust = -0.2) +
  xlab("principal component number") +
  ylab("variance (%)")

#Select non-pool samples -------------------------------------------------------
data_matrix <- data_matrix[,1:18]
pheno_data <- pheno_data[1:18,]

#Single-channel experimental design - Two Groups--------------------------------

des <- copy(pheno_data)
unique(des$ATTRIBUTE_Sample)
des$ATTRIBUTE_Sample <- factor(des$ATTRIBUTE_Sample,
                        levels = c("Ctrl", "Ibr"))

#the first coefficient estimates the mean intensity for control and plays the role
#of an intercept. The second coefficient estimates the difference between ibrutinib and control
des <- model.matrix(~ ATTRIBUTE_Sample, 
                    data = des)

# Model fit 
fit <- lmFit(data_matrix, des)
fit <- eBayes(fit)

#Extract results
coefs <- grep("ATTRIBUTE_Sample", colnames(coef(fit)), value = TRUE)
res <- data.table()

res <- rbind(res, data.table(
  topTable(fit, coef = coefs, adjust.method = "BH", number = nrow(data_matrix)),
  keep.rownames = TRUE,
  coef = gsub("ATTRIBUTE_Sample", "", coefs))
)

res[,direction := ifelse(logFC > 0, "up", "down")]

res[coef == "Ibr"][adj.P.Val < 0.05]
res[coef == "Ibr"][P.Value < 0.05]

# Summarize and Plot results 

# Number of significant hits
res[adj.P.Val < 0.05][,.N, by = c("coef", "direction")]
# Number of tested
res[,.N, by = c("coef", "direction")]

#Calculate row means of technical replicates------------------------------------------
data_matrix_scaled <- data_matrix_scaled[1:18, 3:2797]
data_matrix_scaled$experiment <- pheno_data$vial
#data_matrix_scaled <- data_matrix_scaled[,-c(1:3)]

data_matrix_scaled_means <- data_matrix_scaled %>%
  group_by(experiment) %>%
  summarise_all(mean)

data_matrix_scaled_means$ATTRIBUTE_Sample <- c(rep(c('Ctrl'),3), rep(c('Ibr'),3))

#Final datamatrix and pheno data
data_matrix_scaled_transposed <- t(data_matrix_scaled_means)
data_matrix <- matrix(as.numeric(data_matrix_scaled_transposed[2:2796,]), ncol = 6)


meta_data <- data.frame(vial = data_matrix_scaled_means$experiment,
                        ATTRIBUTE_Sample = data_matrix_scaled_means$ATTRIBUTE_Sample)

meta_data$sample_name <- paste(meta_data$ATTRIBUTE_Sample,
                               meta_data$vial, 
                               sep = "_")
rownames(data_matrix) <- rownames(data_matrix_scaled_transposed[2:2796,])
colnames(data_matrix) <- meta_data$sample_name
rownames(meta_data) <- meta_data$sample_name
rownames(meta_data) == colnames(data_matrix) 

pheno_data <- meta_data
#Wilcoxon signed rank test-----------------------------------------------------

?wilcox.test
res_wt <- 1

for (i in 1:nrow(data_matrix)) {
  res_test <- wilcox.test(data_matrix[i,1:3],
                          data_matrix[i,4:6], 
                          paired = TRUE)
  res_wt[i] <- res_test$p.value
}

#Compare p.values of limma stats and Wilcoxon signed rank tests
stats_res_wt <- data.frame('metabolites' = rownames(data_matrix),
                           'wilcox_pval' = res_wt)
stats_res_limma <- data.frame('metabolites' = res$rn,
                              'limma_pval' = res$P.Value)

stats_res <- stats_res_wt %>% inner_join(stats_res_limma,
                                         by = c('metabolites'))

ggplot(data = stats_res,mapping = aes(x = wilcox_pval, y = limma_pval)) +
  geom_point() +
  geom_smooth() 


































