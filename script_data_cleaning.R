
## Dataset load and explore-------------
load("data/dataset.RDS")


## Blanks removal-----------------
# When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
cutoff = 0.3 # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3

# Getting mean for every feature in blank and Samples # set skipna = False do not exclude NA/null values when computing the result.
blank_avg <- rowMeans(data_matrix[ , pData$type == "blank"],na.rm = FALSE)
sample_avg <- rowMeans(data_matrix[ , pData$type == "sample"],na.rm=FALSE)

# Getting the ratio of blank vs samples
ratio_blank_samples <- (blank_avg+1)/(sample_avg+1)

# Create an array with boolean values: True (is a real feature, ratiocutoff)
is_real_feature <- ratio_blank_samples<cutoff

# Checking if there are any NA values present. Having NA values in the 4 variables will affect the final dataset to be created
print(sum(is.nan(blank_avg)))
print(sum(is.nan(sample_avg)))
print(sum(is.nan(ratio_blank_samples)))
print(sum(is.nan(is_real_feature)))

# Calculating the number of background features and features present 
no_real_features <- sum(is_real_feature == TRUE)
print(no_real_features)
no_background_features <- sum(is_real_feature == FALSE)
print(no_background_features)


#Remove blank features from the dataset
blank_removal <- data_matrix[is_real_feature,]
dim(blank_removal)

#Remove blank samples from the dataset and metadata
data_matrix_blank_removal <- blank_removal[,10:34]
dim(data_matrix_blank_removal)
pData_blank_removal <- pData[10:34,]
dim(pData_blank_removal)

rownames(pData_blank_removal) == colnames(data_matrix_blank_removal)


## Imputation of zero values -----------------------
print(sum((is.nan(data_matrix_blank_removal)==FALSE)))
print(sum((is.nan(data_matrix_blank_removal)==TRUE)))
sum(data_matrix_blank_removal == 0)


## Normalization---------------------
#The following code performs sample-centric (column-wise) normalisation
# Dividing each element of a particular column with its column sum

## Transposing------------------

## Scaling------------------

## Statistics (e.g. limma)-------------------

