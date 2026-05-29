### Load libraries
library(dplyr)
library(readr)

### Load probe list
probe_list <- read.csv("../../FittingPoA/Data/Initial/probelist_450k_epicv1_epicv2.csv")
probe_list <- probe_list$probe

# Load in beta matrix
beta_matrix <- readRDS("../WHI_Data/Initial/GRSet_fully_filtered_bmiq.rds")
whi_pheno_backbone <- read_csv("../WHI_Data/Initial/whi_clean_pheno_dataset.csv")

# Ensure they match the pheno barcodes
all_grset_barcodes <- names(beta_matrix)
all_pheno_barcodes <- whi_pheno_backbone$Barcode

print(length(intersect(all_grset_barcodes, 
                       all_pheno_barcodes)) == length(all_grset_barcodes)) # TRUE

# Filter to 450k, EPICv1, EPICv2, subset
reduced_beta_matrix <- beta_matrix[rownames(beta_matrix) %in% probe_list,]

# Write to CSV
write.csv(reduced_beta_matrix, "../WHI_Data/WHI_betas_450k_EPICv1_EPICv2.csv", row.names = TRUE)


