# Predict FraminghamPACE in WHI

# Load necessary packages
library(readr)
library(glmnet)
library(dplyr)
library(corrplot)
library(preprocessCore)
library(FraminghamPACE)

#####################################
#####################################
## Set up Data

# Load raw beta values and outcomes
WHI_raw_betas <- readRDS("../WHI_Data/Initial/GRSet_fully_filtered_bmiq.rds")
whi_outcomes <- read_csv("../WHI_Data/WHI_complete_pheno_outcome_dataset.csv")

# Switch to wide format
WHI_raw_betas <- as.data.frame(WHI_raw_betas)

# Create betas
# betas <- t(WHI_raw_betas)
betas <- WHI_raw_betas

head(colnames(betas))
tail(colnames(betas))
head(rownames(betas))
tail(rownames(betas))


#####################################
#####################################
# Predict FraminghamPACE

# Create a list of the dataframes with their desired column names
pace_list <- list(
  Ridge = getFraminghamPACE(betas, model_type = "Ridge"),
  ElasticNet = getFraminghamPACE(betas, model_type = "ElasticNet")
)

# Combine into a single dataframe with proper column names
all_Framingham <- do.call(cbind, lapply(names(pace_list), function(name) {
  df <- pace_list[[name]]
  colnames(df) <- name  # Rename the column to match the list name
  return(df)
}))

# Convert to regular dataframe (in case it's a matrix)
all_Framingham <- as.data.frame(all_Framingham)
all_Framingham$BARCODE <- rownames(all_Framingham)

# Automatically scale all PACE columns
pace_columns <- names(pace_list)  # Get all PACE column names
scaled_columns <- paste0(pace_columns, "_scaled")  # Create scaled column names

# Scale all at once
all_Framingham[scaled_columns] <- lapply(all_Framingham[pace_columns], scale)

#####################################
#####################################
## Export Framingham df
write.csv(all_Framingham, "../WHI_Data/WHI_FraminghamPACE_df.csv")

