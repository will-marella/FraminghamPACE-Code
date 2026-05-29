# Create HRS survival data

library(readr)
library(dplyr)
library(corrplot)
library(glmnet)
library(preprocessCore)
library(FraminghamPACE)

HRS_raw_betas <- read_csv("../HRS_Data/HRS_betas_450k_epicv1_epicv2.csv")
hrs_outcomes <- read_csv("../HRS_Data/HRS_outcomes_and_clocks.csv")

# Ensure the subjects are rows
# And the row names are the subject ids
HRS_raw_betas <- as.data.frame(HRS_raw_betas)

samples <- HRS_raw_betas$`...1`
rownames(HRS_raw_betas) <- HRS_raw_betas$`...1`

HRS_raw_betas <- HRS_raw_betas %>%
  dplyr::select(-`...1`)

sugden_CpGs <- names(HRS_raw_betas)
beta_IDs <- rownames(HRS_raw_betas)

betas <- t(HRS_raw_betas)

head(colnames(betas))
head(rownames(betas))

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

########################################################
########################################################
# Export

write.csv(all_Framingham, "../HRS_Data/HRS_FraminghamPACE_df.csv")

