# Create HRS survival data

library(readr)
library(dplyr)
library(corrplot)
library(glmnet)
library(preprocessCore)
library(FraminghamPACE)

Lehne_betas <- readRDS("../Lehne_Data/Lehne_Bvals_Whole_Sample.rds")
Lehne_pheno <- read_csv("../Lehne_Data/Clean_Pheno_FullDataset.csv")

# Ensure the subjects are rows
# And the row names are the subject ids
Lehne_betas <- as.data.frame(Lehne_betas)

head(rownames(Lehne_betas))
tail(rownames(Lehne_betas))
head(colnames(Lehne_betas))
tail(colnames(Lehne_betas))

sugden_CpGs <- rownames(Lehne_betas)
beta_IDs <- names(Lehne_betas)

betas <- Lehne_betas

#####################################
#####################################
# Predict FraminghamPACE

# Create a list of the dataframes with their desired column names
pace_list <- list(
  FraminghamPACE = getFraminghamPACE(betas, model_type = "FraminghamPACE"),
  Alternate1 = getFraminghamPACE(betas, model_type = "Alternate1"),
  Alternate2 = getFraminghamPACE(betas, model_type = "Alternate2"),
  Alternate3 = getFraminghamPACE(betas, model_type = "Alternate3")
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
write.csv(all_Framingham, "../Lehne_Data/Lehne_FraminghamPACE_df.csv")
