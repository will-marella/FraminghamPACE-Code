# estimate_FraminghamPACE.R

# Load necessary packages
library(readr)
library(glmnet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(haven)
library(limma)
library(devtools)
library(RColorBrewer)
library(preprocessCore)
library(FraminghamPACE)

#####################################
## Prepare betas
#####################################

CALERIE_raw_betas <- read_csv("../CALERIE_Data/CALERIE_probefiltered.csv")

CALERIE_raw_betas <- as.data.frame(CALERIE_raw_betas)
rownames(CALERIE_raw_betas) <- CALERIE_raw_betas$...1
beta_barcodes <- CALERIE_raw_betas$`...1`
CALERIE_raw_betas <- CALERIE_raw_betas %>%
  dplyr::select(-`...1`)

head(names(CALERIE_raw_betas))
head(rownames(CALERIE_raw_betas))
tail(names(CALERIE_raw_betas))
tail(rownames(CALERIE_raw_betas))

betas <- t(CALERIE_raw_betas)
betas <- as.data.frame(betas)

#####################################
#####################################
# Predict FraminghamPACE

# Create a list of the dataframes with their desired column names
pace_list <- list(
  ElasticNet = getFraminghamPACE(betas, model_type = "ElasticNet"),
  Ridge = getFraminghamPACE(betas, model_type = "Ridge")
)

# Unique to CALERIE!
# Find common samples across all PACE results
common_samples <- Reduce(intersect, lapply(pace_list, rownames))
cat("Common samples across all clocks:", length(common_samples), "\n")

# Subset all dataframes to common samples and combine
all_Framingham <- do.call(cbind, lapply(names(pace_list), function(name) {
  df <- pace_list[[name]][common_samples, , drop = FALSE]  # Subset to common samples
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

write.csv(all_Framingham, "../CALERIE_Data/CALERIE_FraminghamPACE_df.csv")

