# Predict FraminghamPACEa and b

library(readr)
library(dplyr)
library(FraminghamPACE)

## Read in sugden probes, phenotype info
lehne_sugdens <- read_csv("../Lehne_Data/Lehne_SugdenProbes_replicates.csv")

## Wrangle betas and age_sex_df into correct format

# betas
lehne_sugdens <- as.data.frame(lehne_sugdens)
rownames(lehne_sugdens) <- lehne_sugdens$`...1`
lehne_sugdens <- lehne_sugdens %>%
  dplyr::select(-`...1`)

lehne_sugdens <- t(lehne_sugdens)
betas <- as.matrix(lehne_sugdens)

#####################################
#####################################
# Predict FraminghamPACE

# Create a list of the dataframes with their desired column names
pace_list <- list(
  ElasticNet = getFraminghamPACE(betas, model_type = "ElasticNet"),
  Ridge = getFraminghamPACE(betas, model_type = "Ridge")
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
# Merge and write

write.csv(all_Framingham, "../Lehne_Data/Lehne_combined_FraminghamPACE_replicates.csv")
