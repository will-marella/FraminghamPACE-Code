##############################################
# Predict FraminghamPACE from CARDIA
##############################################

##############################################
# Load libraries, sugden betas, phenotype info
##############################################

library(readr)
library(dplyr)
library(FraminghamPACE)

cardia_sugdens <- read_csv("../CARDIA_Data/CARDIA_SugdenProbes.csv")

##############################################
# Wrangle betas and age_sex_df into correct format
##############################################

# betas
cardia_sugdens <- as.data.frame(cardia_sugdens)
rownames(cardia_sugdens) <- cardia_sugdens$`...1`
cardia_sugdens <- cardia_sugdens %>%
  dplyr::select(-`...1`)

cardia_sugdens <- t(cardia_sugdens)
betas <- as.matrix(cardia_sugdens)

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

########################################################
# Export
########################################################

write.csv(all_Framingham, "../CARDIA_Data/CARDIA_combined_FraminghamPACE.csv")
