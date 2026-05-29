# predict in-sample FraminghamPACEa and b

# Load Packages
library(glmnet)
library(caret)
library(readr)
library(dplyr)
library(FraminghamPACE)

# Load in data
data <- read_csv("../FHS_Data/FraminghamPoA_Training_Set_v2_TruD.csv")
data <- as.data.frame(data)
rownames(data) <- data$...1
data <- data %>%
  dplyr::select(-...1)

########################################
# Clean FraminghamPoA
########################################

# Force sex-specific means to 1
data$FraminghamPoA <- data$FraminghamPoA + (1 - ave(data$FraminghamPoA, data$Female))

print("Sex-specific means after adjustment:")
print(tapply(data$FraminghamPoA, data$Female, mean))

subject_ids <- rownames(data)
PoA_vals <- data$FraminghamPoA

RawFraminghamPoA_df <- data.frame(
  subject_id = subject_ids,
  RawFraminghamPoA = PoA_vals
)


########################################
# Predict FraminghamPACEa and b
########################################

betas_t <- data %>%
  dplyr::select(-c(FraminghamPoA, Age, Female))

betas <- t(betas_t)

############

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

########################################
# Export data
########################################

write.csv(all_Framingham, "../FHS_Data/FraminghamPACE_WTM_noPoA.csv")
