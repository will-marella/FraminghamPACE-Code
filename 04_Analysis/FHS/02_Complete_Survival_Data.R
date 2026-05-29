# Load packages
library(readr)
library(dplyr)
library(haven)

################################
# Load and merge data
################################

# Read in data
outcome_data <- read_csv("../FHS_Data/adjusted_filtered_outcomes.csv")
pheno <- read_csv("../FHS_Data/Initial/framingham_offspring_clean_pheno_dataset.csv")

pheno <- pheno %>%
  dplyr::select("dbGaP_Subject_ID", "Barcode",
         "Age", "Female")


# Add BMI
biomarkers <- read_csv("../FHS_Data/Initial/biomarkers_processed_Furuya.csv")

bmi_subset <- biomarkers %>%
  dplyr::select(c("subject_id", "visit", "bmi_ln_adj_sc"))

bmi_subset_visit8 <- bmi_subset[bmi_subset$visit==8,]

bmi_subset_visit8 <- bmi_subset_visit8 %>%
  dplyr::select(c("subject_id", "bmi_ln_adj_sc"))

pheno <- merge(pheno, bmi_subset_visit8, by.x="dbGaP_Subject_ID", "subject_id")


# Merge
survival_data <- merge(pheno, outcome_data,
                                by.x = "dbGaP_Subject_ID", by.y = "subject_id",
                                all=FALSE)

# Add clocks
FraminghamPACE_df <- read_csv("../FHS_Data/FraminghamPACE_WTM_noPoA.csv")
PC_Clocks <- read_csv("../FHS_Data/Initial/Framingham_PC_Clocks_CPR_Noob_normalized.csv")
PC_Clocks <- PC_Clocks %>%
  dplyr::select(-c("Age", "Female"))
DunedinPACE <- read_csv("../FHS_Data/Initial/Framingham_DunedinPACE.csv")
# FraminghamPoA <- read_csv("../FHS_Data/Initial/FraminghamPoA.csv")
WTM_FraminghamPoA <- read_csv("../FHS_Data/Initial/WTM_FraminghamPoA.csv")
# WTM_FraminghamPoA <- WTM_FraminghamPoA %>%
#   rename(WTM_FraminghamPoA = FraminghamPoA)

complete_survival_data <- merge(survival_data, FraminghamPACE_df, by.x="dbGaP_Subject_ID", by.y="BARCODE")
complete_survival_data <- merge(complete_survival_data, PC_Clocks, by="Barcode")
complete_survival_data <- merge(complete_survival_data, DunedinPACE, by="Barcode")
# complete_survival_data <- merge(complete_survival_data, FraminghamPoA, by.x="dbGaP_Subject_ID", by.y="subject_id")
complete_survival_data <- merge(complete_survival_data, WTM_FraminghamPoA, by="dbGaP_Subject_ID", by.y="subject_id")

# Add cells

cells <- read_csv("../FHS_Data/Initial/Framingham_CPR_SALAS_cell_counts.csv")
complete_survival_data <- merge(complete_survival_data, cells, by="Barcode")


################################
# Define clock variables
################################

# Define all clock variables that need residuals
framingham_vars <- c("Ridge", "ElasticNet")
special_vars <- c("FraminghamPoA", "DunedinPACE")  # Special FraminghamPoA variable
already_resid_vars <- c("PCGrimAgeResid", "PCPhenoAgeResid")  # Already have Resid versions
non_resid_vars <- c("PCGrimAge", "PCPhenoAge")  # Need scaling but no residuals

# All variables that need residuals
vars_needing_residuals <- c(framingham_vars, special_vars)

################################
# Add residuals and scaling automatically
################################

# Create age residuals for variables that need them
resid_names <- paste0(vars_needing_residuals, "Resid")

complete_survival_data <- complete_survival_data %>%
  dplyr::select(-starts_with("..."))
anyDuplicated(names(complete_survival_data))

complete_survival_data <- complete_survival_data %>%
  mutate(across(all_of(vars_needing_residuals), 
                ~ resid(lm(.x ~ Age, data = complete_survival_data)), 
                .names = "{.col}Resid"))

complete_survival_data <- complete_survival_data %>%
  # Create all residuals at once
  mutate(across(all_of(vars_needing_residuals), 
                ~ resid(lm(.x ~ Age, data = complete_survival_data)), 
                .names = "{.col}Resid")) %>%
  # Scale all variables at once
  mutate(
    # Scale original variables that don't get residuals
    across(all_of(c(special_vars, non_resid_vars, already_resid_vars)), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled"),
    # Scale all residual variables  
    across(all_of(resid_names), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled"),
  )

################################
# Clean up and save
################################

# Clean up unnecessary columns
clean_survival_data <- complete_survival_data %>%
  dplyr::select(-matches("^\\.\\.\\.[0-9]+"))  # Remove columns like ...1.x, ...1.y, etc.

# Check final column names
cat("Final columns:\n")
print(names(clean_survival_data))

summary(clean_survival_data)

# Write complete survival data
write.csv(clean_survival_data, "../FHS_Data/complete_survival_data.csv", row.names = FALSE)
