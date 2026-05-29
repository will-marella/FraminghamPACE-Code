# Load packages
library(readr)
library(dplyr)

################################
# Load and merge data
################################

# Load pheno/outcome data
whi_outcomes <- read_csv("../WHI_Data/WHI_complete_pheno_outcome_dataset.csv")

# Load clock data
FraminghamPACE_df <- read_csv("../WHI_Data/WHI_FraminghamPACE_df.csv")
Dunedin_data <- read_csv("../WHI_Data/Initial/WHI_raw_pace_w_attributes.csv")
PC_Clock_Data <- read_csv("../WHI_Data/Initial/WHI_PC_Clocks_Raw_Betas_Assay_Info_CPR.csv")

# Clean Dunedin/PC Clock Data

# Create BARCODE for merge
Dunedin_data$BARCODE <- paste(Dunedin_data$methyl_array, Dunedin_data$array_pos, sep="_")
PC_Clock_Data$BARCODE <- paste(PC_Clock_Data$methyl_array, PC_Clock_Data$array_pos, sep="_")

# Select important columns
Dunedin_data <- Dunedin_data %>%
  dplyr::select(BARCODE, DunedinPACE)

PC_Clock_Data <- PC_Clock_Data %>%
  dplyr::select(BARCODE, PCGrimAgeResid, PCPhenoAgeResid, PCGrimAge, PCPhenoAge)

# Merge
Framingham_outcomes <- merge(FraminghamPACE_df, whi_outcomes, by="BARCODE", all=FALSE)
Framingham_Dunedin_outcomes <- merge(Framingham_outcomes, Dunedin_data, by='BARCODE')

Framingham_clocks <- merge(Framingham_Dunedin_outcomes, PC_Clock_Data, by='BARCODE')

# Add in the additional phenotype information
additional_pheno <- read_csv("../WHI_Data/WHI_additional_pheno.csv")

survival_data <- merge(Framingham_clocks, additional_pheno, by="BARCODE", all=FALSE)

cells <- read_csv("../WHI_Data/Initial/WHI_CPR_SALAS_cell_counts.csv")

# Add them all for simplicity
survival_data <- merge(survival_data, cells, by.x="BARCODE", by.y="Barcode", all=FALSE)

survival_data$RACE <- as.factor(survival_data$RACE)

# Handle the residualization
# survival_data <- survival_data[complete.cases(survival_data),]
sum(is.na(survival_data$DunedinPACE)) # 28


################################
# Define clock variables
################################

# Define all clock variables that need residuals
framingham_vars <- c("Ridge", "ElasticNet")
special_vars <- c("DunedinPACE")  # Special FraminghamPoA variable
already_resid_vars <- c("PCGrimAgeResid", "PCPhenoAgeResid")  # Already have Resid versions
non_resid_vars <- c("PCGrimAge", "PCPhenoAge")  # Need scaling but no residuals

# All variables that need residuals
vars_needing_residuals <- c(framingham_vars, special_vars)

################################
# Add residuals and scaling automatically
################################

# Create age residuals for variables that need them
resid_names <- paste0(vars_needing_residuals, "Resid")

survival_data <- survival_data %>%
  # Create all residuals at once, handling NAs properly
  mutate(across(all_of(vars_needing_residuals), 
                ~ {
                  model <- lm(.x ~ AGE, data = survival_data, na.action = na.exclude)
                  resid(model)
                }, 
                .names = "{.col}Resid")) %>%
  # Scale all variables at once (scale() handles NAs automatically)
  mutate(
    # Scale original variables that don't get residuals
    across(all_of(c(special_vars, non_resid_vars, already_resid_vars)), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled"),
    # Scale all residual variables  
    across(all_of(resid_names), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled"),
    # Scale the framingham variables (if they need original scaled versions)
    across(all_of(framingham_vars), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled")
  )

################################
# Clean up and save
################################

# Clean up unnecessary columns
clean_survival_data <- survival_data %>%
  dplyr::select(-matches("^\\.\\.\\.[0-9]+"))  # Remove columns like ...1.x, ...1.y, etc.

# Check final column names
cat("Final columns:\n")
print(names(clean_survival_data))

#### Write data to csv
write.csv(clean_survival_data, "../WHI_Data/WHI_survival_data.csv")
