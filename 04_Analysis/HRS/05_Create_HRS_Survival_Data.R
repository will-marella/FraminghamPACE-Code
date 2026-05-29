# Create HRS survival data

library(readr)
library(glmnet)
library(dplyr)
library(preprocessCore)
library(FraminghamPACE)

################################
# Load and merge data
################################

hrs_outcomes <- read_csv("../HRS_Data/HRS_outcomes_and_clocks.csv")
FraminghamPACE_df <- read_csv("../HRS_Data/HRS_FraminghamPACE_df.csv")

# Add predictions
pheno_clocks <- merge(hrs_outcomes, FraminghamPACE_df, by.x="IDATid", by.y="...1")

# Add cells
cells <- read_csv("../HRS_Data/Initial/HRS_combined_Salas_cells.csv")
pheno_clocks_cells <- merge(pheno_clocks, cells, by="IDATid")

# Add additional pheno
additional_pheno <- read_csv("../HRS_Data/HRS_clean_additional_pheno.csv")
survival_data <- merge(pheno_clocks_cells, additional_pheno, by="IDATid", all=FALSE)

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
  # Create all residuals at once
  mutate(across(all_of(vars_needing_residuals), 
                ~ resid(lm(.x ~ PAGE, data = survival_data)), 
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
  select(-matches("^\\.\\.\\.[0-9]+"))  # Remove columns like ...1.x, ...1.y, etc.

# Some HRS-Specific stuff
clean_survival_data <- clean_survival_data %>%
  rename(hhid = hhid.x,
         pn = pn.x) %>%
  select(-c("hhid.y", "pn.y"))

# Check final column names
cat("Final columns:\n")
print(names(clean_survival_data))

### Write to CSV
write.csv(clean_survival_data, "../HRS_Data/HRS_complete_survival_data.csv")
