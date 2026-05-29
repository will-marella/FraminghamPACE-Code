# Load packages
library(readr)
library(dplyr)
library(haven)

################################
# Load and merge data
################################

# FraminghamPACE clocks
FraminghamPACE <-read_csv("../CALERIE_Data/CALERIE_FraminghamPACE_df.csv")
FraminghamPACE <- FraminghamPACE %>%
  dplyr::select(-`...1`)
FraminghamPACE$BARCODE <- as.character(FraminghamPACE$BARCODE)

# Load DunedinPACE
DunedinPACE <- read_dta("../CALERIE_Data/Initial/CALERIE_ClocksPoAmPACE.dta")
DunedinPACE <- DunedinPACE %>%
  dplyr::select(barcode, pace) %>%
  rename(DunedinPACE = pace,
         BARCODE = barcode)
DunedinPACE$BARCODE <- as.character(DunedinPACE$BARCODE)

# Load PC Clocks
PC_Clocks <- read_csv("../CALERIE_Data/Initial/CALERIE_PC-clocks_agedwb_CPR_fixed_no_deid_or_age.csv")
PC_Clocks <- PC_Clocks %>%
  dplyr::select(barcode, 
         PCGrimAge, PCPhenoAge,
         PCGrimAgeResid, PCPhenoAgeResid) %>%
  rename(BARCODE = barcode)

# Merge all with phenotype info
pheno_raw <- read_dta("../CALERIE_Data/Initial/CALERIE_Backbone.dta")
pheno_clean <- pheno_raw %>%
  dplyr::select(barcode, deidnum, fu, CR, agevis, female, race, ethnic, bmistrat, deidsite) %>%
  rename(BARCODE = barcode)

PACE_clocks <- merge(FraminghamPACE, DunedinPACE, by="BARCODE")
all_clocks <- merge(PACE_clocks, PC_Clocks, by="BARCODE")

# Merge added vars with pheno
added_vars <- haven::read_dta("../CALERIE_Data/Initial/CALERIENatAgingPhenoData.dta")

clean_added_vars <- added_vars %>%
  dplyr::select(deidnum, fu, race3, agedwb)

pheno <- pheno_clean %>% 
  left_join(clean_added_vars, by = c("deidnum", "fu"))

survival_data_raw <- merge(pheno, all_clocks, by="BARCODE")

## Add cell counts
cells <- read_csv("../CALERIE_Data/Initial/CALERIE_CPR_EPIC_cell_counts.csv")

cells <- cells %>%
  rename(BARCODE = barcode)

survival_data <- merge(survival_data_raw, cells, by="BARCODE")

################################
# Define clock variables
################################

# Define all clock variables that need residuals
framingham_vars <- c("ElasticNet", "Ridge")
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
                ~ resid(lm(.x ~ agevis, data = survival_data)), 
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
survival_data <- survival_data %>%
  dplyr::select(-matches("^\\.\\.\\.[0-9]+"))  # Remove columns like ...1.x, ...1.y, etc.

# Check final column names
cat("Final columns:\n")
print(names(survival_data))
  
## Write survival data
write.csv(survival_data, "../CALERIE_Data/CALERIE_Combined_Data.csv")
