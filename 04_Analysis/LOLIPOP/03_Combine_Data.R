# Load libraries
library(readr)
library(dplyr)

# Load data
FraminghamPACE_df <- read_csv("../Lehne_Data/Lehne_combined_FraminghamPACE_replicates.csv")
Clock_df <- read_csv("../Lehne_Data/Lehne_clocks_and_PC_Clocks_replicates.csv")

# Subset data
clean_Framingham <- FraminghamPACE_df %>%
  dplyr::select(-"...1")

clean_Clocks <- Clock_df %>%
  dplyr::select(c("SampleID", "sample",
                  "Age", "Female",
                  "CD8T", "CD4T",
                  "NK", "Bcell",
                  "Mono", "Gran",
                  "DunedinPoAm_45",
                  "PCPhenoAge", "PCGrimAge",
                  "PCPhenoAgeResid", "PCGrimAgeResid")) %>%
  dplyr::rename(subject_id = sample,
                BARCODE = SampleID,
                DunedinPACE = DunedinPoAm_45)

# Merge data
combined_data <- merge(clean_Clocks, clean_Framingham, by="BARCODE")

########################################################
########################################################
# Residualization


# Define variable groups - exclude already scaled versions
framingham_vars <- c("ElasticNet",
                     "Ridge")

other_pace_vars <- c("DunedinPACE")
already_resid_vars <- c("PCGrimAgeResid", "PCPhenoAgeResid")  # Already have Resid versions
non_resid_vars <- c("PCGrimAge", "PCPhenoAge")  # Need scaling but no residuals

# Create age residuals for FraminghamPACE and other PACE variables
resid_vars <- c(framingham_vars, other_pace_vars)
resid_names <- paste0(resid_vars, "Resid")

combined_data <- combined_data %>%
  # Create all residuals at once
  mutate(across(all_of(resid_vars), 
                ~ resid(lm(.x ~ Age, data = combined_data)), 
                .names = "{.col}Resid")) %>%
  # Scale all variables at once
  mutate(
    # Scale original variables that don't get residuals
    across(all_of(c(other_pace_vars, non_resid_vars, already_resid_vars)), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled"),
    # Scale all residual variables
    across(all_of(resid_names), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled")
  )

########################################################
########################################################
# Export

write.csv(combined_data, "../Lehne_Data/Lehne_Replicate_Combined_Data.csv")
