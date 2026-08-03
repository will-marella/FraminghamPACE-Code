library(haven)
library(readr)
library(dplyr)

# Read in pheno data, outcomes
pheno_data <- read_csv("../WHI_Data/Initial/whi_clean_pheno_dataset.csv")
outcomes <- read.csv("../WHI_Data/Initial/WHI_pheno_outcomes.csv")

# Remove rows with samples that did not pass our QC
pheno_data_passedQC <- pheno_data[pheno_data$passed_cu_core_QC == "pass", ]


# Merge dataframes
pheno_and_outcomes_messy <- merge(pheno_data_passedQC, outcomes,
                                  by.x = "SUBJECT_ID",
                                  by.y = "SUBJID",
                                  all=FALSE)

# Select only relevant columns
# Will match with the probes via Barcode
relevant_columns <- c("Barcode", "dbgap_subject_id",
                      "age", "race",
                      "CHD", "CHDDY",
                      "CHF", "CHFDY",
                      "STROKE", "STROKEDY",
                      "ANYCANCER", "ANYCANCERDY",
                      "ANYFX", "ANYFXDY",
                      "DEATH", "DEATHDY",
                      "ENDFOLLOWDY")

names(pheno_and_outcomes_messy)

# Keep only the relevant columns
# Rename so all columns are capitalized, for consistency
pheno_and_outcomes <- pheno_and_outcomes_messy %>%
  dplyr::select(all_of(relevant_columns)) %>%
  rename(
    BARCODE = Barcode,
    AGE = age,
    RACE = race
  )

# Save complete pheno/outcomes dataset
write.csv(pheno_and_outcomes, "../WHI_Data/WHI_complete_pheno_outcome_dataset.csv")