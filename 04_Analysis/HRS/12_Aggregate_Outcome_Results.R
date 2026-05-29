### Load packages
library(readr)
library(dplyr)

### Load datasets
# Survival
survival_minimal <- read_csv("../HRS_Output/HRS_Survival_Results_Minimal.csv")
survival_cells <- read_csv("../HRS_Output/HRS_Survival_Results_Cells.csv")
survival_bmi <- read_csv("../HRS_Output/HRS_Survival_Results_BMI.csv")
survival_smoking <- read_csv("../HRS_Output/HRS_Survival_Results_Smoking.csv")

# Disease/Disability
disease_disability_minimal <- read_csv("../HRS_Output/HRS_Disease_Disability_Results_Minimal.csv")
disease_disability_cells <- read_csv("../HRS_Output/HRS_Disease_Disability_Results_Cells.csv")
disease_disability_bmi <- read_csv("../HRS_Output/HRS_Disease_Disability_Results_BMI.csv")
disease_disability_smoking <- read_csv("../HRS_Output/HRS_Disease_Disability_Results_Smoking.csv")

# Cognition
cognition_minimal <- read_csv("../HRS_Output/HRS_Cognition_Results_Minimal.csv")
cognition_cells <- read_csv("../HRS_Output/HRS_Cognition_Results_Cells.csv")
cognition_bmi <- read_csv("../HRS_Output/HRS_Cognition_Results_BMI.csv")
cognition_smoking <- read_csv("../HRS_Output/HRS_Cognition_Results_Smoking.csv")


### Adjust datasets for merge

# Remove 'Dataset' column from the survival datasets
survival_minimal$Dataset <- NULL
survival_cells$Dataset <- NULL
survival_bmi$Dataset <- NULL
survival_smoking$Dataset <- NULL

# Remove the `...1` column from each dataset
survival_minimal$`...1` <- NULL
survival_cells$`...1` <- NULL
survival_bmi$`...1` <- NULL
survival_smoking$`...1` <- NULL

disease_disability_minimal$`...1` <- NULL
disease_disability_cells$`...1` <- NULL
disease_disability_bmi$`...1` <- NULL
disease_disability_smoking$`...1` <- NULL

cognition_minimal$`...1` <- NULL
cognition_cells$`...1` <- NULL
cognition_bmi$`...1` <- NULL
cognition_smoking$`...1` <- NULL

# Need to adjust some more column names
disease_disability_minimal <- disease_disability_minimal %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

disease_disability_cells <- disease_disability_cells %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

disease_disability_bmi <- disease_disability_bmi %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

disease_disability_smoking <- disease_disability_smoking %>%
  rename(outcome = Outcome,
         clock = Clock,
         hazard_ratio = rate_ratio)

# Merge datasets
minimal_results <- rbind(survival_minimal, disease_disability_minimal, cognition_minimal)
cells_results <- rbind(survival_cells, disease_disability_cells, cognition_cells)
bmi_results <- rbind(survival_bmi, disease_disability_bmi, cognition_bmi)
smoking_results <- rbind(survival_smoking, disease_disability_smoking, cognition_smoking)

### Save results
write.csv(minimal_results, "../HRS_Output/HRS_Combined_Outcome_Results_Minimal.csv")
write.csv(cells_results, "../HRS_Output/HRS_Combined_Outcome_Results_Cells.csv")
write.csv(bmi_results, "../HRS_Output/HRS_Combined_Outcome_Results_BMI.csv")
write.csv(smoking_results, "../HRS_Output/HRS_Combined_Outcome_Results_Smoking.csv")

