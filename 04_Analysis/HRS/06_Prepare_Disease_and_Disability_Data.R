# Prepare the additional HRS outcomes for further analysis

library(haven)
library(readr)
library(tidyverse)
library(dplyr)

outcome_2020 <- read_dta("../HRS_Data/Initial/HRSRAND_Health2020_DWB.dta")
survival_data <- read_csv("../HRS_Data/HRS_complete_survival_data.csv")

# # Define the fixed columns that always need to be included
# fixed_columns <- c("id", "IDATid", "Ethnicity", "Race", "Female", "PAGE",
#                    "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos",
#                    "Mono", "Neu", "NK", "Treg",
#                    "bmi1213", "smoker1213", "pwavessmoked", "n_smokedata")

# Automatically find all clock-related columns
clock_patterns <- c("^Ridge", "^ElasticNet",
                    "^DunedinPACE", "^PCGrimAge", "^PCPhenoAge")
clock_columns <- grep(paste(clock_patterns, collapse = "|"), names(survival_data), value = TRUE)

# Combine and select
# pheno <- survival_data %>%
#   dplyr::select(all_of(c(fixed_columns, clock_columns)))

chrondxe_all <- outcome_2020 %>%
  dplyr::select(hhidpn, raedyrs, ends_with(c("13", "14","15"))) %>% 
  dplyr::select(hhidpn, raedyrs, starts_with("chrondxe"), -contains("w"))

ADLs_all <- outcome_2020 %>% 
  dplyr::select(hhidpn, raedyrs, ends_with(c("l5a13", "l5a14","l5a15"))) %>% 
  dplyr::select(hhidpn, raedyrs, starts_with("adl"))

iADLs_all <- outcome_2020 %>% 
  dplyr::select(hhidpn, raedyrs, ends_with(c("l5a13", "l5a14","l5a15"))) %>% 
  dplyr::select(hhidpn, raedyrs, starts_with("iadl"))

disease_and_disability <- survival_data %>%
  dplyr::inner_join(chrondxe_all, by = c("id" = "hhidpn")) %>%
  dplyr::inner_join(ADLs_all, by = c("id" = "hhidpn")) %>%
  dplyr::inner_join(iADLs_all, by = c("id" = "hhidpn"))


names(disease_and_disability)

write.csv(disease_and_disability, "../HRS_Data/HRS_disease_disability_data.csv")


