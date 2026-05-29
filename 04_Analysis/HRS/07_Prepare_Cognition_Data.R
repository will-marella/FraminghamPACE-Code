#### Prepare cognition data for analysis

library(haven)
library(tidyverse)
library(readr)

cog_outcome <- read_sas("../HRS_Data/Initial/cogfinalimp_9520wide.sas7bdat", NULL)


cog_outcome <- cog_outcome %>% mutate( hhidpn = paste0(HHID, PN), hhidpn = as.numeric(hhidpn)) %>% 
  select(hhidpn, cogfunction2016, cogfunction2018, cogfunction2020)

survival_data <- read_csv("../HRS_Data/HRS_complete_survival_data.csv")

# Define the fixed columns that always need to be included
# fixed_columns <- c("id", "IDATid", "Ethnicity", "Race", "Female", "PAGE", 
#                    "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos",
#                    "Mono", "Neu", "NK", "Treg",
#                    "bmi1213", "smoker1213", "pwavessmoked", "n_smokedata")

# Automatically find all clock-related columns
clock_patterns <- c("^ElasticNet", "^Ridge",
                    "^DunedinPACE", "^PCGrimAge", "^PCPhenoAge")
clock_columns <- grep(paste(clock_patterns, collapse = "|"), names(survival_data), value = TRUE)

# Combine and select
# pheno <- survival_data %>%
#   dplyr::select(all_of(c(fixed_columns, clock_columns)))

cognition_data <- merge(survival_data, cog_outcome, by.x="id", by.y="hhidpn", all=FALSE)

summary(cognition_data$cogfunction2020)
names(cognition_data)

write.csv(cognition_data, "../HRS_Data/HRS_cognition_data.csv")


