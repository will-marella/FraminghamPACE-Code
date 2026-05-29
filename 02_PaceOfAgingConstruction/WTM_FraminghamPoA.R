## Load packages  --------------------------------------------

library(dplyr)
library(eha)

## Load Data  --------------------------------------------

RawPoA <- read.csv("../Data/WTM_PoA.csv")
backbone <- read.csv("../Data/WTM_backbone.csv")
outcomes <- read.csv("../Data/outcomes.csv")
biomarkers_final <- read.csv("../Data/WTM_biomarkers_final.csv")

## Create Survival Dataset  -----------------------------

age_sex_df <- backbone %>%
  select(subject_id, baseline_age, sex) %>%
  distinct()

RawPoA_age_sex <- merge(RawPoA, age_sex_df, by="subject_id")
complete_survival_data <- merge(RawPoA_age_sex, outcomes, by="subject_id")

## Get Participants with 13+ repeated measurements  -----------------------------

# Count biomarkers per subject that have 2+ measurements
biomarker_counts <- biomarkers_final %>%
  group_by(subject_id, biomarker) %>%
  summarise(n_measurements = n(), .groups = "drop") %>%
  filter(n_measurements >= 2) %>%
  group_by(subject_id) %>%
  summarise(n_biomarkers_with_repeats = n(), .groups = "drop") %>%
  filter(n_biomarkers_with_repeats >= 13)

# Get the list of qualifying subject_ids
qualifying_subjects <- biomarker_counts$subject_id

# Filter your complete_survival_data
complete_survival_data <- complete_survival_data %>%
  filter(subject_id %in% qualifying_subjects)


## Create Sex-Stratified Outcome Datasets  -----------------------------

male_survival_data <- complete_survival_data[complete_survival_data$sex == "Male",]
female_survival_data <- complete_survival_data[complete_survival_data$sex == "Female",]

## Create Slope Ratios  -----------------------------

male_survival_data$slope_ratio <- male_survival_data$quadHM_slopes / mean(male_survival_data$quadHM_slopes, na.rm = TRUE)
female_survival_data$slope_ratio <- female_survival_data$quadHM_slopes / mean(female_survival_data$quadHM_slopes, na.rm = TRUE)

## Fit Gompertz Models to get log hazards  -----------------------------

# 1. Male age model (to get age_hazard for males)
male_age_gompertz <- phreg(Surv(death_date, death_status) ~ baseline_age, 
                           data = male_survival_data, 
                           dist = "gompertz")

# 2. Female age model (to get age_hazard for females)  
female_age_gompertz <- phreg(Surv(death_date, death_status) ~ baseline_age, 
                             data = female_survival_data, 
                             dist = "gompertz")

# 3. Male slope_ratio model (to get slope_ratio_hazard for males)
male_poa_gompertz <- phreg(Surv(death_date, death_status) ~ slope_ratio, 
                           data = male_survival_data, 
                           dist = "gompertz")

# 4. Female slope_ratio model (to get slope_ratio_hazard for females)
female_poa_gompertz <- phreg(Surv(death_date, death_status) ~ slope_ratio, 
                             data = female_survival_data, 
                             dist = "gompertz")

## Extract log hazards  -----------------------------

# Extract age hazards (log hazard ratio per year of age)
male_age_hazard <- male_age_gompertz$coefficients["baseline_age"]
female_age_hazard <- female_age_gompertz$coefficients["baseline_age"]

# Extract slope_ratio hazards (log hazard ratio per unit of slope_ratio)
male_slope_ratio_hazard <- male_poa_gompertz$coefficients["slope_ratio"]
female_slope_ratio_hazard <- female_poa_gompertz$coefficients["slope_ratio"]

## Apply Formula for FraminghamPoA  -----------------------------

# For males
male_survival_data$FraminghamPoA <- ((male_survival_data$slope_ratio - 1) / 
                                       ((30 * male_age_hazard) / male_slope_ratio_hazard)) + 1

# For females  
female_survival_data$FraminghamPoA <- ((female_survival_data$slope_ratio - 1) / 
                                         ((30 * female_age_hazard) / female_slope_ratio_hazard)) + 1

## Create FraminghamPoA Dataset  -----------------------------

FraminghamPoA <- rbind(male_survival_data, female_survival_data)

FraminghamPoA <- FraminghamPoA %>%
  select(subject_id, FraminghamPoA)

## Export FraminghamPoA Dataset  -----------------------------

write.csv(FraminghamPoA, "../Data/WTM_FraminghamPoA.csv")
