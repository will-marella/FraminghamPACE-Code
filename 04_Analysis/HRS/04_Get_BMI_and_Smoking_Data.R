# Make additional pheno

# Packages
library(haven)
library(readr)
library(dplyr)
library(stringr)

# Load raw data
HRSRand <- read_dta("../HRS_Data/Initial/HRSRAND_Health2020_DWB.dta")

names(HRSRand)

raw_pheno <- HRSRand %>%
  dplyr::select(hhidpn, bmi12, bmi13,
                smoker12, smoker13,
                pwavessmoked, n_smokedata) %>%
  # Add bmi1213 - use bmi13 if available, otherwise use bmi12
  mutate(bmi1213 = ifelse(!is.na(bmi13), bmi13, bmi12)) %>%
  # Add smoker1213 - use smoker13 if available, otherwise use smoker12
  mutate(smoker1213 = ifelse(!is.na(smoker13), smoker13, smoker12))

# Get pheno backbone
pheno_backbone <- read_csv("../HRS_Data/Initial/HRS_combined_clean_pheno_w_hhid.csv")

# First, if an hhid starts with `0` in pheno_backbone, remove it
pheno_backbone <- pheno_backbone %>%
  mutate(hhid = str_remove(hhid, "^0+"))

pheno_backbone$hhidpn <- paste0(pheno_backbone$hhid, pheno_backbone$pn)


# Combine pheno
combined_pheno <- merge(pheno_backbone, raw_pheno, by="hhidpn", all=FALSE)

clean_additional_pheno <- combined_pheno %>%
  dplyr::select(hhidpn, SubjectID, hhid, pn, IDATid, bmi12, bmi13, smoker13, 
                pwavessmoked, n_smokedata, bmi1213, smoker1213)

write.csv(clean_additional_pheno, "../HRS_Data/HRS_clean_additional_pheno.csv")
