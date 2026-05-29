# Get BMI, Pack Years for WHI
library(haven)
library(dplyr)

has_bmi <- read_dta("../WHI_Data/Initial/WHIPhysMeasures.dta")
has_pack_years <- read_dta("../WHI_Data/Initial/WHIBackbone_220930.dta")

length(unique(has_bmi$dbgap_subject_id))
length(unique(has_pack_years$dbgap_subject_id))

# Get the first bmix value for each individual
first_bmi <- has_bmi %>%
  # Ensure data is arranged by ID and follow-up day
  arrange(dbgap_subject_id, f80days) %>%
  # Group by individual
  group_by(dbgap_subject_id) %>%
  # Take the first entry
  dplyr::slice(1) %>%
  # Select only necessary columns
  dplyr::select(dbgap_subject_id, subjid, bmix, f80days) %>%
  ungroup()

# Check the resulting dataset
head(first_bmi)
summary(first_bmi)


# Get the pack years for each individual
pack_years <- has_pack_years %>%
  dplyr::select(dbgap_subject_id, packyrs)

head(pack_years)
summary(pack_years)

######

first_bmi <- first_bmi %>%
  dplyr::select(dbgap_subject_id, subjid, bmix)

additional_pheno <- merge(first_bmi, pack_years, by="dbgap_subject_id")


##### 
# Try to get it to have the barcode
pheno_backbone <- read_csv("../WHI_Data/WHI_complete_pheno_outcome_dataset.csv")

all_pheno <- merge(additional_pheno, pheno_backbone, by="dbgap_subject_id", all=FALSE)

clean_additional_pheno <- all_pheno %>%
  dplyr::select(BARCODE, dbgap_subject_id, bmix, packyrs)

summary(clean_additional_pheno$bmix)
summary(clean_additional_pheno$packyrs)

write.csv(clean_additional_pheno, "../WHI_Data/WHI_additional_pheno.csv")
