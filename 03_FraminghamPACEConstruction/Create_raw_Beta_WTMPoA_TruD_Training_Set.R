# create_raw_beta_PoA_TruD_training_set.R

# Load packages
library(dplyr)

# Read in initial data
betas <- readRDS("../FHS_Data/Initial/GRSet_fully_filtered_bmiq.rds")
pheno_data <- read.csv("../FHS_Data/Initial/framingham_offspring_clean_pheno_dataset.csv")
probe_list_data <- read.csv("../FHS_Data/probelist_sugden_TruD.csv")
FraminghamPoA <- read.csv("../FHS_Data/Initial/WTM_FraminghamPoA.csv")

# Filter betas by reliable (Sugden) probes contained in 450k, EPICv1, and EPICv2
probe_list <- probe_list_data$probe
betas <- betas[rownames(betas) %in% probe_list, ]

# Filter betas by samples that passed Columbia QC (passed_cu_core_QC)
# Transposing beforehand, as it will be required for training
betas <- t(betas)
pheno_passed_QC <- pheno_data[pheno_data$passed_cu_core_QC == "pass",]
betas <- betas[rownames(betas) %in% pheno_passed_QC$Barcode, ]

# Switch to dbGaP_Subject_ID
id_map <- setNames(pheno_data$dbGaP_Subject_ID, pheno_data$Barcode)
betas <- as.data.frame(betas)
rownames(betas) <- id_map[rownames(betas)]


# Filter subjects with PACE
PACE_subjects <- FraminghamPoA$subject_id

betas <- betas[rownames(betas) %in% PACE_subjects , ] # 2069 subjects remaining
pheno_data <- pheno_data[pheno_data$dbGaP_Subject_ID %in% PACE_subjects , ] # 2099 subjects remaining

# The final training set should have Age and Sex as Features
# PACE will be added as the final 'feature'
age_sex_df <- pheno_data %>%
  select(dbGaP_Subject_ID, Age, Female)

# Merge columns, then remove dbGaP_SubjectID as a column, rename Female to Sex
matched_age_sex <- age_sex_df[match(rownames(betas), age_sex_df$dbGaP_Subject_ID), ]
complete_X <- cbind(betas, matched_age_sex)

complete_X <- complete_X %>%
  select(-dbGaP_Subject_ID)

complete_data <- cbind(complete_X, FraminghamPoA[match(rownames(betas), FraminghamPoA$subject_id), ])
complete_data <- complete_data %>%
  select(-c(X, subject_id))

# Quick check
head(colnames(complete_data))
tail(colnames(complete_data))
head(rownames(complete_data))
tail(rownames(complete_data))

# Additional quick validation
# Scatterplot with correlation
plot(complete_data$Age, complete_data$FraminghamPoA, 
     xlab = "Age", ylab = "FraminghamPoA")
text(x = min(complete_data$Age), y = max(complete_data$FraminghamPoA), 
     labels = paste("r =", round(cor(complete_data$Age, complete_data$FraminghamPoA), 3)),
     pos = 4)

# Sex-specific histograms with means
library(ggplot2)
ggplot(complete_data, aes(x = FraminghamPoA)) +
  geom_histogram(bins = 30) +
  geom_vline(data = aggregate(FraminghamPoA ~ Female, complete_data, mean),
             aes(xintercept = FraminghamPoA), color = "red") +
  facet_wrap(~Female) +
  theme_minimal()

# Save the data
write.csv(complete_data, "../FHS_Data/WTM_FraminghamPoA_Training_Set_v2_TruD.csv", row.names = TRUE)

