### Calculate Pace of Aging with Hierarchical Modeling

## Set Up  --------------------------------------------

library(readr)
library(nlme)

## Load Data  --------------------------------------------

biomarkers_final <- read_csv("../Data/WTM_biomarkers_final.csv")

## Fit Models  --------------------------------------------

quadHM <- lme(
  fixed = value ~ time_5year + time_5year_squared,
  random = ~ time_5year + time_5year_squared | subject_id/biomarker,
  weights = varIdent(form = ~1 | biomarker),
  data = biomarkers_final,
  method = "REML",
  control = lmeControl(opt = "nlminb", msMaxIter = 10000, msMaxEval = 10000)
)


## Extract Slopes  --------------------------------------------

fixed_quadHM_slope <- fixed.effects(quadHM)[2]
random_quadHM_slopes <- random.effects(quadHM)[["subject_id"]][["time_5year"]]
quadHM_slopes <- fixed_quadHM_slope + random_quadHM_slopes


## Create and Export Slopes --------------------------------------------

PoA <- data.frame(
  subject_id = rownames(random.effects(quadHM)[["subject_id"]]),
  quadHM_slopes = quadHM_slopes
)

write.csv(PoA, file="../Data/WTM_PoA.csv", row.names = FALSE)


## Create Dataset of all Weights --------------------------------------------

# Extract Fixed Effects
fixed_effects <- fixed.effects(quadHM)
fixed_intercept <- fixed_effects[1]  # Fixed_0
fixed_linear <- fixed_effects[2]     # Fixed_1  
fixed_quadratic <- fixed_effects[3]  # Fixed_2

# Extract Subject-level Random Effects
subject_random <- random.effects(quadHM)[["subject_id"]]

# Extract Biomarker-level Random Effects, nested within subjects
biomarker_random <- random.effects(quadHM)[["biomarker"]]

# Get all unique biomarkers and subjects
all_biomarkers <- unique(biomarkers_final$biomarker)
all_subjects <- unique(biomarkers_final$subject_id)

# Create base dataframe with fixed effects
model_weights <- data.frame(
  subject_id = all_subjects,
  Fixed_0 = fixed_intercept,
  Fixed_1 = fixed_linear,
  Fixed_2 = fixed_quadratic,
  Theta_0 = NA,
  Theta_1 = NA,
  Theta_2 = NA
)

# Fill in subject random effects
for(i in 1:nrow(model_weights)) {
  subj_id <- as.character(model_weights$subject_id[i])
  if(subj_id %in% rownames(subject_random)) {
    model_weights$Theta_0[i] <- subject_random[subj_id, "(Intercept)"]
    model_weights$Theta_1[i] <- subject_random[subj_id, "time_5year"]
    model_weights$Theta_2[i] <- subject_random[subj_id, "time_5year_squared"]
  }
}

# Add columns for each biomarker's delta values
for(biomarker in all_biomarkers) {
  # Initialize columns with NAs
  model_weights[[paste0("d0_", biomarker)]] <- NA
  model_weights[[paste0("d1_", biomarker)]] <- NA
  model_weights[[paste0("d2_", biomarker)]] <- NA
  
  # Fill in values where subject-biomarker combinations exist
  for(i in 1:nrow(model_weights)) {
    subj_id <- model_weights$subject_id[i]
    biomarker_key <- paste(subj_id, biomarker, sep = "/")
    
    if(biomarker_key %in% rownames(biomarker_random)) {
      model_weights[[paste0("d0_", biomarker)]][i] <- biomarker_random[biomarker_key, "(Intercept)"]
      model_weights[[paste0("d1_", biomarker)]][i] <- biomarker_random[biomarker_key, "time_5year"]
      model_weights[[paste0("d2_", biomarker)]][i] <- biomarker_random[biomarker_key, "time_5year_squared"]
    }
  }
}

# Export the complete weights dataset
write.csv(model_weights, file="../Data/WTM_model_weights.csv", row.names = FALSE)

