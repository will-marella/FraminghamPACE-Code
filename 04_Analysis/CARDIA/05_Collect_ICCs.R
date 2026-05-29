# Load packages
library(readr)
library(lme4)
library(performance)
library(dplyr)
library(tibble)

# Load data
CARDIA_combined_data <- read_csv("../CARDIA_Data/CARDIA_combined_data.csv")

#########################################################
# Function to compute ICCs
#########################################################

compute_clock_iccs <- function(data, clock_names, id_var = "dbgap_subject_id") {
  results <- lapply(clock_names, function(clock) {
    # Fit model
    formula <- as.formula(paste0(clock, " ~ 1 + (1 | ", id_var, ")"))
    model <- lmer(formula, data = data, REML = TRUE)
    
    # Compute ICCs
    icc_df <- performance::icc(model, ci = TRUE)
    
    # Extract values correctly by row
    tibble(
      clock = clock,
      icc = icc_df$ICC_adjusted[1],
      ci_low = icc_df$ICC_adjusted[2],
      ci_high = icc_df$ICC_adjusted[3],
      n_subjects = length(unique(data[[id_var]])),
      n_obs = nrow(data)
    )
  })
  
  bind_rows(results)
}




#########################################################
# Get results
#########################################################

CARDIA_combined_data$PCGrimAge_noagecoeff <- CARDIA_combined_data$PCGrimAge - (0.139 * CARDIA_combined_data$age_at_collection)

cor(CARDIA_combined_data$age_at_collection, CARDIA_combined_data$PCGrimAge)
cor(CARDIA_combined_data$age_at_collection, CARDIA_combined_data$PCGrimAge_noagecoeff)

clock_names <- c("ElasticNet", 
                 "Ridge", 
                 "DunedinPACE",
                 "PCGrimAge")

icc_results <- compute_clock_iccs(CARDIA_combined_data, clock_names)
print(icc_results)

#########################################################
# Save results
#########################################################

write.csv(icc_results, "../CARDIA_Output/ICC_Results.csv")
