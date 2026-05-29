# Analyze cognition data

library(readr)
library(glmnet)
library(ordinal)
library(nnet)
library(MASS)

cognition_data <- read_csv("../HRS_Data/HRS_cognition_data.csv")
cognition_data_long <- read_csv("../HRS_Data/HRS_cognition_data_long.csv")

cognition_data$Race <- as.factor(cognition_data$Race)
cognition_data$Female <- as.factor(cognition_data$Female)


#######################################
# Multinomial logit function
#######################################

multi_logit_restructured <- function(data, clocks, control_covariates){
  
  all_results = data.frame(
    clock = character(),
    outcome = character(),
    hazard_ratio = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    p_value = numeric(),
    n_subjects = numeric(),
    n_events = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get the total number of subjects
  n_subjects <- nrow(data)
  
  for (clock in clocks){
    
    CP_formula <- as.formula(paste("cogfunction2020", " ~ cogfunction2016 + ",
                                   paste(c(clock, control_covariates), collapse = " + ")))
    
    CP_model <- suppressWarnings(multinom(CP_formula, data = data))
    
    model_data <- model.frame(CP_formula, data = data)
    n_subjects <- nrow(model_data)
    
    complete_data <- model_data[complete.cases(model_data), ]
    
    # Count transitions from lower to higher states (worsening)
    # For CIND: count cases where cogfunction2020 is 2 and cogfunction2016 was 1
    n_events_cind <- sum(complete_data$cogfunction2020 == 2 & complete_data$cogfunction2016 == 1, na.rm = TRUE)
    
    # For Dementia: count cases where cogfunction2020 is 3 and cogfunction2016 was 1 or 2
    n_events_dementia <- sum(complete_data$cogfunction2020 == 3 & 
                               (complete_data$cogfunction2016 == 1 | complete_data$cogfunction2016 == 2), 
                             na.rm = TRUE)
    
    # Extract coefficients and SEs for Level 2 vs Level 1
    CP_coeff_2 <- coef(CP_model)[1, clock]  # First row is Level 2 vs reference
    CP_se_2 <- summary(CP_model)$standard.errors[1, clock]
    
    # Calculate odds ratios and confidence intervals
    CP_OR_2 <- exp(CP_coeff_2)
    CP_lower_2 <- exp(CP_coeff_2 - (1.96 * CP_se_2))
    CP_upper_2 <- exp(CP_coeff_2 + (1.96 * CP_se_2))
    CP_z_2 <- CP_coeff_2 / CP_se_2
    CP_pval_2 <- 2 * pnorm(-abs(CP_z_2))
    
    # Add the Level 2 vs Level 1 results to the dataframe
    level2_results <- data.frame(
      clock = clock,
      outcome = "CIND",
      hazard_ratio = CP_OR_2,
      lower_ci = CP_lower_2,
      upper_ci = CP_upper_2,
      p_value = CP_pval_2,
      n_subjects = n_subjects,
      n_events = n_events_cind,
      stringsAsFactors = FALSE
    )
    
    all_results <- rbind(all_results, level2_results)
    
    # Extract coefficients and SEs for Level 3 vs Level 1
    CP_coeff_3 <- coef(CP_model)[2, clock]  # Second row is Level 3 vs reference
    CP_se_3 <- summary(CP_model)$standard.errors[2, clock]
    CP_z_3 <- CP_coeff_3 / CP_se_3
    CP_pval_3 <- 2 * pnorm(-abs(CP_z_3))
    
    # Calculate odds ratios and confidence intervals
    CP_OR_3 <- exp(CP_coeff_3)
    CP_lower_3 <- exp(CP_coeff_3 - (1.96 * CP_se_3))
    CP_upper_3 <- exp(CP_coeff_3 + (1.96 * CP_se_3))
    
    
    # Add the Level 3 vs Level 1 results to the dataframe
    level3_results <- data.frame(
      clock = clock,
      outcome = "Dementia",
      hazard_ratio = CP_OR_3,
      lower_ci = CP_lower_3,
      upper_ci = CP_upper_3,
      p_value = CP_pval_3,
      n_subjects = n_subjects,
      n_events = n_events_dementia,
      stringsAsFactors = FALSE
    )
    
    all_results <- rbind(all_results, level3_results)
    
  }
  
  return(all_results)
  
}


#################################
# Get multinomial results
#################################

clocks <- c("RidgeResid_scaled",
            "ElasticNetResid_scaled",
            "DunedinPACEResid_scaled", 
            "PCGrimAgeResid_scaled", 
            "PCPhenoAgeResid_scaled")


## With minimal covariates
minimal_covariates <- c("PAGE", "Female", "Race", "Ethnicity")
multi_results_overall <- multi_logit_restructured(cognition_data, clocks, minimal_covariates)
write.csv(multi_results_overall, "../HRS_Output/HRS_Cognition_Results_Minimal.csv")


## With cell controlss
minimal_plus_cells <- c("Female", "PAGE", "Race", "Ethnicity",
                        "Bas", "Bmem", "Bnv", "CD4mem",
                        "CD4nv", "CD8mem", "CD8nv", "Eos",
                        "Mono", "Neu", "NK", "Treg")
multi_results_with_cells <- multi_logit_restructured(cognition_data, clocks, minimal_plus_cells)
write.csv(multi_results_with_cells, "../HRS_Output/HRS_Cognition_Results_Cells.csv")


## With bmi
minimal_plus_bmi <- c("Female", "PAGE", "Race", "Ethnicity",
                      "bmi1213")
multi_results_with_bmi <- multi_logit_restructured(cognition_data, clocks, minimal_plus_bmi)
write.csv(multi_results_with_bmi, "../HRS_Output/HRS_Cognition_Results_BMI.csv")


## With smoking
minimal_plus_smoking <- c("Female", "PAGE", "Race", "Ethnicity",
                          "smoker1213", "pwavessmoked", "n_smokedata")
multi_results_with_smoking <- multi_logit_restructured(cognition_data, clocks, minimal_plus_smoking)
write.csv(multi_results_with_smoking, "../HRS_Output/HRS_Cognition_Results_Smoking.csv")

