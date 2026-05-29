# Assess HRS results

library(readr)
library(survival)
library(dplyr)
library(eha)

survival_data <- read_csv("../HRS_Data/HRS_complete_survival_data.csv")

survival_data$Race <- as.factor(survival_data$Race)
survival_data$Female <- as.factor(survival_data$Female)
survival_data$Ethnicity <- as.factor(survival_data$Ethnicity)


#####################################

#####################################
# Survival Analyses
#####################################

#####################################

run_survival_analysis <- function(data, clocks, outcome_time, outcome_status, control_covariates){
  
  
  all_results <- data.frame(
    Dataset = character(),
    outcome = character(),
    clock = character(),
    hazard_ratio = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    p_value = numeric(),
    n_subjects = numeric(),
    n_events = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (clock in clocks){
    
    formula <- as.formula(paste("Surv(", outcome_time, ",", outcome_status, ") ~", 
                                paste(c(clock, control_covariates), collapse = " + ")))
    
    model <- suppressWarnings(coxph(formula, data = data))
    
    n_subjects <- model$n[1]  # The first element of model$n gives the number of subjects
    n_events <- model$nevent  # This gives the number of events in the model
    
    coeffs <- summary(model)[["coefficients"]]
    
    hr <- coeffs[clock, "exp(coef)"]
    p_val <- coeffs[clock, "Pr(>|z|)"]
    
    conf_ints <- confint(model)
    lower_ci <- exp(conf_ints[clock, "2.5 %"])
    upper_ci <- exp(conf_ints[clock, "97.5 %"])
    
    # Add to results dataframe
    clock_results <- data.frame(
      Dataset = "HRS",
      outcome = "Death",
      clock = clock,
      hazard_ratio = hr,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      p_value = p_val,
      n_subjects = n_subjects,
      n_events = n_events,
      stringsAsFactors = FALSE
    )
    
    all_results <- rbind(all_results, clock_results)
    
  }
  
  # Print formatted table at the end
  for(i in 1:nrow(all_results)) {
    cat(sprintf("%-25s HR: %.4f (95%% CI: %.4f-%.4f)\n", 
                all_results$clock[i], 
                all_results$hazard_ratio[i], 
                all_results$lower_ci[i], 
                all_results$upper_ci[i]))
  }
  cat("\n")
  
  return(all_results)
  
}
#####################################
# Outcome: Death / Population: All
#####################################

## Minimal covariates
minimal_covariates <- c("Female", "PAGE", "Race", "Ethnicity")

clocks <- c("RidgeResid_scaled",
            "ElasticNetResid_scaled",
            "DunedinPACEResid_scaled", 
            "PCGrimAgeResid_scaled", 
            "PCPhenoAgeResid_scaled")


results_death_all_nocells <- run_survival_analysis(
  data = survival_data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = minimal_covariates
)

write.csv(results_death_all_nocells, "../HRS_Output/HRS_Survival_Results_Minimal.csv")


## Minimal covariates + cells

minimal_plus_cells <- c("Female", "PAGE", "Race", "Ethnicity",
                        "Bas", "Bmem", "Bnv", "CD4mem",
                        "CD4nv", "CD8mem", "CD8nv", "Eos",
                        "Mono", "Neu", "NK", "Treg")

results_death_cells <- run_survival_analysis(
  data = survival_data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = minimal_plus_cells
)

write.csv(results_death_cells, "../HRS_Output/HRS_Survival_Results_Cells.csv")


## Minimal covariates + bmi

minimal_plus_bmi <- c("Female", "PAGE", "Race", "Ethnicity",
                        "bmi1213")

results_death_bmi <- run_survival_analysis(
  data = survival_data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = minimal_plus_bmi
)

write.csv(results_death_bmi, "../HRS_Output/HRS_Survival_Results_BMI.csv")


## Minimal covariates + smoking

minimal_plus_smoking <- c("Female", "PAGE", "Race", "Ethnicity",
                      "smoker1213", "pwavessmoked", "n_smokedata")

results_death_smoking <- run_survival_analysis(
  data = survival_data,
  clocks = clocks,
  outcome_time = "death_time",
  outcome_status = "death_status",
  control_covariates = minimal_plus_smoking
)

write.csv(results_death_smoking, "../HRS_Output/HRS_Survival_Results_Smoking.csv")

