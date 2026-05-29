# FHS Survival analysis

library(readr)
library(stringr)
library(survival)
library(ggplot2)
library(dplyr)

survival_data <- read_csv("../FHS_Data/complete_survival_data.csv")

survival_data$Female <- as.factor(survival_data$Female)
summary(survival_data$death_date / 365.25)

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
    
    coeffs <- summary(model)$coefficients
    conf_ints <- confint(model)
    
    # Get sample size info
    n_subjects <- nrow(data)
    n_events <- sum(data[[outcome_status]], na.rm = TRUE)
    
    # Create row of results
    clock_results <- data.frame(
      Dataset = "FHS",
      outcome = outcome_status,
      clock = clock,
      hazard_ratio = coeffs[clock, "exp(coef)"],
      lower_ci = exp(conf_ints[clock, "2.5 %"]),
      upper_ci = exp(conf_ints[clock, "97.5 %"]),
      p_value = coeffs[clock, "Pr(>|z|)"],
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

run_survival_analyses <- function(data, clocks, outcomes, control_covariates){
  
  combined_results <- data.frame(
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
  
  for (outcome in outcomes){
    
    outcome_status = str_c(outcome, "_status")
    outcome_time = str_c(outcome, "_date")
    
    outcome_results <- run_survival_analysis(
      data = data,
      clocks = clocks,
      outcome_time = outcome_time,
      outcome_status = outcome_status,
      control_covariates = control_covariates
    )
    
    combined_results <- rbind(combined_results, outcome_results)
  }
  
  return(combined_results)
}


clocks <- c("FraminghamPoAResid_scaled", 
            "RidgeResid_scaled", "ElasticNetResid_scaled",
            "DunedinPACEResid_scaled", "PCGrimAgeResid_scaled")

control_covariates <- c("Female", "Age")

outcomes <- c("death", "chd", "chf",
              "cvd", "stroke", "dementia")

results_all_nocells <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = control_covariates
)

write.csv(results_all_nocells, "../FHS_Output/fhs_overall_survival_results.csv")
