# Assess WHI Survival Data

# Load packages
library(readr)
library(survival)
library(dplyr)
library(stringr)

# Load Data
survival_data <- read_csv("../WHI_Data/WHI_survival_data.csv")

# survival_data <- survival_data[complete.cases(survival_data),]
survival_data$RACE <- as.factor(survival_data$RACE)


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
      Dataset = "WHI",
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
                all_results$std_error[i],
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
    
    outcome_status = outcome
    outcome_time = str_c(outcome, "DY")
    
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


#####################################
# Survival results
# Minimal controls (Age, Race)
#####################################

outcomes <- c("CHD", "CHF", "STROKE", "ANYCANCER", "ANYFX", "DEATH")

clocks <- c("RidgeResid_scaled",
            "ElasticNetResid_scaled",
            "DunedinPACEResid_scaled", 
            "PCGrimAgeResid_scaled", 
            "PCPhenoAgeResid_scaled")

all_noncell_covariates <- c("AGE", "RACE")

results_all_nocells <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = all_noncell_covariates
)

write.csv(results_all_nocells, "../WHI_Output/whi_overall_survival_results.csv")


#####################################
# Survival results
# + cell controls
#####################################

all_cell_covariates <- c("AGE", "RACE", 
                         "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem",
                         "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")


results_all_cells <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = all_cell_covariates
)

write.csv(results_all_cells, "../WHI_Output/whi_survival_results_w_cells.csv")

#####################################
# Survival results
# + smoking history
#####################################

smoking_covariates <- c("AGE", "RACE",
                        "packyrs")

results_all_smoking <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = smoking_covariates
)

write.csv(results_all_smoking, "../WHI_Output/whi_survival_results_w_smoking.csv")

#####################################
# Survival results
# + bmi
#####################################

bmi_covariates <- c("AGE", "RACE",
                        "bmix")

results_all_bmi <- run_survival_analyses(
  data=survival_data,
  clocks=clocks,
  outcomes=outcomes,
  control_covariates = bmi_covariates
)

write.csv(results_all_bmi, "../WHI_Output/whi_survival_results_w_bmi.csv")
