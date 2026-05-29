# Analyze disease/disability data

library(readr)
library(glmnet)

disease_and_disability_long <- read_csv("../HRS_Data/HRS_disease_disability_long_data.csv")
disease_and_disability <- read_csv("../HRS_Data/HRS_disease_disability_data.csv")

disease_and_disability$Race <- as.factor(disease_and_disability$Race)
disease_and_disability$Female <- as.factor(disease_and_disability$Female)

disease_and_disability_long$Race <- as.factor(disease_and_disability_long$Race)
disease_and_disability_long$Female <- as.factor(disease_and_disability_long$Female)


######################################################
# Function to get all poission results
######################################################

get_pois_results <- function(data, clocks, control_covariates, conf_level = 0.95){
  
  all_results = data.frame(
    Outcome = character(),
    Clock = character(),
    rate_ratio = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    p_value = numeric(),
    n_subjects = numeric(),
    n_events = numeric()
  )
  
  for (clock in clocks){
    
    #################################################
    # Chronic Disease - Conditional Prospective only
    #################################################
    
    chron_CP_formula <- as.formula(paste("chrondxe15", " ~ ", "chrondxe13 + ", 
                                         paste(c(clock, control_covariates), collapse = " + ")))
    
    chron_CP_model <- suppressWarnings(glm(chron_CP_formula, family = poisson(link = "log"), data = data))
    
    # Get coefficient and CI
    chron_CP_coeff <- exp(chron_CP_model[["coefficients"]][[clock]])
    chron_CP_ci <- suppressWarnings(exp(confint(chron_CP_model, parm = clock, level = conf_level)))
    chron_CP_pval <- summary(chron_CP_model)$coefficients[clock, "Pr(>|z|)"]
    
    n_subjects_chron <- nrow(chron_CP_model$model)
    
    complete_data_chron <- chron_CP_model$model
    n_events_chron <- sum(complete_data_chron$chrondxe15 > complete_data_chron$chrondxe13, na.rm = TRUE)
    
    # Create chronic disease row
    clock_CD_results <- data.frame(
      Outcome = "Chronic Disease",
      Clock = clock,
      rate_ratio = chron_CP_coeff,
      lower_ci = chron_CP_ci[1],
      upper_ci = chron_CP_ci[2],
      p_value = chron_CP_pval,
      n_subjects = n_subjects_chron,
      n_events = n_events_chron
    )
    
    all_results <- rbind(all_results, clock_CD_results)
    
    #################################################
    # ADLs - Conditional Prospective only
    #################################################
    
    ADL_CP_formula <- as.formula(paste("adl5a15", " ~ ", "adl5a13 + ", 
                                       paste(c(clock, control_covariates), collapse = " + ")))
    
    ADL_CP_model <- suppressWarnings(glm(ADL_CP_formula, family = poisson(link = "log"), data = data))
    
    # Get coefficient and CI
    ADL_CP_coeff <- exp(ADL_CP_model[["coefficients"]][[clock]])
    ADL_CP_ci <- suppressWarnings(exp(confint(ADL_CP_model, parm = clock, level = conf_level)))
    ADL_CP_pval <- summary(ADL_CP_model)$coefficients[clock, "Pr(>|z|)"]
    
    n_subjects_ADL <- nrow(ADL_CP_model$model)
    
    complete_data_ADL <- ADL_CP_model$model
    n_events_ADL <- sum(complete_data_ADL$adl5a15 > complete_data_ADL$adl5a13, na.rm = TRUE)
    
    # Create ADL row
    clock_ADL_results <- data.frame(
      Outcome = "ADLs",
      Clock = clock,
      rate_ratio = ADL_CP_coeff,
      lower_ci = ADL_CP_ci[1],
      upper_ci = ADL_CP_ci[2],
      p_value = ADL_CP_pval,
      n_subjects = n_subjects_ADL,
      n_events = n_events_ADL
    )
    
    all_results <- rbind(all_results, clock_ADL_results)
    
    #################################################
    # IADLs - Conditional Prospective only
    #################################################
    
    IADL_CP_formula <- as.formula(paste("iadl5a15", " ~ ", "iadl5a13 + ", 
                                        paste(c(clock, control_covariates), collapse = " + ")))
    
    IADL_CP_model <- suppressWarnings(glm(IADL_CP_formula, family = poisson(link = "log"), data = data))
    
    # Get coefficient and CI
    IADL_CP_coeff <- exp(IADL_CP_model[["coefficients"]][[clock]])
    IADL_CP_ci <- suppressWarnings(exp(confint(IADL_CP_model, parm = clock, level = conf_level)))
    IADL_CP_pval <- summary(IADL_CP_model)$coefficients[clock, "Pr(>|z|)"]
    
    n_subjects_IADL <- nrow(IADL_CP_model$model)
    
    complete_data_IADL <- IADL_CP_model$model
    n_events_IADL <- sum(complete_data_IADL$iadl5a15 > complete_data_IADL$iadl5a13, na.rm = TRUE)
    
    # Create IADL row
    clock_IADL_results <- data.frame(
      Outcome = "IADLs",
      Clock = clock,
      rate_ratio = IADL_CP_coeff,
      lower_ci = IADL_CP_ci[1],
      upper_ci = IADL_CP_ci[2],
      p_value = IADL_CP_pval,
      n_subjects = n_subjects_IADL,
      n_events = n_events_IADL
    )
    
    all_results <- rbind(all_results, clock_IADL_results)
  }
  
  return(all_results)
}


###########################
# Get results, minimal
###########################

clocks <- c("RidgeResid_scaled",
            "ElasticNetResid_scaled",
            "DunedinPACEResid_scaled", 
            "PCGrimAgeResid_scaled", 
            "PCPhenoAgeResid_scaled")



## Minimal covariates
minimal_covariates <- c("PAGE", "Female", "Race", "Ethnicity")

results_minimal <- get_pois_results(disease_and_disability, clocks, minimal_covariates)

write.csv(results_minimal, "../HRS_Output/HRS_Disease_Disability_Results_Minimal.csv")


###########################
# Get results, with cells
###########################

minimal_plus_cells <- c("Female", "PAGE", "Race", "Ethnicity",
                        "Bas", "Bmem", "Bnv", "CD4mem",
                        "CD4nv", "CD8mem", "CD8nv", "Eos",
                        "Mono", "Neu", "NK", "Treg")

results_cells <- get_pois_results(disease_and_disability, clocks, minimal_plus_cells)

write.csv(results_cells, "../HRS_Output/HRS_Disease_Disability_Results_Cells.csv")

###########################
# Get results, with BMI
###########################

minimal_plus_bmi <- c("Female", "PAGE", "Race", "Ethnicity",
                      "bmi1213")

results_bmi <- get_pois_results(disease_and_disability, clocks, minimal_plus_bmi)

write.csv(results_bmi, "../HRS_Output/HRS_Disease_Disability_Results_BMI.csv")

###########################
# Get results, with smoking
###########################

minimal_plus_smoking <- c("Female", "PAGE", "Race", "Ethnicity",
                          "smoker1213", "pwavessmoked", "n_smokedata")

results_smoking <- get_pois_results(disease_and_disability, clocks, minimal_plus_smoking)

write.csv(results_smoking, "../HRS_Output/HRS_Disease_Disability_Results_Smoking.csv")
