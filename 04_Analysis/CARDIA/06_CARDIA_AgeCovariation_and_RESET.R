# Load packages
library(readr)
library(dplyr)
library(purrr)
library(lmtest)

# Load data
CARDIA_combined_data <- read_csv("../CARDIA_Data/CARDIA_combined_data.csv")

#########################################################
# Function to get RESET test results
#########################################################

test_linearity <- function(data, clock_names, age_var = "age_at_collection") {
  
  # Load required package
  library(lmtest)
  
  results_list <- list()
  
  for(i in seq_along(clock_names)) {
    clock <- clock_names[i]
    
    # Extract the clock values and age variable
    clock_values <- data[[clock]]
    age <- data[[age_var]]
    
    # Remove missing values
    complete_cases <- !is.na(clock_values) & !is.na(age)
    clock_clean <- clock_values[complete_cases]
    age_clean <- age[complete_cases]
    
    # Fit linear model
    lm_model <- lm(clock_clean ~ age_clean)
    
    # Get basic linear model stats
    model_summary <- summary(lm_model)
    r_squared <- model_summary$r.squared
    age_coef <- model_summary$coefficients[2, 1]  # slope
    age_pvalue <- model_summary$coefficients[2, 4]  # p-value for age
    
    # RESET test
    reset_result <- resettest(lm_model)
    
    # Calculate correlation
    correlation <- cor(clock_clean, age_clean)
    
    results_list[[i]] <- data.frame(
      clock = clock,
      n_obs = length(clock_clean),
      correlation = round(correlation, 3),
      r_squared = round(r_squared, 3),
      age_slope = round(age_coef, 4),
      age_p_value = round(age_pvalue, 6),
      reset_statistic = round(reset_result$statistic, 3),
      reset_p_value = round(reset_result$p.value, 4),
      linearity_supported = ifelse(reset_result$p.value > 0.05, "Yes", "No"),
      age_significant = ifelse(age_pvalue < 0.05, "Yes", "No")
    )
  }
  
  results <- do.call(rbind, results_list)
  return(results)
}

#########################################################
# Get results
#########################################################

clock_names <- c("ElasticNet", 
                 "Ridge",
                 "DunedinPACE", 
                 "PCGrimAge")

linearity_results <- test_linearity(CARDIA_combined_data, clock_names)

#########################################################
# Export results
#########################################################

write.csv(sex_diff_results, "../CARDIA_Output/CARDIA_sex_difference_test_results.csv")
write.csv(linearity_results, "../CARDIA_Output/CARDIA_linearity_test_results.csv")
