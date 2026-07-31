# Load packages
library(readr)
library(dplyr)
library(plm)
library(broom)
library(ggplot2)
library(tidyr)
library(stringr)

# Load data
full_data <- read_csv("../CARDIA_Data/CARDIA_combined_data.csv")

# Create a proper time variable
full_data <- full_data %>%
  group_by(dbgap_subject_id) %>%
  mutate(time_since_baseline = collection_year - min(collection_year)) %>%
  ungroup()

# Remove technical replicates
clock_cols <- c("ElasticNet", 
                "Ridge",
                "DunedinPACE",
                "PCGrimAge",
                "Bas",                       
                "Bmem", "Bnv",                
                "CD4mem", "CD4nv",
                "CD8mem", "CD8nv",                    
                "Eos", "Mono",                       
                "Neu", "NK",                      
                "Treg")

unique_cols <- c("sample_name", "barcode")

clean_data <- full_data %>%
  group_by(dbgap_subject_id, collection_visit) %>%
  summarise(
    # Take first instance of unique identifier columns
    across(all_of(unique_cols), first),
    
    # Take mean of clock columns across technical replicates (ignoring NA values)
    across(all_of(clock_cols), ~ mean(.x, na.rm = TRUE)),
    
    # For any remaining columns not specified above, take the first value
    # (assuming they should be the same across technical replicates)
    across(everything(), first),
    
    # Count number of technical replicates for this subject-visit combination
    n_tech_replicates = n(),
    
    .groups = "drop"
  )


#########################################################
# Get complete timepoint data
#########################################################

# Find individuals who have all 4 timepoints
complete_individuals <- clean_data %>%
  group_by(dbgap_subject_id) %>%
  summarise(
    timepoints = n_distinct(collection_visit),
    visits = paste(sort(unique(collection_visit)), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(timepoints == 4) %>%
  pull(dbgap_subject_id)

# Filter to only those individuals
complete_timepoint_data <- clean_data %>%
  filter(dbgap_subject_id %in% complete_individuals)

# Quick check
cat("Original data:", length(unique(clean_data$dbgap_subject_id)), "individuals\n")
cat("Complete timepoint data:", length(unique(complete_timepoint_data$dbgap_subject_id)), "individuals\n")
cat("Rows in complete data:", nrow(complete_timepoint_data), "\n")

# Verify each person has all 4 visits
complete_timepoint_data %>%
  group_by(dbgap_subject_id) %>%
  summarise(visits = n_distinct(collection_visit)) %>%
  pull(visits) %>%
  table()


#########################################################
# Function to run fixed effects
#########################################################

run_clock_analysis <- function(data, clock_column, quadratic=TRUE) {
  
  # Create formula based on quadratic argument
  if (quadratic) {
    formula_str <- paste(clock_column, "~ time_since_baseline + I(time_since_baseline^2)")
  } else {
    formula_str <- paste(clock_column, "~ time_since_baseline")
  }
  formula_obj <- as.formula(formula_str)
  
  # Run fixed effects model
  model <- plm(formula_obj,
               data = data,
               index = c("dbgap_subject_id", "collection_year"),
               model = "within")
  
  # Extract results
  model_summary <- summary(model)
  
  # Get coefficient info for both terms
  linear_coef <- model_summary$coefficients["time_since_baseline", ]
  
  # Calculate SD from first wave (Y15) data
  first_wave_sd <- data %>%
    filter(collection_visit == "Y15") %>%
    pull(!!sym(clock_column)) %>%
    sd(na.rm = TRUE)
  
  # Scale the linear slope by first wave SD
  scaled_linear_slope <- linear_coef["Estimate"] / first_wave_sd
  scaled_ci_lower <- (linear_coef["Estimate"] - 1.96 * linear_coef["Std. Error"]) / first_wave_sd
  scaled_ci_upper <- (linear_coef["Estimate"] + 1.96 * linear_coef["Std. Error"]) / first_wave_sd
  
  
  # Base results tibble
  results <- tibble(
    clock = clock_column,
    linear_slope = linear_coef["Estimate"],
    linear_se = linear_coef["Std. Error"],
    linear_t = linear_coef["t-value"],
    linear_p = linear_coef["Pr(>|t|)"],
    linear_ci_lower = linear_slope - 1.96 * linear_se,
    linear_ci_upper = linear_slope + 1.96 * linear_se,
    scaled_linear_slope = scaled_linear_slope,
    scaled_ci_lower = scaled_ci_lower,
    scaled_ci_upper = scaled_ci_upper,
    first_wave_sd = first_wave_sd,
    n_obs = nobs(model),
    n_individuals = length(unique(data$dbgap_subject_id))
  )
  
  # Add quadratic results if requested
  if (quadratic) {
    quad_coef <- model_summary$coefficients["I(time_since_baseline^2)", ]
    results <- results %>%
      mutate(
        quad_slope = quad_coef["Estimate"],
        quad_se = quad_coef["Std. Error"],
        quad_t = quad_coef["t-value"],
        quad_p = quad_coef["Pr(>|t|)"],
        quad_ci_lower = quad_slope - 1.96 * quad_se,
        quad_ci_upper = quad_slope + 1.96 * quad_se
      )
  } else {
    results <- results %>%
      mutate(
        quad_slope = NA,
        quad_se = NA,
        quad_t = NA,
        quad_p = NA,
        quad_ci_lower = NA,
        quad_ci_upper = NA
      )
  }
  
  return(results)
}

#########################################################
# Function to run all clocks
#########################################################

run_multiple_clocks <- function(data, clock_names, quadratic=TRUE) {
  results_list <- list()
  
  for(i in 1:length(clock_names)) {
    results_list[[i]] <- run_clock_analysis(data, clock_names[i], quadratic=quadratic)
  }
  
  all_results <- do.call(rbind, results_list)
  return(all_results)
}

#########################################################
# Get results
#########################################################

# Usage:
clock_names <- c("ElasticNet", 
                 "Ridge",
                 "DunedinPACE",
                 "PCGrimAge")

# All measurements, no quadratic term
all_measurements_NOquadratic <- run_multiple_clocks(clean_data, clock_names, quadratic=FALSE)

# Complete timepoint measurements, no quadratic term
complete_timepoints_NOquadratic <- run_multiple_clocks(complete_timepoint_data, clock_names, quadratic=FALSE)

# All measurements, with quadratic term
all_measurements_quadratic <- run_multiple_clocks(clean_data, clock_names, quadratic=TRUE)

# Complete timepoint measurements, with quadratic term
complete_timepoints_quadratic <- run_multiple_clocks(complete_timepoint_data, clock_names, quadratic=TRUE)


## Export
write.csv(all_measurements_NOquadratic, "../CARDIA_Output/CARDIA_FixEf_AllMeasurements_NoQuadratic.csv")
write.csv(complete_timepoints_NOquadratic, "../CARDIA_Output/CARDIA_FixEf_CompleteTimepoints_NoQuadratic.csv")
write.csv(all_measurements_quadratic, "../CARDIA_Output/CARDIA_FixEf_AllMeasurements_WithQuadratic.csv")
write.csv(complete_timepoints_quadratic, "../CARDIA_Output/CARDIA_FixEf_CompleteTimepoints_WithQuadratic.csv")