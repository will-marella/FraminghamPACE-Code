# Load packages
library(readr)
library(dplyr)
library(purrr)
library(lmtest)

# Clock variables included in the CARDIA analysis
clock_names <- c(
  "ElasticNet",
  "Ridge",
  "DunedinPACE",
  "PCGrimAge"
)

# Return NA rather than NaN when every replicate value is missing
mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Load data and average technical replicates within participant and visit
clean_data <- read_csv(
  "../CARDIA_Data/CARDIA_combined_data.csv",
  show_col_types = FALSE
) %>%
  group_by(dbgap_subject_id, collection_visit) %>%
  summarise(
    age_at_collection = first(age_at_collection),
    across(all_of(clock_names), mean_or_na),
    .groups = "drop"
  )

#########################################################
# Pearson correlation, OLS model, and Ramsey RESET test
#########################################################

test_linearity <- function(
    data,
    clock,
    age_var = "age_at_collection"
) {
  
  complete_cases <- !is.na(data[[clock]]) & !is.na(data[[age_var]])
  analysis_data <- data[complete_cases, ]
  
  lm_model <- lm(
    reformulate(age_var, response = clock),
    data = analysis_data
  )
  
  model_summary <- summary(lm_model)
  age_result <- model_summary$coefficients[age_var, ]
  reset_result <- resettest(lm_model)
  
  tibble(
    clock = clock,
    n_participants = n_distinct(analysis_data$dbgap_subject_id),
    n_observations = nrow(analysis_data),
    
    correlation = cor(
      analysis_data[[clock]],
      analysis_data[[age_var]]
    ),
    
    r_squared = model_summary$r.squared,
    
    age_slope = unname(age_result["Estimate"]),
    age_se = unname(age_result["Std. Error"]),
    age_t = unname(age_result["t value"]),
    age_df = df.residual(lm_model),
    age_p_value = unname(age_result["Pr(>|t|)"]),
    
    reset_statistic = unname(reset_result$statistic),
    reset_df1 = unname(reset_result$parameter[1]),
    reset_df2 = unname(reset_result$parameter[2]),
    reset_p_value = reset_result$p.value
  )
}

linearity_results <- map_dfr(
  clock_names,
  ~ test_linearity(clean_data, .x)
)

write_csv(
  linearity_results,
  "../CARDIA_Output/CARDIA_linearity_test_results.csv"
)

print(linearity_results)