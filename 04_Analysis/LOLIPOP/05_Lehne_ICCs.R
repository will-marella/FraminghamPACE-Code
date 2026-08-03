library(tidyverse)
library(rptR)

combined_data <- read_csv("../Lehne_Data/Lehne_Replicate_Combined_Data.csv")

length(unique(combined_data$subject_id))

selected_clocks = c("ElasticNet",
                    "Ridge",
                    "DunedinPACE", 
                    "PCGrimAge",
                    "PCPhenoAge"
)


calculate_icc <- function(df, clocks, n_boot = 100) {
  
  # Validate inputs
  if (!all(clocks %in% names(df))) {
    missing_clocks <- clocks[!clocks %in% names(df)]
    stop("Missing clock columns: ", paste(missing_clocks, collapse = ", "))
  }
  
  if (!"subject_id" %in% names(df)) {
    stop("Missing required column: subject_id")
  }
  
  # Transform data to long format
  icc_data <- df |> 
    dplyr::select(subject_id, all_of(clocks)) |>
    pivot_longer(
      cols = all_of(clocks),
      names_to = "clock_name",
      values_to = "clock_value"
    )
  
  # Initialize results tibble
  results <- tibble(
    clock_name = character(length(clocks)),
    icc = numeric(length(clocks)),
    lower_ci = numeric(length(clocks)),
    upper_ci = numeric(length(clocks))
  )
  
  # Calculate ICC for each clock
  for (i in seq_along(clocks)) {
    
    clock_data <- icc_data |>
      filter(clock_name == clocks[i])
    
    # Calculate ICC using rptR
    icc_result <- rpt(
      formula = clock_value ~ (1|subject_id),
      grname = "subject_id",
      data = clock_data,
      datatype = "Gaussian",
      nboot = n_boot,
      npermut = 0
    )
    
    # Extract results
    icc_summary <- summary(icc_result)
    
    results$clock_name[i] <- clocks[i]
    results$icc[i] <- round(as.numeric(icc_summary$R$subject_id), 3)
    results$lower_ci[i] <- round(as.numeric(icc_summary$CI_emp$`2.5%`), 3)
    results$upper_ci[i] <- round(as.numeric(icc_summary$CI_emp$`97.5%`), 3)
  }
  
  return(results)
}

# Example usage:
selected_clocks_resid <- paste0(selected_clocks, "Resid")

icc_results <- calculate_icc(
  df = combined_data,
  clocks = selected_clocks_resid,
  n_boot = 1000
)

# Write to CSV
write.csv(icc_results, "../Lehne_Output/Lehne_ICC_Result_Summary.csv")

