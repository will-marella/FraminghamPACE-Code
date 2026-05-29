##########################################
## Load packages
##########################################

library(glmnet)
library(readr)
library(dplyr)

##########################################
## Load and format data
##########################################

data <- read_csv("../FHS_Data/WTM_FraminghamPoA_Training_Set_v2_TruD.csv")
data <- as.data.frame(data)
rownames(data) <- data$...1
data$`...1` <- NULL
data$Age <- data$Age - 65

# Force sex-specific means to 1
data$FraminghamPoA <- data$FraminghamPoA + (1 - ave(data$FraminghamPoA, data$Female))

print("Sex-specific means after adjustment:")
print(tapply(data$FraminghamPoA, data$Female, mean))

# Check data
head(colnames(data))
tail(colnames(data))
head(rownames(data))
tail(rownames(data))

# Basic data prep
age_sex_cols <- c(14466, 14467)
cpg_cols <- 1:14465

X <- as.matrix(data[, c(cpg_cols, age_sex_cols)])
y <- data[["FraminghamPoA"]]

penalty_factors <- c(rep(1, length(cpg_cols)), rep(0, length(age_sex_cols)))

##########################################
## Initialize training summary
##########################################

# Initialize training summary
training_summary <- data.frame(
  Metric = character(),
  Value = character(),
  stringsAsFactors = FALSE
)

##########################################
## Hyperparameter tuning on 100% of data (10-fold CV)
##########################################

set.seed(123)  # Fixed seed for reproducible hyperparameter selection

cv_results <- cv.glmnet(
  x = X,
  y = y,
  alpha = 0,
  penalty.factor = penalty_factors,
  nfolds = 10,
  standardize = TRUE
)

# Get best lambda, and CV R2
best_lambda <- cv_results$lambda.min
cv_r2 <- 1 - (min(cv_results$cvm) / mean((y - mean(y))^2))

# Get features selected
tuning_model <- glmnet(x = X, y = y, alpha = 0, lambda = best_lambda, 
                       penalty.factor = penalty_factors, standardize = TRUE)
tuning_n_features <- sum(coef(tuning_model) != 0) - 1  # subtract intercept


print(paste("Best lambda:", best_lambda))
print(paste("CV R²:", round(cv_r2, 3)))
print(paste("Features Selected (Tuning):", tuning_n_features))

training_summary <- rbind(training_summary, 
                          data.frame(Metric = "Hyperparameter: Best Lambda", Value = sprintf("%.6f", best_lambda)),
                          data.frame(Metric = "Hyperparameter: 10-Fold CV R²", Value = sprintf("%.3f", cv_r2)),
                          data.frame(Metric = "Hyperparameter: Features Selected", Value = as.character(tuning_n_features))
)

##########################################
## Subsample-based ensemble training
##########################################

n_seeds <- 100
subsample_prop <- 0.9
subsample_size <- floor(nrow(X) * subsample_prop)
ensemble_coefficients <- matrix(0, nrow = ncol(X) + 1, ncol = n_seeds)

# Track ensemble metrics
held_out_r2_values <- numeric(n_seeds)
n_features_per_model <- numeric(n_seeds)

print(paste("Training subsampling ensemble models with", subsample_size, "samples each..."))
for(seed in 1:n_seeds) {
  if(seed %% 5 == 0) print(paste("Training model", seed, "of", n_seeds))
  
  set.seed(seed)
  
  # Subsample (90% without replacement)
  sub_indices <- sample(nrow(X), size = subsample_size, replace = FALSE)
  held_out_indices <- setdiff(1:nrow(X), sub_indices)
  
  X_sub <- X[sub_indices, ]
  y_sub <- y[sub_indices]
  X_held_out <- X[held_out_indices, ]
  y_held_out <- y[held_out_indices]
  
  final_fit <- glmnet(
    x = X_sub,
    y = y_sub,
    alpha = 0,
    lambda = best_lambda,
    penalty.factor = penalty_factors,
    standardize = TRUE
  )
  
  ensemble_coefficients[, seed] <- as.vector(coef(final_fit))
  
  # Count features
  n_features_per_model[seed] <- sum(coef(final_fit) != 0) - 1
  
  # Evaluate on held-out set
  held_out_predictions <- predict(final_fit, newx = X_held_out, s = best_lambda)
  held_out_r2_values[seed] <- cor(held_out_predictions, y_held_out)^2
}

# Calculate ensemble metrics
mean_held_out_r2 <- mean(held_out_r2_values)
sd_held_out_r2 <- sd(held_out_r2_values)
mean_features <- mean(n_features_per_model)
sd_features <- sd(n_features_per_model)
min_features <- min(n_features_per_model)
max_features <- max(n_features_per_model)

print(paste("Held-out R² - Mean:", round(mean_held_out_r2, 3), "± SD:", round(sd_held_out_r2, 3)))
print(paste("Features per model - Mean:", round(mean_features, 1), "± SD:", round(sd_features, 1), 
            "Range:", min_features, "-", max_features))

# Add ensemble training results to summary
training_summary <- rbind(training_summary,
                          data.frame(Metric = "Ensemble: Mean Held-out R²", Value = sprintf("%.3f ± %.3f", mean_held_out_r2, sd_held_out_r2)),
                          data.frame(Metric = "Ensemble: Mean Features per Model", Value = sprintf("%.1f ± %.1f", mean_features, sd_features)),
                          data.frame(Metric = "Ensemble: Feature Range", Value = sprintf("%d - %d", min_features, max_features))
)

##########################################
## Average coefficients and finalize training summary
##########################################

final_coefficients <- rowMeans(ensemble_coefficients)
final_intercept <- final_coefficients[1]
final_betas <- final_coefficients[-1]

# Number of non-zero coefficients in the final ensemble model
n_nonzero_final <- sum(final_betas != 0)

# Train a single model on full data for comparison
single_model <- glmnet(
  x = X,
  y = y,
  alpha = 0,
  lambda = best_lambda,
  penalty.factor = penalty_factors,
  standardize = TRUE
)

# Get predictions from both models
single_predictions <- predict(single_model, newx = X, s = best_lambda)
ensemble_predictions <- final_intercept + X %*% final_betas

# Calculate in-sample R²
single_r2 <- cor(single_predictions, y)^2
ensemble_r2 <- cor(ensemble_predictions, y)^2

# Correlation between single model and ensemble predictions
prediction_correlation <- cor(single_predictions, ensemble_predictions)

# Add final model metrics to summary
training_summary <- rbind(training_summary,
                          data.frame(Metric = "Final Model: Features", Value = as.character(n_nonzero_final)),
                          data.frame(Metric = "Final Model: R²", Value = sprintf("%.3f", ensemble_r2)),
                          data.frame(Metric = "Comparison: Single Model R²", Value = sprintf("%.3f", single_r2)),
                          data.frame(Metric = "Comparison: Model Correlation", Value = sprintf("%.4f", prediction_correlation))
)

##########################################
## Export weights, training summary
##########################################

# Create coefficient dataframe with non-zero features only
coef_names <- c("(Intercept)", colnames(X))
nonzero_indices <- which(final_coefficients != 0)

ensemble_results_df <- data.frame(
  feature = coef_names[nonzero_indices],
  coefficient = final_coefficients[nonzero_indices],
  selection_frequency = rowSums(ensemble_coefficients != 0)[nonzero_indices] / n_seeds
)

colnames <- colnames(data)
all_sds <- apply(data, 2, sd, na.rm = TRUE)

sd_df <- data.frame(
  feature = colnames,
  sd = all_sds
)

importance_df <- merge(ensemble_results_df, sd_df, by="feature")

importance_df$importance <- abs(importance_df$coefficient) * importance_df$sd * 100

non_features <- c("Age", "Female")
importance_df <- importance_df[!importance_df$feature %in% non_features,  ]

##########################################
## No Export
## CHECK IMPORTANCE DISTRIBUTION
##########################################

############################################################
# Calculate model importance scores (add after model_weights calculation)
############################################################

model_probes <- importance_df$feature
model_sds <- importance_df$sd
model_weights <- importance_df$coefficient
model_importance <- importance_df$importance

# Inspection of importance scores
cat("=== Model Importance Score Inspection ===\n")
cat(sprintf("Total importance: %.4f\n", sum(model_importance)))
cat(sprintf("Mean importance: %.4f\n", mean(model_importance)))
cat(sprintf("Median importance: %.4f\n", median(model_importance)))
cat(sprintf("Range: %.4f to %.4f\n", min(model_importance), max(model_importance)))

model_importance_over01 <- model_importance[model_importance > 0.1]
hist(model_importance_over01,
     breaks=50)
hist(model_sds)
hist(model_weights, breaks=50)
hist(model_importance, breaks=50)


# Check for outliers (values > 3 SDs from mean)
importance_z <- abs(model_importance - mean(model_importance)) / sd(model_importance)
outliers <- which(importance_z > 3)
if(length(outliers) > 0) {
  cat(sprintf("\nFound %d potential outliers (>3 SD from mean):\n", length(outliers)))
  outlier_info <- data.frame(
    probe = model_probes[outliers],
    coefficient = model_weights[outliers],
    sd = model_sds[outliers],
    importance = model_importance[outliers],
    z_score = importance_z[outliers]
  )
  print(outlier_info)
}

# Show top 10 most important probes
top_probes <- order(model_importance, decreasing = TRUE)[1:10]
cat("\nTop 10 most important probes:\n")
top_info <- data.frame(
  probe = model_probes[top_probes],
  coefficient = model_weights[top_probes],
  sd = model_sds[top_probes],
  importance = model_importance[top_probes]
)
print(top_info)


# Calculate cumulative importance percentages
sorted_importance <- sort(model_importance, decreasing = TRUE)
cumulative_importance <- cumsum(sorted_importance)
total_importance <- sum(model_importance)

# Function to find number of probes for given percentage
find_probes_for_percentage <- function(percentage) {
  target <- total_importance * (percentage / 100)
  which.max(cumulative_importance >= target)
}

# Calculate for your desired percentages
percentages <- c(20, 50, 80, 90, 95)
results <- data.frame(
  percentage = percentages,
  num_probes = sapply(percentages, find_probes_for_percentage),
  proportion_of_total = sapply(percentages, function(p) {
    find_probes_for_percentage(p) / length(model_importance)
  })
)

cat("Cumulative Importance Analysis:\n")
print(results)

# Calculate outlier contribution
outlier_importance <- sum(model_importance[outliers])
outlier_percentage <- outlier_importance / sum(model_importance) * 100

cat(sprintf("Outlier Analysis:\n"))
cat(sprintf("160 outliers contribute: %.4f importance (%.1f%% of total)\n", 
            outlier_importance, outlier_percentage))


#
#
#
#

write.csv(ensemble_results_df, "../Models/WTM_FraminghamPACEa_TruDx_RIDGE_Ensemble_Model_Weights.csv")
write.csv(training_summary, "../FHS_Output/WTM_FraminghamPACEa_TruDx_RIDGE_Ensemble_Training_Summary.csv")
