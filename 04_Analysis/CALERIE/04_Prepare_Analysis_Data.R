################################
# Load libraries
################################

library(tidyverse)
library(sjlabelled)
library(haven)
library(sjPlot)
library(labelled)
library(lme4)
library(rmdformats)
library(ggpubr)

################################
# Define clock variables
################################

# Define all FraminghamPACE variables
framingham_vars <- c(
  "ElasticNet", "Ridge"
)

# Other clock variables
other_clock_vars <- c("DunedinPACE", "PCGrimAge")

# All clock variables combined
all_clock_vars <- c(framingham_vars, other_clock_vars)

################################
# Load data
################################

combined_data <- read_csv("../CALERIE_Data/CALERIE_Combined_Data.csv")
survey <- haven::read_dta("../CALERIE_Data/Initial/CALERIE_DWB220314.dta")

# ## Quick check for MS
# baseline <- combined_data[combined_data$fu == 0,]
# treatment <- sum(baseline$CR == 1)
# control <- sum(baseline$CR == 0)

added_vars <- haven::read_dta("../CALERIE_Data/Initial/CALERIENatAgingPhenoData.dta")

clean_added_vars <- added_vars %>%
  dplyr::select(deidnum, fu, race3, agedwb)

# Merge with survey data on deidnum + fu
survey <- survey %>%
  left_join(clean_added_vars, by = c("deidnum", "fu"))

# Select clock variables dynamically
clocks <- combined_data %>%
  dplyr::select(BARCODE, all_of(all_clock_vars)) %>%
  dplyr::rename(barcode = BARCODE)

survey <-survey %>% 
  mutate(Visit.ID = if_else(fu == 0, "Baseline", 
                            if_else(fu == 1, "12 Month", "24 Month"))) %>% 
  mutate(Treatment = if_else(CR == 0, "Ad Libitum", "Caloric Restriction")) %>% 
  dplyr::select(barcode, deidnum, Treatment, Visit.ID, cbage, agedwb,
                fu, Treatment, race3, deidsite, sex, bmistrat, asample,
                epiccd4t_br, epiccd8t_br, epicnk_br, epicbcell_br, epicmono_br, 
                epicneu_br, methpc1:methpc24) %>% filter(barcode != "")

df <- left_join(survey, clocks, by="barcode") %>%
  sjlabelled::as_factor(Visit.ID) %>%
  mutate(Visit.ID = fct_relevel(Visit.ID, "Baseline", "12 Month", "24 Month"))

# Filter for complete cases across ALL clock variables
clock_vars <- names(clocks)[-1]  # Exclude 'barcode' column
# df <- df[complete.cases(df[clock_vars]), ]

# Filter for asample==1 to match STATA pipeline
# Remove?
# df <- df %>% 
#   filter(asample == 1)

##############################################
# Process the clocks
##############################################

# Get baseline age
df <-df %>% 
  mutate(baseline_age = (cbage*5)+38) 


# 3. Batch correct the clocks using random effects approach
methpc_vars <- colnames(df %>% dplyr::select(contains("methpc")))
f_fixed <- paste(methpc_vars, collapse = ' + ')

for (var in all_clock_vars) {
  bc_var <- paste0(var, "_bc")
  if (all(is.na(df[[var]]))) { df[[bc_var]] <- NA_real_; next }
  
  re_model <- lme4::lmer(
    as.formula(paste0(var, " ~ ", f_fixed, " + (1|deidnum)")),
    data = df, na.action = na.exclude
  )
  
  xb_fixed   <- predict(re_model, re.form = NA)                 # intercept + methPCs
  intercept  <- fixef(re_model)[["(Intercept)"]]
  df[[bc_var]] <- df[[var]] - xb_fixed + intercept              # = y - methPC effects
}
##############################################
# Calculate baseline values, and changes from baseline
##############################################

# Create batch corrected variable names (no age residualization)
bc_vars <- paste0(all_clock_vars, "_bc")

# --- Baseline values per ID (no asample restriction here; mirrors Stata b_<y>)
baseline_values <- df %>%
  dplyr::filter(Visit.ID == "Baseline") %>%
  dplyr::select(deidnum, dplyr::all_of(bc_vars)) %>%
  dplyr::rename_with(~ paste0(., "_baseline"), dplyr::all_of(bc_vars))

# Join baseline values back to the original dataframe
df <- df %>%
  dplyr::left_join(baseline_values, by = "deidnum")

# Define all variables that need change calculations
change_vars <- c("epiccd4t_br", "epiccd8t_br", "epicnk_br", "epicbcell_br", 
                 "epicmono_br", "epicneu_br", bc_vars)

# Changes vs each subject's own baseline (if baseline exists)
df <- df %>%
  dplyr::group_by(deidnum) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(change_vars),
      ~ if (any(.data$Visit.ID == "Baseline")) {
        .x - .x[.data$Visit.ID == "Baseline"][1]
      } else {
        NA_real_
      },
      .names = "d_{.col}"
    )
  ) %>%
  dplyr::ungroup()

# ---------- KEY CHANGE ----------
# Compute baseline means/SDs used for scaling ONLY on (Baseline & asample==1),
# and compute them from the original bc_vars (not the *_baseline columns).
baseline_stats_df <- df %>%
  dplyr::filter(Visit.ID == "Baseline", asample == 1)

baseline_sds_named <- baseline_stats_df %>%
  dplyr::summarize(dplyr::across(
    dplyr::all_of(bc_vars), ~ sd(.x, na.rm = TRUE), .names = "{.col}_sd"
  )) %>%
  dplyr::slice_head(n = 1) %>%
  unlist()

baseline_means_named <- baseline_stats_df %>%
  dplyr::summarize(dplyr::across(
    dplyr::all_of(bc_vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"
  )) %>%
  dplyr::slice_head(n = 1) %>%
  unlist()
# --------------------------------

## Scale numericals (non-clock)
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

df_scaled_w_baseline <- df %>%
  # Scale non-clock continuous variables normally (across the full analysis set)
  dplyr::mutate_at(
    c("cbage", "agedwb",
      "d_epiccd4t_br", "d_epiccd8t_br", "d_epicnk_br", 
      "d_epicbcell_br", "d_epicmono_br", "d_epicneu_br"),
    scale2, na.rm = TRUE
  ) %>%
  # Keep categorical variables as factors
  sjlabelled::as_factor(deidsite) %>%
  dplyr::mutate(
    race3   = as.factor(race3),
    bmistrat = as.factor(bmistrat),
    sex     = as.factor(sex)
  )

# STATA-style scaling for clock baseline values:
# (baseline - baseline_mean) / baseline_sd, where mean/sd come from Baseline & asample==1
baseline_sd_vars <- paste0(bc_vars, "_baseline")

for (i in seq_along(baseline_sd_vars)) {
  baseline_col <- baseline_sd_vars[i]            # e.g., "DunedinPACE_bc_baseline"
  bc_name      <- sub("_baseline$", "", baseline_col)  # e.g., "DunedinPACE_bc"
  sd_name      <- paste0(bc_name, "_sd")         # from baseline_sds_named
  mean_name    <- paste0(bc_name, "_mean")       # from baseline_means_named
  
  sd0  <- as.numeric(baseline_sds_named[sd_name])
  mu0  <- as.numeric(baseline_means_named[mean_name])
  
  # Safe guard for zero/NA SD
  if (is.na(sd0) || sd0 == 0) {
    df_scaled_w_baseline[[baseline_col]] <- NA_real_
  } else {
    df_scaled_w_baseline[[baseline_col]] <-
      (df_scaled_w_baseline[[baseline_col]] - mu0) / sd0
  }
}

# STATA-style scaling for clock change variables: change / baseline_sd (no mean centering)
d_bc_vars <- paste0("d_", bc_vars)

for (i in seq_along(d_bc_vars)) {
  dcol   <- d_bc_vars[i]              # e.g., "d_DunedinPACE_bc"
  bc_name <- bc_vars[i]               # e.g., "DunedinPACE_bc"
  sd_name <- paste0(bc_name, "_sd")   # from baseline_sds_named
  
  sd0 <- as.numeric(baseline_sds_named[sd_name])
  
  if (is.na(sd0) || sd0 == 0) {
    df_scaled_w_baseline[[dcol]] <- NA_real_
  } else {
    df_scaled_w_baseline[[dcol]] <- df_scaled_w_baseline[[dcol]] / sd0
  }
}

# Filter out Baseline (match Stata's "if fu>0")
df_scaled <- df_scaled_w_baseline %>%
  dplyr::filter(Visit.ID != "Baseline")

# (Optional) sanity checks / counts — do NOT subset analysis to require both timepoints
df_scaled %>%
  dplyr::group_by(deidnum) %>%
  dplyr::filter(dplyr::n() >= 2) %>%
  dplyr::distinct(deidnum) %>% nrow()

## Factor conversions
df_scaled_w_baseline <- df_scaled_w_baseline %>%
  sjlabelled::as_factor(Visit.ID) %>%
  dplyr::mutate(Visit.ID = forcats::fct_relevel(Visit.ID, "Baseline", "12 Month", "24 Month")) %>%
  sjlabelled::as_factor(deidsite)

df_scaled <- df_scaled %>%
  sjlabelled::as_factor(Visit.ID) %>%
  dplyr::mutate(Visit.ID = forcats::fct_relevel(Visit.ID, "12 Month", "24 Month")) %>%  # No baseline since filtered out
  sjlabelled::as_factor(deidsite)


# ##############################################
# # For plotting
# ##############################################
# 
# df_plotting <- df_scaled_w_baseline %>% 
#   dplyr::select(barcode, deidnum, Visit.ID, Treatment, all_of(bc_vars))

##############################################
# Export data
##############################################

saveRDS(df_scaled, "../CALERIE_Data/CALERIE_Analysis_Data.rds")
saveRDS(df_scaled_w_baseline, "../CALERIE_Data/CALERIE_Analysis_Data_for_Plotting.rds")
