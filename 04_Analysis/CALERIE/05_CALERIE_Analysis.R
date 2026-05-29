################################
# Packages & options
################################
library(tidyverse)
library(lme4)
library(emmeans)
library(clubSandwich)
library(sjPlot)

# emmeans: avoid lmerTest / pbkrtest; use asymptotic df
emm_options(lmer.df = "asymptotic", pbkrtest.limit = 0, lmerTest.limit = 0)

# Theme (optional)
My_Theme <- theme(
  axis.title.x = element_text(size = 16),
  axis.text.x  = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y  = element_text(size = 14)
)

################################
# Load ONE dataset for both analysis & plotting
################################

dat <- readRDS("../CALERIE_Data/CALERIE_Analysis_Data.rds")


# Stata-matching race recode (from your cross-tab: 3=White, 1=Black, 2=Other)
dat <- dat %>%
  mutate(
    race3_std = recode(race3,
                       `3` = "White",
                       `1` = "Black",
                       `2` = "Other",
                       .default = NA_character_
    ),
    race3_std = factor(race3_std, levels = c("White","Black","Other"))
  )

# Build analysis frame = follow-up only (fu > 0), with consistent references
analysis_df <- dat %>%
  filter(Visit.ID != "Baseline") %>%
  mutate(
    Treatment = forcats::fct_relevel(Treatment, "Ad Libitum", "Caloric Restriction"),
    Visit.ID  = forcats::fct_relevel(Visit.ID, "12 Month", "24 Month"),
    fu_numeric = as.integer(Visit.ID) - 1L,         # 0=12m, 1=24m for random slope
    # ensure common covariates are factors
    deidsite  = as.factor(deidsite),
    bmistrat  = as.factor(bmistrat),
    sex       = as.factor(sex)
  )


################################
# Detour: Confirm n_subj, n_obs
# By CR/AL
################################

## N subjects
dat_asamp1 <- dat[dat$asample == 1,]
n_sub_asamp1 <- length(unique(dat_asamp1$deidnum))

asamp1_CR <- dat_asamp1[dat_asamp1$Treatment == "Caloric Restriction",]
n_sub_asamp1_CR <- length(unique(asamp1_CR$deidnum))

asamp1_AL <- dat_asamp1[dat_asamp1$Treatment == "Ad Libitum",]
n_sub_asamp1_AL <- length(unique(asamp1_AL$deidnum))

## N obs
dat_CR <- dat[dat$Treatment == "Caloric Restriction",]
dat_AL <- dat[dat$Treatment == "Ad Libitum",]

n_obs_CR <- length(unique(dat_CR$barcode))
n_obs_AL <- length(unique(dat_AL$barcode))


################################
# Core modeling (shared by analysis & plotting)
################################
# Build & fit model for a given clock; return model + model_data actually used
fit_clock <- function(clock_name, include_cells = FALSE, data = analysis_df) {
  outcome_var  <- paste0("d_", clock_name, "_bc")
  baseline_var <- paste0(clock_name, "_bc_baseline")
  
  # Fixed effects: Visit.ID * Treatment + covariates + baseline value
  base_terms <- c("Visit.ID * Treatment",
                  "deidsite", "bmistrat", "race3_std", "sex", "cbage",
                  baseline_var)
  if (include_cells) {
    base_terms <- c(base_terms,
                    "d_epiccd4t_br","d_epiccd8t_br","d_epicnk_br",
                    "d_epicbcell_br","d_epicmono_br","d_epicneu_br")
  }
  fixed_part <- paste(base_terms, collapse = " + ")
  
  # Random effects: independent intercept & time slope (stable vs unstructured)
  form <- as.formula(
    paste0(outcome_var, " ~ ", fixed_part, " + (1 + fu_numeric || deidnum)")
  )
  
  model <- lmer(
    form, data = data,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  
  list(model = model, model_data = data)
}

################################
# Robust effect extraction (12m, 12→24m, 0→24m)
################################
extract_effects_cr2 <- function(model, clock_name) {
  mf <- model.frame(model)
  cl <- factor(mf$deidnum)
  
  V <- vcovCR(model, cluster = cl, type = "CR2")
  tt <- coef_test(model, vcov = V, test = "Satterthwaite") %>% as.data.frame()
  
  # Rows as named by factor coding (12m ref; AL ref)
  tr_row <- rownames(tt) == "TreatmentCaloric Restriction"
  ix_row <- rownames(tt) == "Visit.ID24 Month:TreatmentCaloric Restriction"
  
  # safety fallback with regex if names differ slightly
  if (!any(tr_row)) tr_row <- grepl("^Treatment.*Caloric.?Restriction$", rownames(tt))
  if (!any(ix_row)) ix_row <- grepl("^Visit\\.ID24 Month:Treatment.*Caloric.?Restriction$", rownames(tt))
  
  stopifnot(sum(tr_row) == 1L, sum(ix_row) == 1L)
  
  eff_12     <- tt$beta[tr_row]
  se_12      <- tt$SE[tr_row]
  df_12      <- tt$df_Satt[tr_row]
  
  eff_inc    <- tt$beta[ix_row]  # increment 12->24
  se_inc     <- tt$SE[ix_row]
  df_inc     <- tt$df_Satt[ix_row]
  
  ci12_lo <- eff_12 - qt(0.975, df_12) * se_12
  ci12_hi <- eff_12 + qt(0.975, df_12) * se_12
  ciin_lo <- eff_inc - qt(0.975, df_inc) * se_inc
  ciin_hi <- eff_inc + qt(0.975, df_inc) * se_inc
  
  # 0–24 cumulative = sum(main + increment), robust SE via linear combo
  coef_names <- names(fixef(model))
  L <- rep(0, length(coef_names))
  L[which(coef_names == "TreatmentCaloric Restriction")] <- 1
  L[which(coef_names == "Visit.ID24 Month:TreatmentCaloric Restriction")] <- 1
  cum_var <- as.numeric(t(L) %*% V %*% L)
  cum_se  <- sqrt(cum_var)
  df_app  <- min(df_12, df_inc)
  eff_cum <- eff_12 + eff_inc
  ci_c_lo <- eff_cum - qt(0.975, df_app) * cum_se
  ci_c_hi <- eff_cum + qt(0.975, df_app) * cum_se
  p12     <- 2 * pt(abs(eff_12 / se_12), df_12, lower.tail = FALSE)
  pinc    <- 2 * pt(abs(eff_inc / se_inc), df_inc, lower.tail = FALSE)
  pcum    <- 2 * pt(abs(eff_cum / cum_se), df_app, lower.tail = FALSE)
  
  tibble(
    clock = clock_name,
    effect_0_12   = eff_12,    ci_lower_0_12 = ci12_lo,  ci_upper_0_12 = ci12_hi,  p_value_0_12 = p12,
    effect_12_24  = eff_inc,   ci_lower_12_24 = ciin_lo, ci_upper_12_24 = ciin_hi, p_value_12_24 = pinc,
    effect_0_24   = eff_cum,   ci_lower_0_24 = ci_c_lo,  ci_upper_0_24 = ci_c_hi,  p_value_0_24 = pcum
  )
}

################################
# Robust margins (predicted means) for plotting
################################
# Returns CR2-robust emmeans at 12 & 24 months (as observed weighting) + Baseline=0 anchor.
predictions_cr2 <- function(model, clock_name) {
  V <- vcovCR(model, cluster = factor(model.frame(model)$deidnum), type = "CR2")
  
  em <- emmeans(
    model,
    ~ Visit.ID * Treatment,
    weights = "proportional",
    vcov. = V
  )
  
  # Handle version-specific CI column names
  em_df <- as.data.frame(summary(em, infer = c(TRUE, TRUE), level = 0.95, adjust = "none"))
  lower_name <- intersect(names(em_df), c("lower.CL", "asymp.LCL", "LCL"))[1]
  upper_name <- intersect(names(em_df), c("upper.CL", "asymp.UCL", "UCL"))[1]
  if (is.na(lower_name) || is.na(upper_name)) {
    stop("Could not find CI columns in emmeans output. Columns were: ", paste(names(em_df), collapse = ", "))
  }
  em_df <- em_df %>% mutate(ci_lower = .data[[lower_name]], ci_upper = .data[[upper_name]])
  
  preds <- em_df %>%
    mutate(
      Visit.ID  = forcats::fct_relevel(Visit.ID, "12 Month", "24 Month"),
      Treatment = forcats::fct_relevel(Treatment, "Ad Libitum", "Caloric Restriction")
    ) %>%
    select(Visit.ID, Treatment,
           predicted = emmean, se = SE,
           ci_lower, ci_upper) %>%
    # Baseline anchor (0; no CI)
    bind_rows(
      tibble(
        Visit.ID  = factor(c("Baseline","Baseline"),
                           levels = c("Baseline","12 Month","24 Month")),
        Treatment = factor(c("Ad Libitum","Caloric Restriction"),
                           levels = c("Ad Libitum","Caloric Restriction")),
        predicted = 0, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_
      )
    ) %>%
    mutate(
      clock     = clock_name,
      Visit.ID  = factor(Visit.ID, levels = c("Baseline","12 Month","24 Month")),
      Treatment = factor(Treatment, levels = c("Ad Libitum","Caloric Restriction"))
    ) %>%
    arrange(Treatment, Visit.ID)
  
  preds
}

################################
# Plot helper (no CI bars at Baseline)
################################
plot_predicted_trajectories <- function(predictions, clock_name,
                                        colors = c("navy","firebrick3"),
                                        show_ci = TRUE) {
  p <- ggplot(predictions, aes(x = Visit.ID, y = predicted, color = Treatment, group = Treatment)) +
    geom_point(size = 3, shape = 17) +
    geom_line(linewidth = 2) +
    theme_bw() + My_Theme +
    theme(legend.position = "top") +
    scale_color_manual(values = colors, name = "Treatment Group") +
    xlab("Visit") +
    ylab("Change from Baseline (SD units)") +
    ggtitle(paste("Predicted Change in", clock_name))
  
  if (show_ci) {
    p <- p + geom_errorbar(
      data = predictions %>% filter(Visit.ID != "Baseline"),
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.1, linewidth = 1
    )
  }
  p
}

################################
# Define clocks
################################

clock_names <- c(
  "ElasticNet",
  "Ridge",
  "DunedinPACE",
  "PCGrimAge"
)

################################
# Collect results (no cells)
################################

results_no_cells <- list()

for (clk in clock_names) {
  cat("\n", paste(rep("=", 50), collapse = ""), "\nAnalyzing:", clk, "\n", paste(rep("=", 50), collapse = ""), "\n")
  
  fit <- fit_clock(clk, include_cells = FALSE, data = analysis_df)
  model <- fit$model
  
  # Optional: pretty (non-robust) model table for reference
  suppressMessages(print(
    sjPlot::tab_model(
      model,
      terms = c("Visit.ID [24 Month]",
                "Treatment [Caloric Restriction]",
                "Visit.ID [24 Month:TreatmentCaloric Restriction]")
    )
  ))
  
  # Robust effects (12m, 12→24, 0→24)
  eff <- extract_effects_cr2(model, clk)
  results_no_cells[[clk]] <- eff

}

results_no_cells_df <- bind_rows(results_no_cells)
print(results_no_cells_df)

results_no_cells_df <- as.data.frame(results_no_cells_df)

################################
# Collect results (with cells)
################################


results_w_cells <- list()

for (clk in clock_names) {
  cat("\n", paste(rep("=", 50), collapse = ""), "\nAnalyzing:", clk, "\n", paste(rep("=", 50), collapse = ""), "\n")
  
  fit <- fit_clock(clk, include_cells = TRUE, data = analysis_df)
  model <- fit$model
  
  # Optional: pretty (non-robust) model table for reference
  suppressMessages(print(
    sjPlot::tab_model(
      model,
      terms = c("Visit.ID [24 Month]",
                "Treatment [Caloric Restriction]",
                "Visit.ID [24 Month:TreatmentCaloric Restriction]")
    )
  ))
  
  # Robust effects (12m, 12→24, 0→24)
  eff <- extract_effects_cr2(model, clk)
  results_w_cells[[clk]] <- eff
  
}

results_with_cells_df <- bind_rows(results_w_cells)
print(results_with_cells_df)

results_with_cells_df <- as.data.frame(results_with_cells_df)

################################
# Save results
################################

write.csv(results_no_cells_df, "../CALERIE_Output/CALERIE_Clock_Treatment_Effect.csv", row.names = FALSE)
write.csv(results_with_cells_df, "../CALERIE_Output/CALERIE_Clock_Treatment_Effect_with_cells.csv", row.names = FALSE)
