# Imports
library(dplyr)
library(tidyr)
library(haven)

# Utility to track data loss through processing pipeline.
# Returns counts of non-missing, finite values for each biomarker at a given stage.
# Use: counts0 <- count_nonmissing(df, bio_cols, "0_post_transform")
count_nonmissing <- function(d, cols, label) {
  tibble(
    stage = label,
    biomarker = cols,
    n = sapply(cols, \(v) sum(!is.na(d[[v]]) & is.finite(d[[v]])))
  )
}

# ----------------------------
# Import
# ----------------------------

# Define biomarker columns in the dataset
biomarker_cols <- c("glucose", "insulin", "bun", "creatinine",
                    "uricacid", "albumin", "totprot", "wbc",
                    "mcv", "apolipo_ratio", "crp", "il6", "fibrinogen",
                    "homocysteine", "fev1", "fvc", "fev1_fvc",
                    "fef25", "hip_waist_ratio")

# Raw biomarker dataset
# Need subject/visit/biomarker_cols
bio_raw <- read.csv("../Data/WTM_biomarkers_raw.csv") %>%
  select(subject_id, visit, all_of(biomarker_cols)) 

# Backbone dataset for ages and dates
backbone <- read.csv("../Data/WTM_backbone.csv") %>%
  mutate(
    c_age = as.numeric(scale(age, scale = FALSE)),
    c_age2 = c_age^2
  )

df <- backbone %>%
  inner_join(bio_raw, by = c("subject_id", "visit")) %>%
  mutate(
    visit = factor(visit),
    sex = as.factor(sex)
  )

stopifnot(all(biomarker_cols %in% names(df)))

df <- df %>%
  mutate(across(all_of(biomarker_cols), ~ suppressWarnings(as.numeric(.x))))

df <- df %>%
  mutate(
    age = if_else(is.na(age) & !is.na(baseline_age) & !is.na(time),
                  baseline_age + time / 365.25,  # Add years since baseline
                  age)
  )

bio_cols_raw <- biomarker_cols

# ----------------------------
# Choose transforms (ONLY these will be used downstream)
#   - if a var is in log_map, we will process var_ln and NOT var
# ----------------------------

# Define biomarkers to be logged
# (Be sure to justify which variables were logged)
# log: for naturally positive quantities (ratios)
# log1p: for biomarkers that may include zeros/noise
# Use log1p when in doubt
log_map <- c(
  crp = "log1p",
  insulin = "log1p",
  glucose = "log1p",
  il6 = "log1p",
  homocysteine = "log1p",
  bun = "log1p",
  creatinine = "log1p",
  wbc = "log1p",
  apolipo_ratio = "log"   # ratio > 0
)

log_map <- log_map[names(log_map) %in% names(df)]  # only what exists


for (v in names(log_map)) {
  fn  <- log_map[[v]]
  new <- paste0(v, "_ln")
  
  df[[new]] <- if (fn == "log1p") {
    x <- df[[v]]
    x[x < 0] <- NA_real_
    log1p(x)
  } else {
    x <- df[[v]]
    x[x <= 0] <- NA_real_
    log(x)
  }
}

bio_cols <- setdiff(bio_cols_raw, names(log_map))
bio_cols <- unique(c(bio_cols, paste0(names(log_map), "_ln")))
bio_cols <- bio_cols[bio_cols %in% names(df)]

# Check sum of biomarker counts post log
counts0 <- count_nonmissing(df, bio_cols, "0_post_transform") 
sum(counts0$n)

# ----------------------------
# Remove 5 SD outliers within sex
# ----------------------------

df <- df %>%
  group_by(sex) %>%
  mutate(across(all_of(bio_cols), ~{
    m <- mean(.x, na.rm = TRUE)
    s <- sd(.x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(.x)
    ifelse(.x < m - 5*s | .x > m + 5*s, NA_real_, .x)
  })) %>%
  ungroup()

# Check counts post outlier removal
counts1 <- count_nonmissing(df, bio_cols, "1_post_outliers") 
sum(counts1$n)

# ----------------------------
# Residualize visit effects while controlling for age and sex
# For each biomarker, fit a linear model: biomarker ~ visit + c_age + c_age2 + sex + interactions
# Then subtract the visit coefficient from observations at that visit.
# This removes systematic shifts between visits while holding age/sex constant.
# ----------------------------

# Define the visits to adjust 
# (All visits besides the first visit they are measured at)
# glucose is measured at visits 1-8, so we list 2-8.
visit_adjust_map <- list(
  albumin          = c(2, 3),
  apolipo_ratio_ln = c(4, 5),
  bun_ln           = c(2),
  creatinine_ln    = c(2, 3, 5,6,7,8),
  crp_ln           = c(2, 3, 6,7,8),
  glucose_ln       = c(2,3,4,5,6,7,8),
  homocysteine_ln  = c(5, 6,7),
  il6_ln           = c(7, 8),
  insulin_ln       = c(7,8),
  mcv              = c(2),
  totprot          = c(2),
  uricacid         = c(2),
  wbc_ln           = c(2),
  hip_waist_ratio  = c(3, 4, 5,6,7),
  fev1             = c(5,6,7,8),
  fvc              = c(5,6,7,8),
  fibrinogen       = c(5,6,7)
)

# safety: only keep ones that exist in your bio_cols
visit_adjust_map <- visit_adjust_map[names(visit_adjust_map) %in% bio_cols]

# ensure consistent reference visit (even if it has no visit 1 measure)
df <- df %>% mutate(visit = relevel(factor(visit), ref = "1"))

for (y in bio_cols) {
  yy <- df[[y]]
  y_adj <- yy  # start with original values
  
  if (y %in% names(visit_adjust_map)) {
    # Fit model to get coefficients
    f <- as.formula(paste0(y, " ~ visit + c_age + c_age2 + sex + c_age:sex + c_age2:sex"))
    fit <- lm(f, data = df, na.action = na.exclude)
    coefs <- coef(fit)
    
    # Get which visits to adjust
    vset <- as.character(visit_adjust_map[[y]])
    
    # Subtract coefficient for each specified visit
    for (i in seq_along(vset)) {
      v <- vset[i]
      coef_name <- paste0("visit", v)
      if (!is.na(coefs[coef_name])) {
        y_adj[df$visit == v & !is.na(yy)] <- yy[df$visit == v & !is.na(yy)] - coefs[coef_name]
      }
    }
  }
  
  df[[paste0(y, "_adj")]] <- y_adj
}

# Nice check:
# how many values were changed per biomarker
changed_counts <- sapply(bio_cols, function(y) {
  a <- df[[paste0(y, "_adj")]]
  o <- df[[y]]
  sum(!is.na(a) & !is.na(o) & a != o)
})
changed_counts
# looks exactly right

# Check for any lost samples
counts2 <- count_nonmissing(df, paste0(bio_cols, "_adj"), "2_post_visitadj")
sum(counts2$n)

# ----------------------------
# Standardize within sex using age 40–50 reference window + debugging
# ----------------------------
adj_cols <- paste0(bio_cols, "_adj")

# Count before
before_counts <- sapply(adj_cols, function(v) sum(!is.na(df[[v]])))
cat("Total observations before standardization:", sum(before_counts), "\n")

# Compute reference stats (age 40-50)
ref <- df %>%
  filter(between(age, 40, 50)) %>% # testing
  group_by(sex) %>%
  summarise(
    across(
      all_of(adj_cols),
      list(mu = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# Apply standardization
for (v in adj_cols) {
  mu_col <- paste0(v, "_mu")
  sd_col <- paste0(v, "_sd")
  sc_col <- paste0(v, "_sc")
  
  df <- df %>%
    left_join(ref %>% select(sex, all_of(c(mu_col, sd_col))), by = "sex") %>%
    mutate(
      !!sc_col := (!!sym(v) - !!sym(mu_col)) / !!sym(sd_col)
    ) %>%
    select(-all_of(c(mu_col, sd_col)))
}


# Count after
after_counts <- sapply(paste0(adj_cols, "_sc"), function(v) sum(!is.na(df[[v]])))
cat("Total observations after standardization:", sum(after_counts), "\n")
cat("Lost:", sum(before_counts) - sum(after_counts), "observations\n")
cat("Percent lost:", round(100 * (sum(before_counts) - sum(after_counts)) / sum(before_counts), 2), "%\n")


# Some quick validation
df %>%
  filter(between(age, 40, 50)) %>%
  group_by(sex) %>%
  summarise(across(ends_with("_sc"), ~mean(.x, na.rm = TRUE))) %>%
  print()

df %>%
  filter(between(age, 40, 50)) %>%
  group_by(sex) %>%
  summarise(across(ends_with("_sc"), ~sd(.x, na.rm = TRUE))) %>%
  print()

counts3 <- count_nonmissing(df, paste0(bio_cols, "_adj_sc"), "3_post_standardize")
cat("\nFinal count:", sum(counts3$n), "\n")
cat("Rows with missing sex in current df:", sum(is.na(df$sex)), "\n")

# ----------------------------
# Reverse biomarkers as needed
# (Biomarkers that decrease with aging need to be reversed)
# ----------------------------

df <- df %>%
  mutate(
    albumin_adj_sc      = -albumin_adj_sc,
    fev1_adj_sc         = -fev1_adj_sc,
    fvc_adj_sc          = -fvc_adj_sc,
    fev1_fvc_adj_sc     = -fev1_fvc_adj_sc,
    fef25_adj_sc        = -fef25_adj_sc,
    totprot_adj_sc      = -totprot_adj_sc
  )

# ----------------------------
# Validate distributions
# ----------------------------

hist(df$uricacid_adj_sc) # Good
hist(df$albumin_adj_sc) # Decent
hist(df$totprot_adj_sc) # Good
hist(df$mcv_adj_sc) # Decent
hist(df$fibrinogen_adj_sc) # Good
hist(df$fev1_adj_sc) # Good
hist(df$fvc_adj_sc) # Good
hist(df$fev1_fvc_adj_sc) # Pretty skewed
hist(df$fef25_adj_sc) # Good
hist(df$hip_waist_ratio_adj_sc) # Good
hist(df$crp_ln_adj_sc) # Skewed.
hist(df$insulin_ln_adj_sc) # Decent
hist(df$glucose_ln_adj_sc) # Decent
hist(df$il6_ln_adj_sc) # Manageable skew
hist(df$homocysteine_ln_adj_sc) # Good
hist(df$bun_ln_adj_sc) # Good
hist(df$creatinine_ln_adj_sc) # Decent
hist(df$wbc_ln_adj_sc) # Good
hist(df$apolipo_ratio_ln_adj_sc) # Good

# ----------------------------
# Export
# ----------------------------

# Count summary
sum(counts0$n)
sum(counts1$n)
sum(counts2$n)
sum(counts3$n)


write.csv(df, "../Data/WTM_biomarkers_processed.csv", row.names = FALSE)
