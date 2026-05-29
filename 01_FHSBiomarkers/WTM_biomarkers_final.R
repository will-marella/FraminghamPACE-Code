### Finalize biomarker data for input into `nlme`

library(dplyr)
library(tidyr)

# ----------------------------
# Import
# ----------------------------
biomarkers_processed <- read.csv("../Data/WTM_biomarkers_processed.csv")
backbone <- read.csv("../Data/WTM_backbone.csv")

# ----------------------------
# Basic hygiene (types) + bring in sex if needed
# ----------------------------
df <- biomarkers_processed %>%
  mutate(
    subject_id = as.character(subject_id),
    visit      = as.integer(visit)
  )

sex_lookup <- backbone %>%
  transmute(subject_id = as.character(subject_id), sex) %>%
  distinct()

df <- df %>%
  left_join(sex_lookup, by = "subject_id", suffix = c("", "_bb")) %>%
  mutate(sex = coalesce(sex, sex_bb)) %>%
  select(-sex_bb) %>%
  relocate(sex, .after = subject_id)

# ----------------------------
# Time scaling (use backbone `time` already in processed file)
#   - assumes `time` is days since baseline/first visit (as in your backbone)
# ----------------------------
stopifnot("time" %in% names(df))

df <- df %>%
  mutate(
    days_since_first_visit = as.numeric(time),
    time_5year = days_since_first_visit / (365.25 * 5),
    time_5year_squared = time_5year^2
  )

# ----------------------------
# Long format
# ----------------------------
bio_cols <- names(df)[grepl("_adj_sc$", names(df))]
bio_cols

biomarkers_final <- df %>%
  pivot_longer(
    cols = all_of(bio_cols),
    names_to = "biomarker",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  # Keep only these columns:
  select(subject_id, visit, biomarker, value, time_5year, time_5year_squared)

unique(biomarkers_final$biomarker)

# ----------------------------
# Drop confounded biomarkers (only if they exist)
# And the unwanted pulmonary variables
# ----------------------------
confounded_biomarkers <- c(
  "bmi_ln_adj", "ldl_adj", "totchol_adj", "dbp_adj",
  "hdl_adj", "sbp_adj", "trig_ln_adj",
  "fvc_adj_sc", "fef25_adj_sc"
)

biomarkers_final <- biomarkers_final %>%
  filter(!biomarker %in% confounded_biomarkers)

# ----------------------------
# Validation
# ----------------------------
scale_check <- biomarkers_final %>%
  group_by(biomarker) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    p05  = quantile(value, 0.05, na.rm = TRUE),
    p95  = quantile(value, 0.95, na.rm = TRUE),
    n    = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  arrange(desc(sd))

scale_check

# ----------------------------
# Save
# ----------------------------
write.csv(biomarkers_final, "../Data/WTM_biomarkers_final.csv", row.names = FALSE)

