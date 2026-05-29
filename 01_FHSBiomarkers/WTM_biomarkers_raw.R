# ==============================================================================
# Framingham Biomarker Data Processing
# Author: WTM
# Date: 2025
# ==============================================================================

# Setup ------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(janitor)
library(purrr)

# Helper Functions -------------------------------------------------------------

# Load Framingham files with standard format
load_framingham_file <- function(file_c1, file_c2) {
  # Load consent group 1
  data_1 <- read.delim(file_c1, header = FALSE, stringsAsFactors = FALSE)
  data_1 <- data_1[-1, ]                 # Remove metadata row
  colnames(data_1) <- as.character(data_1[1, ])  # First row becomes column names
  data_1 <- data_1[-1, ]                 # Remove header row
  
  # Load consent group 2
  data_2 <- read.delim(file_c2, header = FALSE, stringsAsFactors = FALSE)
  data_2 <- data_2[-1, ]
  colnames(data_2) <- as.character(data_2[1, ])
  data_2 <- data_2[-1, ]
  
  # Combine both consent groups
  bind_rows(data_1, data_2)
}

# Reshape variable to long format by visit
reshape_to_long <- function(df, var_pattern, var_name) {
  df %>%
    select(matches(paste0("subject_id|", var_pattern))) %>%
    pivot_longer(
      cols = -subject_id,
      names_to = "visit",
      values_to = var_name
    ) %>%
    mutate(visit = as.numeric(gsub("[^0-9.-]", "", visit)))
}

# ==============================================================================
# Load Anthropometry Data ------------------------------------------------------

vr_wkthru <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht006027.v4.p14.c1.vr_wkthru_ex09_1_1001s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht006027.v4.p14.c2.vr_wkthru_ex09_1_1001s.HMB-IRB-NPU-MDS.txt"
)

# Clean and filter to original cohort only
vr_wkthru <- vr_wkthru %>%
  clean_names() %>%
  filter(idtype == 1) %>%  # idtype==1 is original cohort
  select(-idtype, -shareid) %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  # Keep only variables of interest, exclude fasting indicators
  select(matches("subject_id|bg|hip|waist|creat")) %>%
  select(!matches("fasting"))

# Reshape each variable from wide to long format -------------------------------
creatinine <- reshape_to_long(vr_wkthru, "creat", "creatinine")
glucose    <- reshape_to_long(vr_wkthru, "bg",    "glucose")

hip   <- reshape_to_long(vr_wkthru, "hip",   "hip")
waist <- reshape_to_long(vr_wkthru, "waist", "waist")

# Combine all anthropometry measures -------------------------------------------
anthropometry <- list(
  creatinine, glucose,
  hip, waist
) %>%
  reduce(full_join, by = c("subject_id", "visit"))

anthropometry <- anthropometry %>%
  filter(visit != 9)

# Convert to numeric and calculate derived measures ----------------------------
anthropometry <- anthropometry %>%
  mutate(across(c(hip, waist), as.numeric)) %>%
  mutate(
    # Waist-to-hip ratio (cardiovascular risk marker)
    hip_waist_ratio = waist / hip
  ) %>%
  select(-hip, -waist)  # Remove individual measurements, keep ratio
# ==============================================================================
# Quick validation:
mean(as.numeric(anthropometry$creatinine), na.rm=TRUE) # 1.027071
mean(as.numeric(anthropometry$glucose), na.rm=TRUE) # 100.9453
mean(anthropometry$hip_waist_ratio, na.rm=TRUE) # 0.923033

# ==============================================================================
# Clinic Biomarker Import  --------------------------------------------

# Visit 1
clinic_1 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000030.v10.p14.c1.ex1_1s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000030.v10.p14.c2.ex1_1s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 1,
    wbc      = as.numeric(a138),
    mcv      = as.numeric(a142),
    bun      = as.numeric(a32),
    uricacid = as.numeric(a33),
    totprot  = as.numeric(a34),
    albumin  = as.numeric(a35)
  )

# Visit 2
clinic_2 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000031.v10.p14.c1.ex1_2s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000031.v10.p14.c2.ex1_2s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 2,
    bun      = as.numeric(b738),
    uricacid = as.numeric(b739),
    totprot  = as.numeric(b741),
    albumin  = as.numeric(b742),
    wbc      = as.numeric(b761),
    hgb      = as.numeric(b763),
    mcv      = as.numeric(b765)
  )

# Visit 3
clinic_3 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000032.v9.p14.c1.ex1_3s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000032.v9.p14.c2.ex1_3s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 3
    # add other non-apo vars here if you need them
  )

clinic <- bind_rows(clinic_1, clinic_2, clinic_3) %>%
  group_by(subject_id, visit) %>%
  summarise(across(everything(), ~ dplyr::coalesce(.[1], NA_real_)), .groups = "drop") %>%
  select(-hgb)

# ==============================================================================
# Quick validation:
mean(as.numeric(clinic$wbc), na.rm=TRUE) # 63.45021
mean(as.numeric(clinic$mcv), na.rm=TRUE) # 91.26976
mean(as.numeric(clinic$bun), na.rm=TRUE) # 15.15735
mean(as.numeric(clinic$uricacid), na.rm=TRUE) # 54.21726
mean(as.numeric(clinic$totprot), na.rm=TRUE) # 72.45566
mean(as.numeric(clinic$albumin), na.rm=TRUE) # 45.71692
mean(as.numeric(clinic$apolipo_ratio), na.rm=TRUE) # 0.6527008

# ==============================================================================
# CRP:

# Visit 2
crp_2 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000080.v9.p14.c1.crp1_2s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000080.v9.p14.c2.crp1_2s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 2,
    crp = as.numeric(crp)
  )

# Visit 6
crp_6 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000081.v9.p14.c1.crp1_6s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000081.v9.p14.c2.crp1_6s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 6,
    crp = as.numeric(crp)
  )

# Visit 7 (only c1 file in your legacy snippet)
crp_7_raw <- read.delim(
  "../Data/RAW/phs000007.v32.pht000082.v8.p13.c1.crp1_7s.HMB-IRB-MDS.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Drop metadata row, then promote row 1 to header, then drop it
crp_7 <- crp_7_raw %>%
  slice(-1) %>%                    # removes "## phv...." row
  { setNames(.[-1, ], as.character(.[1, ])) }  # sets colnames from row 1, drops it

# Now clean + select like the others
crp_7 <- crp_7 %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 7,
    crp = as.numeric(crp)
  )


# Visit 8 (only c1 file in your legacy snippet)
crp_8_raw <- read.delim(
  "../Data/RAW/phs000007.v32.pht002888.v6.p13.c1.l_crp_2008_m_0477s.HMB-IRB-MDS.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

crp_8 <- crp_8_raw %>%
  slice(-1) %>%  # drop "## phv..." row
  { setNames(.[-1, ], as.character(.[1, ])) } %>%  # promote header row, drop it
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 8,
    crp = as.numeric(crp)
  )

head(crp_8)

# Combine (aligned to subject_id + visit)
crp <- bind_rows(crp_2, crp_6, crp_7, crp_8) %>%
  group_by(subject_id, visit) %>%
  summarise(crp = dplyr::coalesce(crp[1], NA_real_), .groups = "drop")

# ==============================================================================
# Quick validation:

mean(crp$crp, na.rm=TRUE) # 3.53
quantile(crp$crp, probs = c(0.5, 0.75, 0.9, 0.95, 0.99), na.rm = TRUE)
# Outliers making it very high


# ==============================================================================
# Insulin  ----------------------------------------------------------------------
#  - Visit 5, 7: dbtlab table (pht010725...) filtered by EXAM
#  - Visit 8: 2008 insulin table (pht003901...)
# ==============================================================================

# --- Visits 5 & 7 (dbtlab) ---
dbtlab <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht010725.v2.p14.c1.l_dbtlab_ex07_1b_1237s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht010725.v2.p14.c2.l_dbtlab_ex07_1b_1237s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  select(-shareid, -idtype)

insulin_5_7 <- dbtlab %>%
  mutate(exam = as.numeric(exam)) %>%
  filter(exam %in% c(5, 7)) %>%
  transmute(
    subject_id,
    visit = exam,
    insulin = dplyr::case_when(
      exam == 5 ~ as.numeric(insln_pf_dpc),
      exam == 7 ~ as.numeric(insln_pf_lnc),
      TRUE ~ NA_real_
    )
  )

# --- Visit 8 (2008 insulin table) ---
insulin_8_raw <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht003901.v5.p14.c1.l_insulin_2008_m_0704s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht003901.v5.p14.c2.l_insulin_2008_m_0704s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  select(-shareid, -idtype) %>%
  mutate(visit = 8)

# Inspect once to confirm the insulin column name:
# names(insulin_8_raw)

# After you confirm (very likely `insulin`):
insulin_8 <- insulin_8_raw %>%
  transmute(subject_id, visit, insulin = as.numeric(insulin))

# --- Combine all insulin visits ---
insulin <- bind_rows(insulin_5_7, insulin_8) %>%
  group_by(subject_id, visit) %>%
  summarise(insulin = dplyr::coalesce(insulin[1], NA_real_), .groups = "drop")

# ==============================================================================
# Quick validation
mean(insulin_5_7$insulin, na.rm=TRUE)
mean(insulin_8$insulin, na.rm=TRUE)

quantile(insulin_5_7$insulin, probs = c(0.1, 0.5, 0.9), na.rm = TRUE) / 10
quantile(insulin_8$insulin, probs = c(0.1, 0.5, 0.9), na.rm = TRUE) / 10

# ==============================================================================
# Fibrinogen  -------------------------------------------------------------------
# Sources:
#  - Visit 5: fib1_5s
#  - Visit 6: fib1_6s
#  - Visit 7: fibrin1_7s
# ==============================================================================

# --- Visit 5 ---
fibrin_5 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000084.v9.p14.c1.fib1_5s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000084.v9.p14.c2.fib1_5s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  select(-shareid, -idtype) %>%
  transmute(
    subject_id,
    visit = 5,
    fibrinogen = as.numeric(fibrinog)
  )

# --- Visit 6 ---
fibrin_6 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000085.v9.p14.c1.fib1_6s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000085.v9.p14.c2.fib1_6s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  select(-shareid, -idtype) %>%
  transmute(
    subject_id,
    visit = 6,
    fibrinogen = as.numeric(fibrin)
  )

# --- Visit 7 (single-file format) ---
fibrin_7_raw <- read.delim(
  "../Data/RAW/phs000007.v32.pht000087.v8.p13.c1.fibrin1_7s.HMB-IRB-MDS.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

fibrin_7 <- fibrin_7_raw %>%
  slice(-1) %>%  # drop phv metadata row
  { setNames(.[-1, ], as.character(.[1, ])) } %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 7,
    fibrinogen = as.numeric(fibrngen)
  )

# --- Combine ---
fibrinogen <- bind_rows(fibrin_5, fibrin_6, fibrin_7) %>%
  group_by(subject_id, visit) %>%
  summarise(fibrinogen = dplyr::coalesce(fibrinogen[1], NA_real_), .groups = "drop")

# ==============================================================================
# Quick validation

mean(fibrin_5$fibrinogen, na.rm=TRUE)
mean(fibrin_6$fibrinogen, na.rm=TRUE)
mean(fibrin_7$fibrinogen, na.rm=TRUE)

# ==============================================================================
# Pulmonary (PFTs)  -------------------------------------------------------------
# Final variables:
#   - fev1
#   - fvc
#   - fev1_fvc
#   - fef25_75
# ==============================================================================

pulmonary_standardize <- function(df, visit,
                                  fev1_col,
                                  fvc_col,
                                  ratio_col,
                                  fef_col) {
  
  # clean names of requested columns to match clean_names(df)
  fev1_col  <- janitor::make_clean_names(fev1_col)
  fvc_col   <- janitor::make_clean_names(fvc_col)
  ratio_col <- janitor::make_clean_names(ratio_col)
  fef_col   <- janitor::make_clean_names(fef_col)
  
  df %>%
    clean_names() %>%
    rename(subject_id = db_ga_p_subject_id) %>%
    select(-shareid, -idtype) %>%
    transmute(
      subject_id,
      visit = visit,
      fev1     = as.numeric(.data[[fev1_col]]),
      fvc      = as.numeric(.data[[fvc_col]]),
      fev1_fvc = as.numeric(.data[[ratio_col]]),
      fef25 = as.numeric(.data[[fef_col]])
    )
}

# --- Visit 3 ---
pulmonary_3 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000105.v9.p14.c1.pft1_3s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000105.v9.p14.c2.pft1_3s.HMB-IRB-NPU-MDS.txt"
) %>%
  pulmonary_standardize(
    visit     = 3,
    fev1_col  = "FV1_3_1",
    fvc_col   = "FVC_3_1",
    ratio_col = "RAT_3_1",
    fef_col   = "FF1_3_1"
  )

names(pulmonary_3)

# --- Visit 5 ---
pulmonary_5 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000101.v9.p14.c1.pft1_5s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000101.v9.p14.c2.pft1_5s.HMB-IRB-NPU-MDS.txt"
) %>%
  pulmonary_standardize(
    visit     = 5,
    fev1_col  = "FV1_5_1",
    fvc_col   = "FVC_5_1",
    ratio_col = "RAT_5_1",
    fef_col   = "FF1_5_1"
  )

# --- Visit 6 ---
pulmonary_6 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000104.v9.p14.c1.pft1_6s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000104.v9.p14.c2.pft1_6s.HMB-IRB-NPU-MDS.txt"
) %>%
  pulmonary_standardize(
    visit     = 6,
    fev1_col  = "FV1_6_1",
    fvc_col   = "FVC_6_1",
    ratio_col = "RAT_6_1",
    fef_col   = "FF1_6_1"
  )

# --- Visit 7 ---
pulmonary_7 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000106.v9.p14.c1.pft1_7s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000106.v9.p14.c2.pft1_7s.HMB-IRB-NPU-MDS.txt"
) %>%
  pulmonary_standardize(
    visit     = 7,
    fev1_col  = "FV1_7_1",
    fvc_col   = "FVC_7_1",
    ratio_col = "RAT_7_1",
    fef_col   = "FF1_7_1"
  )

# --- Visit 8 ---
pulmonary_8 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000832.v8.p14.c1.pft1_8s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000832.v8.p14.c2.pft1_8s.HMB-IRB-NPU-MDS.txt"
) %>%
  pulmonary_standardize(
    visit     = 8,
    fev1_col  = "fv1_8_1",
    fvc_col   = "fvc_8_1",
    ratio_col = "rbb_8_1",
    fef_col   = "ff1_8_1"
  )

# --- Combine ---
pulmonary <- bind_rows(
  pulmonary_3,
  pulmonary_5,
  pulmonary_6,
  pulmonary_7,
  pulmonary_8
) %>%
  group_by(subject_id, visit) %>%
  summarise(
    across(everything(), ~ dplyr::coalesce(.[1], NA_real_)),
    .groups = "drop"
  )

# ==============================================================================
# Quick validation  -------------------------------------------------------------
mean(pulmonary$fev1) # 2.8
mean(pulmonary$fvc) # 3.8
mean(pulmonary$fev1_fvc) # 0.73
mean(pulmonary$fef25, na.rm=TRUE) # 6.2

# ==============================================================================
# Interleukin-6 (IL-6)  ----------------------------------------------------------
# Visits: 7, 8
# ==============================================================================

# --- Visit 7 ---
il6_7 <- read.delim(
  "../Data/RAW/phs000007.v32.pht000161.v8.p13.c1.il61_7s.HMB-IRB-MDS.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

il6_7 <- il6_7 %>%
  slice(-1) %>%                          # drop phv metadata row
  { setNames(.[-1, ], as.character(.[1, ])) } %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 7,
    il6 = as.numeric(il6)
  )

# --- Visit 8 ---
il6_8 <- read.delim(
  "../Data/RAW/phs000007.v32.pht002891.v6.p13.c1.l_il6_2008_m_0433s.HMB-IRB-MDS.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

il6_8 <- il6_8 %>%
  slice(-1) %>%
  { setNames(.[-1, ], as.character(.[1, ])) } %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 8,
    il6 = as.numeric(il6)
  )

# --- Combine ---
il6 <- bind_rows(il6_7, il6_8) %>%
  group_by(subject_id, visit) %>%
  summarise(il6 = dplyr::coalesce(il6[1], NA_real_), .groups = "drop")

# ==============================================================================
# Quick validation  -------------------------------------------------------------
mean(il6$il6)
summary(il6$il6)


# ==============================================================================
# Apolipoproteins (Exam 5)  ------------------------------------------------------
# Source: lipids1_5s
# ==============================================================================

# Apolipoproteins (Visits 3 + 5) -----------------------------------------------

apolipo_3 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000032.v9.p14.c1.ex1_3s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000032.v9.p14.c2.ex1_3s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  transmute(
    subject_id = db_ga_p_subject_id,
    visit = 3,
    apoa1 = as.numeric(c443),
    apob  = as.numeric(c444)
  ) %>%
  mutate(
    apoa1 = ifelse(apoa1 <= 0, NA_real_, apoa1),
    apob  = ifelse(apob  <= 0, NA_real_, apob),
    apolipo_ratio = apob / apoa1
  ) %>%
  select(subject_id, visit, apolipo_ratio)

apolipo_5 <- load_framingham_file(
  "../Data/RAW/phs000007.v33.pht000205.v9.p14.c1.lipids1_5s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000205.v9.p14.c2.lipids1_5s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  select(-shareid, -idtype) %>%
  transmute(
    subject_id,
    visit = 5,
    apoa1 = as.numeric(apoa14),
    apob  = as.numeric(apob4)
  ) %>%
  mutate(
    apoa1 = na_if(apoa1, 0),
    apob  = na_if(apob, 0),
    apoa1 = ifelse(apoa1 < 0, NA_real_, apoa1),
    apob  = ifelse(apob  < 0, NA_real_, apob),
    apolipo_ratio = apob / apoa1
  ) %>%
  select(subject_id, visit, apolipo_ratio)

apolipo <- bind_rows(apolipo_3, apolipo_5) %>%
  group_by(subject_id, visit) %>%
  summarise(apolipo_ratio = dplyr::coalesce(apolipo_ratio[1], NA_real_), .groups = "drop")


# ==============================================================================
# Homocysteine  ---------------------------------------------------------------
# Visits: 5, 6, 7, 8
# Source: laba1_* / fhslab1_* chemistry panels
# ==============================================================================

extract_hcy <- function(file_c1, file_c2, visit, hcy_col) {
  load_framingham_file(file_c1, file_c2) %>%
    clean_names() %>%
    rename(subject_id = db_ga_p_subject_id) %>%
    select(-shareid, -idtype) %>%
    transmute(
      subject_id,
      visit = visit,
      homocysteine = as.numeric(.data[[janitor::make_clean_names(hcy_col)]])
    )
}

# --- Visit 5 ---
homocysteine_5 <- extract_hcy(
  "../Data/RAW/phs000007.v33.pht000202.v9.p14.c1.laba1_5s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000202.v9.p14.c2.laba1_5s.HMB-IRB-NPU-MDS.txt",
  visit = 5,
  hcy_col = "HCYST"
)

# --- Visit 6 ---
homocysteine_6 <- extract_hcy(
  "../Data/RAW/phs000007.v33.pht000203.v9.p14.c1.laba1_6s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000203.v9.p14.c2.laba1_6s.HMB-IRB-NPU-MDS.txt",
  visit = 6,
  hcy_col = "HCYST1"
)

# --- Visit 7 ---
homocysteine_7 <- extract_hcy(
  "../Data/RAW/phs000007.v33.pht000204.v9.p14.c1.laba1_7s.HMB-IRB-MDS.txt",
  "../Data/RAW/phs000007.v33.pht000204.v9.p14.c2.laba1_7s.HMB-IRB-NPU-MDS.txt",
  visit = 7,
  hcy_col = "LB_HCYST"
)

# --- Combine ---
homocysteine <- bind_rows(
  homocysteine_5,
  homocysteine_6,
  homocysteine_7
) %>%
  group_by(subject_id, visit) %>%
  summarise(
    homocysteine = dplyr::coalesce(homocysteine[1], NA_real_),
    .groups = "drop"
  )

# ==============================================================================
# Quick validation  -------------------------------------------------------------
mean(homocysteine$homocysteine, na.rm=TRUE)
summary(homocysteine$homocysteine)


# ==============================================================================
# Sanity Check  ---------------------------------------------------------------
# ==============================================================================
check_keys <- function(df, name) {
  stopifnot(!any(is.na(df$subject_id)))
  stopifnot(!any(is.na(df$visit)))
  message(name, ": ", nrow(df), " rows")
}

check_keys(anthropometry, "anthropometry")
check_keys(clinic, "clinic")
check_keys(crp, "crp")
check_keys(insulin, "insulin")
check_keys(fibrinogen, "fibrinogen")
check_keys(pulmonary, "pulmonary")
check_keys(il6, "il6")
check_keys(homocysteine, "homocysteine")

# ==============================================================================
# Merge  ---------------------------------------------------------------
# ==============================================================================

biomarker <- list(
  anthropometry,
  clinic,
  crp,
  insulin,
  fibrinogen,
  pulmonary,
  il6,
  homocysteine,
  apolipo
) %>%
  reduce(full_join, by = c("subject_id", "visit"))

biomarker <- biomarker %>%
  arrange(subject_id, visit) %>%
  mutate(across(-c(subject_id, visit), as.numeric))

biomarker <- biomarker %>%
  select(
    subject_id,
    visit,
    
    # Metabolic / renal
    glucose,
    insulin,
    bun,
    creatinine,
    uricacid,
    
    # Proteins / hematologic
    albumin,
    totprot,
    wbc,
    mcv,
    
    # Lipids
    apolipo_ratio,
    
    # Inflammation
    crp,
    il6,
    fibrinogen,
    homocysteine,
    
    # Pulmonary
    fev1,
    fvc,
    fev1_fvc,
    fef25,
    
    # Anthropometry
    hip_waist_ratio
  )

# ==============================================================================
# Global Validation  ---------------------------------------------------------------
# ==============================================================================

summary(biomarker)

table(biomarker$visit)

colMeans(is.na(biomarker[ , -c(1,2)]))

# List of biomarker columns (everything except keys)
biomarker_vars <- setdiff(names(biomarker), c("subject_id", "visit"))

# Long format: one row per biomarker per subject/visit
biomarker_long <- biomarker %>%
  pivot_longer(
    cols = all_of(biomarker_vars),
    names_to = "biomarker",
    values_to = "value"
  )

# Presence table: count non-missing per visit
biomarker_visit_table <- biomarker_long %>%
  group_by(biomarker, visit) %>%
  summarise(
    n_non_missing = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  arrange(biomarker, visit)

biomarker_visit_summary <- biomarker_visit_table %>%
  filter(n_non_missing > 0) %>%
  group_by(biomarker) %>%
  summarise(
    visits = paste(sort(unique(visit)), collapse = ", "),
    .groups = "drop"
  )

biomarker_visit_summary


# ==============================================================================
# Load in additional datasets  (WTM_RAW) ---------------------------------------
# ==============================================================================

load_framingham_file <- function(file_c1, file_c2) {
  # dbGaP files have ~10-11 comment lines starting with #, then column names
  data_1 <- read.delim(file_c1, comment.char = "#", stringsAsFactors = FALSE)
  data_2 <- read.delim(file_c2, comment.char = "#", stringsAsFactors = FALSE)
  
  bind_rows(data_1, data_2)
}

homocys1_2s <- load_framingham_file(
  "../Data/WTM_RAW/phs000007.v33.pht000092.v9.p14.c1.homocys1_2s.HMB-IRB-MDS.txt",
  "../Data/WTM_RAW/phs000007.v33.pht000092.v9.p14.c2.homocys1_2s.HMB-IRB-NPU-MDS.txt"
) %>%
  clean_names() %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  select(-shareid, -idtype) %>%
  transmute(
    subject_id,
    visit = 2,
    homocysteine = as.numeric(thc),
  )

lipids1_5s <- load_framingham_file(
  "../Data/WTM_RAW/phs000007.v33.pht000205.v9.p14.c1.lipids1_5s.HMB-IRB-MDS.txt",
  "../Data/WTM_RAW/phs000007.v33.pht000205.v9.p14.c2.lipids1_5s.HMB-IRB-NPU-MDS.txt"
) 

apolipo_ex4 <- lipids1_5s %>%
  clean_names() %>%
  rename(subject_id = db_ga_p_subject_id) %>%
  select(-shareid, -idtype) %>%
  transmute(
    subject_id,
    visit = 4,
    apoa1 = as.numeric(apoa14),
    apob = as.numeric(apob4)
  ) %>%
  mutate(
    apoa1 = na_if(apoa1, 0),
    apob  = na_if(apob, 0),
    apoa1 = ifelse(apoa1 < 0, NA_real_, apoa1),
    apob  = ifelse(apob  < 0, NA_real_, apob),
    apolipo_ratio = apob / apoa1
  ) %>%
  select(subject_id, visit, apolipo_ratio)

# ==============================================================================
# Validate New Datasets Against Existing Data ---------------------------------
# ==============================================================================

# 1. Compare homocysteine scales
cat("\n=== HOMOCYSTEINE VALIDATION ===\n")
cat("Existing data (visits 5-7):\n")
print(summary(homocysteine$homocysteine))

cat("\nNew visit 2 data:\n")
print(summary(homocys1_2s$homocysteine))

# 3. Compare apolipo_ratio scales
cat("\n=== APOLIPOPROTEIN RATIO VALIDATION ===\n")
cat("Existing data (visits 3, 5):\n")
print(summary(apolipo$apolipo_ratio))

cat("\nNew visit 4 data:\n")
print(summary(apolipo_ex4$apolipo_ratio))

# ==============================================================================
# Merge New Data ---------------------------------------------------------------
# ==============================================================================
# Strategy: Use coalesce to prefer existing data when there are conflicts
# (assuming existing data is already validated)

# Homocysteine - add visit 2
homocysteine_updated <- bind_rows(
  homocysteine,
  homocys1_2s %>% mutate(subject_id = as.character(subject_id))
) %>%
  group_by(subject_id, visit) %>%
  summarise(homocysteine = coalesce(homocysteine[1], NA_real_), .groups = "drop")


# Apolipoprotein ratio - add visit 4
apolipo_updated <- bind_rows(
  apolipo,
  apolipo_ex4 %>% mutate(subject_id = as.character(subject_id))
) %>%
  group_by(subject_id, visit) %>%
  summarise(apolipo_ratio = coalesce(apolipo_ratio[1], NA_real_), .groups = "drop")

# Now re-merge everything
biomarker_updated <- list(
  anthropometry,
  clinic,
  crp,
  insulin,
  fibrinogen,
  pulmonary,
  il6,
  homocysteine_updated,
  apolipo_updated
) %>%
  reduce(full_join, by = c("subject_id", "visit"))

# Apply same final formatting as before
biomarker_updated <- biomarker_updated %>%
  arrange(subject_id, visit) %>%
  mutate(across(-c(subject_id, visit), as.numeric)) %>%
  select(
    subject_id, visit,
    glucose, insulin, bun, creatinine, uricacid,
    albumin, totprot, wbc, mcv,
    apolipo_ratio,
    crp, il6, fibrinogen, homocysteine,
    fev1, fvc, fev1_fvc, fef25,
    hip_waist_ratio
  )

# ==============================================================================
# Validation -------------------------------------------------------------------
# ==============================================================================

cat("\n=== UPDATED COVERAGE ===\n")
biomarker_visit_summary_updated <- biomarker_updated %>%
  pivot_longer(cols = -c(subject_id, visit), names_to = "biomarker", values_to = "value") %>%
  group_by(biomarker, visit) %>%
  summarise(n_non_missing = sum(!is.na(value)), .groups = "drop") %>%
  filter(n_non_missing > 0) %>%
  group_by(biomarker) %>%
  summarise(visits = paste(sort(unique(visit)), collapse = ", "), .groups = "drop")

print(biomarker_visit_summary_updated)


# ==============================================================================
# Output  ---------------------------------------------------------------
# ==============================================================================
write.csv(biomarker_updated, "../Data/WTM_biomarkers_raw.csv")

