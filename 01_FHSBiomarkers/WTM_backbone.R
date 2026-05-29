### Backbone Data Set Up
# Set Up  --------------------------------------------
library(dplyr)
library(tidyr)
library(janitor)
library(haven)

# Import Backbone Data --------------------------------------------
ages_dates <- read_dta("../Data/RAW/Age&Dates_AllCohorts_v33.dta") %>%
  clean_names() %>%
  rename(subject_id = dbgap_subject_id,
         id_type = idtype) %>%
  filter(id_type == 1) %>%  # Original cohort only
  select(-id_type, -shareid)

# Extract sex
sex <- ages_dates %>%
  distinct(subject_id, sex) %>%
  mutate(sex = ifelse(sex == 1, "Male", "Female"))

# Extract age (visits 1-8 only)
age <- ages_dates %>%
  select(subject_id, starts_with("age")) %>%
  select(subject_id, age1:age8) %>%  # Only visits 1-8
  pivot_longer(cols = -subject_id, names_to = "visit", values_to = "age") %>%
  mutate(visit = gsub("[^0-9.-]", "", visit)) %>%
  filter(!is.na(age))

# Extract attendance (visits 1-8 only)
att <- ages_dates %>%
  select(subject_id, starts_with("att")) %>%
  select(subject_id, att1:att8) %>%  # Only visits 1-8
  pivot_longer(cols = -subject_id, names_to = "visit", values_to = "att") %>%
  mutate(visit = gsub("[^0-9.-]", "", visit)) %>%
  filter(!is.na(att))

# Extract date (visits 2-8 only)
date <- ages_dates %>%
  select(subject_id, starts_with("date")) %>%
  select(subject_id, date2:date8) %>%  # Only visits 2-8
  pivot_longer(cols = -subject_id, names_to = "visit", values_to = "time") %>%
  mutate(visit = gsub("[^0-9.-]", "", visit)) %>%
  filter(!is.na(time))

# Baseline age
baseline_age <- age %>%
  filter(visit == "1") %>%
  select(subject_id, age) %>%
  rename(baseline_age = age)

# Merge all
backbone <- att %>%
  full_join(age, by = c("subject_id", "visit")) %>%
  full_join(date, by = c("subject_id", "visit")) %>%
  left_join(sex, by = "subject_id") %>%
  left_join(baseline_age, by = "subject_id") %>%
  filter(att != 0) %>%  # Remove non-attendees
  mutate(time = ifelse(visit == "1", 0, time))  # Set baseline time to 0

# Data Export --------------------------------------------
write.csv(backbone, file="../Data/WTM_backbone.csv", row.names = FALSE)