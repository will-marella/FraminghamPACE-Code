# Load packages
library(data.table)
library(dplyr)

# Load raw data files
raw1 <- fread("../CARDIA_Data/Raw/phs000285.v3.pht001559.v2.p2.c2.A4F01.HMB-IRB-NPU.txt.gz",
              fill = TRUE,
              skip=10)

raw2 <- fread("../CARDIA_Data/Raw/phs000285.v3.pht001559.v2.p2.c1.A4F01.HMB-IRB.txt.gz",
              fill = TRUE,
              skip=10)

# Get important info
raw1_df <- as.data.frame(raw1)
raw2_df <- as.data.frame(raw2)

race1_df <- raw1_df %>%
  select(dbGaP_Subject_ID, Individual_ID, A01RACE1) %>%
  rename(Race = A01RACE1) %>%
  mutate(Race_char = case_when(
    Race == 5 ~ "White, not Hispanic",
    Race == 4 ~ "Black, not Hispanic", 
    Race == 3 ~ "Hispanic",
    TRUE ~ NA_character_  # for any other values
  ))

race2_df <- raw2_df %>%
  select(dbGaP_Subject_ID, Individual_ID, A01RACE1) %>%
  rename(Race = A01RACE1) %>%
  mutate(Race_char = case_when(
    Race == 5 ~ "White, not Hispanic",
    Race == 4 ~ "Black, not Hispanic", 
    Race == 3 ~ "Hispanic",
    TRUE ~ NA_character_  # for any other values
  ))

# Merge into one dataframe
race_df <- rbind(
  race1_df,
  race2_df)

# Export
write.csv(race_df, "../CARDIA_Data/Raw/CARDIA_race_data.csv")

