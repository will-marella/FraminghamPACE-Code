# prepare_hrs_outcomes_and_clocks.R
library(dplyr)
library(haven)
library(readr)
library(lubridate)


raw_outcomes_df <- read_dta("../HRS_Data/Initial/trk2022tr_r.dta")
pheno_df <- read_csv("../HRS_Data/Initial/HRS_combined_clean_pheno_w_hhid.csv")
PC_Clock_df <- read_csv("../HRS_Data/Initial/HRS_combined_PC_clocks.csv")
Dunedin_df <- read_csv("../HRS_Data/Initial/HRS_combined_DunedinPACE.csv")

# Select the features to take from the large tracker file
filter_features <- c("HHID", "PN", "BIRTHMO", "BIRTHYR",
                     'KNOWNDECEASEDMO', 'KNOWNDECEASEDYR',
                     'LASTALIVEMO', 'LASTALIVEYR',
                     'PIWMONTH', 'PIWYEAR', 'PAGE')

raw_outcomes_df_filtered <- raw_outcomes_df %>%
  select(all_of(filter_features))

# Function to clean up the month/year columns using lubridate
clean_month_year <- function(month, year, min_year = 1900, max_year = 2024) {
  # Clean months and years
  clean_month <- ifelse(month %in% 1:12, month, NA)
  clean_year <- ifelse(year >= min_year & year <= max_year, year, NA)
  
  # Create dates
  make_date(year = clean_year, month = clean_month, day = 1)
}

# Clean up birth date, death date, and lastalive date
raw_outcomes_df_filtered$birth_date <- clean_month_year(
  raw_outcomes_df_filtered$BIRTHMO,
  raw_outcomes_df_filtered$BIRTHYR,
  min_year = 1900
)

raw_outcomes_df_filtered$death_date <- clean_month_year(
  raw_outcomes_df_filtered$KNOWNDECEASEDMO,
  raw_outcomes_df_filtered$KNOWNDECEASEDYR,
  min_year = 1900
)

raw_outcomes_df_filtered$lastalive_date <- clean_month_year(
  raw_outcomes_df_filtered$LASTALIVEMO,
  raw_outcomes_df_filtered$LASTALIVEYR,
  min_year = 1900
)

# Create the hhid_pn for merging
raw_outcomes_df_filtered$hhid_pn <- paste0(as.character(raw_outcomes_df_filtered$HHID), 
                                           "0", 
                                           as.character(raw_outcomes_df_filtered$PN))


# Read in the VBS Data for more detailed collection dates
collection_df <- read_dta("../HRS_Data/Initial/HRS2016VBS.dta")
collection_df <- collection_df %>%
  select("HHID", "PN", "PVBS_N_DAYS")

# Create the hhid_pn
collection_df$hhid_pn <- paste0(as.character(collection_df$HHID), 
                                           "0", 
                                           as.character(collection_df$PN))

# Merge
raw_outcomes_w_collection_df <- merge(raw_outcomes_df_filtered, collection_df,
                                      by="hhid_pn", all=FALSE)

# Filter important variables
key_columns <- c("hhid_pn", "HHID.x", "PN.x", 'death_date', 'lastalive_date',
                     'PIWMONTH', 'PIWYEAR', 'PAGE',
                 'PVBS_N_DAYS')
outcomes_collection_filtered <- raw_outcomes_w_collection_df %>%
  select(all_of(key_columns))

# Create collection date
# First, make it into a date with clean_month_year()
# Then, add the number of days from the VBS data
outcomes_collection_filtered$collection_date <- clean_month_year(
  outcomes_collection_filtered$PIWMONTH,
  outcomes_collection_filtered$PIWYEAR,
  min_year = 1900
)
outcomes_collection_filtered$collection_date <- outcomes_collection_filtered$collection_date + outcomes_collection_filtered$PVBS_N_DAYS

# From here:
# Create the death variable (if there is a death date)
# Create time to end of follow-up which is either death date - collection date or lastalive date - collection date
outcomes_collection_filtered <- outcomes_collection_filtered |>
  mutate(death_status = ifelse(!is.na(death_date), 1, 0))

# Now, if death status == 1
  # death_date == death date - collection date
# If death status == 0
  # death_date == lastalive date - collection date
outcomes_collection_filtered <- outcomes_collection_filtered |>
  mutate(death_time = case_when(
    death_status == 1 ~ death_date - collection_date,
    death_status == 0 ~ lastalive_date - collection_date
  ))
outcomes_collection_filtered$death_time <- as.numeric(outcomes_collection_filtered$death_time)

outcomes_collection_filtered <- outcomes_collection_filtered %>%
  mutate(death_time = pmax(death_time, 0))

## Now it's time to make the change
# It's really pretty simple -- if death_time < 0, set it to 0
hist(outcomes_collection_filtered$death_time)

outcomes_df <- outcomes_collection_filtered %>%
  select(c("HHID.x", "PN.x", "PAGE", "death_status", "death_time"))

# Create the hhid_pn correctly
# First, we need to remove the 0 when the hhid starts with 0
outcomes_df <- outcomes_df |>
  mutate(HHID.x = sub("^0", "", HHID.x))

pheno_df <- pheno_df |>
  mutate(hhid = sub("^0", "", hhid))

outcomes_df$id <- paste0(as.character(outcomes_df$HHID.x), 
                                           as.character(outcomes_df$PN))

pheno_df$id <- paste0(as.character(pheno_df$hhid),
                         as.character(pheno_df$pn))

# Now merge by id
pheno_and_outcomes_df <- merge(pheno_df, outcomes_df, by="id", all=FALSE)

# Select key columns:
pheno_outcome_columns <- c("id", "hhid", "pn", "IDATid", "Ethnicity", "Race", 
                           "Female", "PAGE", "death_status", "death_time")
pheno_and_outcomes_df <- pheno_and_outcomes_df %>%
  select(all_of(pheno_outcome_columns))

# Merge with the clocks
class(pheno_and_outcomes_df$IDATid)
outcomes_and_clocks <- merge(x=pheno_and_outcomes_df, y=Dunedin_df, by="IDATid", all=TRUE)
outcomes_and_clocks <- merge(outcomes_and_clocks, PC_Clock_df, by="IDATid")

write.csv(outcomes_and_clocks, "../HRS_Data/HRS_outcomes_and_clocks.csv")

