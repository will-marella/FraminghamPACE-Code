# Load packages
library(readr)
library(dplyr)
library(haven)

########################################################
# Load and merge data
########################################################

# FraminghamPACE
FraminghamPACE <- read_csv("../CARDIA_Data/CARDIA_combined_FraminghamPACE.csv")
FraminghamPACE <- FraminghamPACE %>%
  rename(sample_name = `...1`)
FraminghamPACE$sample_name <- as.character(FraminghamPACE$sample_name)

# Load DunedinPACE
DunedinPACE <- read_csv("../CARDIA_Data/CARDIA_combined_DunedinPACE.csv")
DunedinPACE <- DunedinPACE %>%
  dplyr::select(filenames, DunedinPACE) %>%
  rename(sample_name = filenames)
DunedinPACE$sample_name <- as.character(DunedinPACE$sample_name)

# Load PC Clocks
PC_Clocks <- read_csv("../CARDIA_Data/CARDIA_combined_PC_clocks.csv")
PC_Clocks <- PC_Clocks %>%
  dplyr::select(filenames, 
         PCPhenoAge, PCGrimAge,
         PCPhenoAgeResid, PCGrimAgeResid,
         PCHorvath1, PCHannum,
         PCHorvath1Resid, PCHannumResid) %>%
  rename(sample_name = filenames)
PC_Clocks$sample_name <- as.character(PC_Clocks$sample_name)

# Prepare phenotype info
pheno <- read_csv("../CARDIA_Data/Raw/CARDIA_combined_clean_pheno.csv")
pheno <- pheno %>%
  dplyr::select(filenames, barcode, dbgap_subject_id, collection_visit, collection_year,
         age_at_collection, sex) %>%
  rename(sample_name = filenames)
pheno$sample_name <- as.character(pheno$sample_name)

# Add race
race <- read_csv("../CARDIA_Data/Raw/CARDIA_race_data.csv")

race <- race %>%
  rename(dbgap_subject_id = dbGaP_Subject_ID)

# Then do a left join to add race data to pheno
pheno <- pheno %>%
  left_join(race %>% dplyr::select(dbgap_subject_id, Race, Race_char), 
            by = "dbgap_subject_id")

# Merge all
PACE_clocks <- merge(FraminghamPACE, DunedinPACE, by="sample_name")
all_clocks <- merge(PACE_clocks, PC_Clocks, by="sample_name")

combined_data_raw <- merge(pheno, all_clocks, by="sample_name")

## Add cell counts
cells <- read_csv("../CARDIA_Data/CARDIA_combined_Salas_cells.csv")
cells <- cells %>%
  rename(sample_name=filenames)

combined_data <- merge(combined_data_raw, cells, by="sample_name")

################################
# Define clock variables
################################

# Define all clock variables that need residuals
framingham_vars <- c("ElasticNet", "Ridge")
special_vars <- c("DunedinPACE")  # Special FraminghamPoA variable
already_resid_vars <- c("PCGrimAgeResid", "PCPhenoAgeResid")  # Already have Resid versions
non_resid_vars <- c("PCGrimAge", "PCPhenoAge")  # Need scaling but no residuals

# All variables that need residuals
vars_needing_residuals <- c(framingham_vars, special_vars)

########################################################
# Add resid columns
########################################################

# Create age residuals for variables that need them
resid_names <- paste0(vars_needing_residuals, "Resid")

combined_data <- combined_data %>%
  # Create all residuals at once
  mutate(across(all_of(vars_needing_residuals), 
                ~ resid(lm(.x ~ age_at_collection, data = combined_data)), 
                .names = "{.col}Resid")) %>%
  # Scale all variables at once
  mutate(
    # Scale original variables that don't get residuals
    across(all_of(c(special_vars, non_resid_vars, already_resid_vars)), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled"),
    # Scale all residual variables  
    across(all_of(resid_names), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled"),
    # Scale the framingham variables (if they need original scaled versions)
    across(all_of(framingham_vars), 
           ~ as.vector(scale(.x)), 
           .names = "{.col}_scaled")
  )

########################################################
# Create a proper time variable
# Remove technical replicates
# And make female column
########################################################

full_data <- combined_data %>%
  group_by(dbgap_subject_id) %>%
  mutate(time_since_baseline = collection_year - min(collection_year)) %>%
  ungroup()

# Automatically identify clock columns by prefix
clock_prefixes <- c("^Framingham", "^Dunedin", "^PC")
clock_cols <- grep(paste(clock_prefixes, collapse = "|"), 
                   names(full_data), value = TRUE)

# Hardcoded cell columns
cell_cols <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", 
               "Eos", "Mono", "Neu", "NK", "Treg")

# Combine all columns that should be averaged
avg_cols <- c(clock_cols, cell_cols)

# Unique identifier columns (take first)
unique_cols <- c("sample_name", "barcode")

clean_data <- full_data %>%
  group_by(dbgap_subject_id, collection_visit) %>%
  summarise(
    # Take first instance of unique identifier columns
    across(all_of(unique_cols), first),
    
    # Take mean of clock and cell columns across technical replicates (ignoring NA values)
    across(all_of(avg_cols), ~ mean(.x, na.rm = TRUE)),
    
    # For any remaining columns not specified above, take the first value
    # (assuming they should be the same across technical replicates)
    across(everything(), first),
    
    # Count number of technical replicates for this subject-visit combination
    n_tech_replicates = n(),
    
    .groups = "drop"
  )

clean_data$Female <- clean_data$sex - 1

########################################################
# Export
########################################################

## Write combined data
write.csv(clean_data, "../CARDIA_Data/CARDIA_combined_data.csv")

