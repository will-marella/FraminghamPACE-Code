###############################################
# Load libraries
###############################################

library(readr)
library(dplyr)

###############################################
# Load, merge, and clean datasets
###############################################

Framingham_df <- read_csv("../Lehne_Data/Lehne_FraminghamPACE_df.csv")
pheno_df <- read_csv("../Lehne_Data/Clean_Pheno_FullDataset.csv")

combined_df <- merge(pheno_df, Framingham_df, by.x="SampleID", by.y="...1")

clean_df <- combined_df %>%
  dplyr::select(-c("...1", "BARCODE")) %>%
  dplyr::rename(BARCODE = SampleID)

###############################################
# Export clean dataset
###############################################

write.csv(clean_df, "../Lehne_Data/Lehne_Combined_Data.csv")

