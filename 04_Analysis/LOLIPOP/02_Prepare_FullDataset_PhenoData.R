library(readr)
library(dplyr)

raw_data <- read_csv("../Lehne_Data/GSE55763_BioAge_MK05042021.csv")

names(raw_data)

clean_data <- raw_data %>%
  dplyr::select(SampleID, age.ch1, gender.ch1) %>%
  dplyr::rename(Age = age.ch1,
                Gender = gender.ch1)

write.csv(clean_data, "../Lehne_Data/Clean_Pheno_FullDataset.csv")



