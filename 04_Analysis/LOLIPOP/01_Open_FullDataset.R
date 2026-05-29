library(data.table)
library(dplyr)

Lehne_data <- fread("../Lehne_Data/GSE55763_normalized_betas.txt", 
                    header = TRUE)

# Remove all columns named exactly "Detection Pval"
Lehne_data <- Lehne_data[, !names(Lehne_data) %in% "Detection Pval", with = FALSE]

head(rownames(Lehne_data))
tail(rownames(Lehne_data))
head(colnames(Lehne_data))
tail(colnames(Lehne_data))

# Convert to df
Lehne_data <- as.data.frame(Lehne_data)

# Convert ID_REF to rownames, remove column
rownames(Lehne_data) <- Lehne_data$ID_REF
Lehne_data <- Lehne_data %>%
  dplyr::select(-ID_REF)

# Save as RDS
saveRDS(Lehne_data, "../Lehne_Data/Lehne_Bvals_Whole_Sample.rds")

