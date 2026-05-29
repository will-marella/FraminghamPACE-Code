# create_CALERIE_sugden_probes.R

library(readr)
library(glmnet)

data <- load("../CALERIE_Data/Initial/CALERIE_EPIC_data.Rdata")
probe_list <- read_csv("../../FittingPoA/Data/Initial/probelist_450k_epicv1_epicv2.csv")
probe_list <- as.list(probe_list)


betas <- as.data.frame(betas)
beta_matrix <- betas[rownames(betas) %in% probe_list$probe,]

CALERIE_probefiltered <- t(beta_matrix)
CALERIE_probefiltered <- as.data.frame(CALERIE_probefiltered)

write.csv(CALERIE_probefiltered, "../CALERIE_Data/CALERIE_probefiltered.csv")
