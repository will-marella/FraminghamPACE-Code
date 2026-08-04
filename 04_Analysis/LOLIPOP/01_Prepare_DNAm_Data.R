# Probe filter Lehne

# Load probe list
probe_list <- read.csv("../../FittingPoA/Data/Initial/probelist_450k_epicv1_epicv2_sugden.csv")
probe_list <- probe_list$probe

### Load in chunks

load("../Lehne_Data/Lehne_Bvals_Replicates_Only.RData")
betas <- as.data.frame(dat)


betas_t <- t(betas)
betas_t <- as.data.frame(betas_t)
betas_probefiltered <- betas_t[, colnames(betas_t) %in% probe_list]


# Save as .csv
write.csv(betas_probefiltered, "../Lehne_Data/Lehne_SugdenProbes_replicates.csv")