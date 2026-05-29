# Load probe list
probe_list <- read.csv("../../FittingPoA/Data/Initial/probelist_450k_epicv1_epicv2_sugden.csv")
probe_list <- probe_list$probe

### Load in chunks

# Chunk 1
beta_chunk1 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk1.rds")
beta_chunk1 <- t(beta_chunk1)
beta_chunk1 <- as.data.frame(beta_chunk1)
beta_chunk1 <- beta_chunk1[, colnames(beta_chunk1) %in% probe_list]

gc()

# Chunk 2
beta_chunk2 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk2.rds")
beta_chunk2 <- t(beta_chunk2)
beta_chunk2 <- as.data.frame(beta_chunk2)
beta_chunk2 <- beta_chunk2[, colnames(beta_chunk2) %in% probe_list]

gc()

# Chunk 3
beta_chunk3 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk3.rds")
beta_chunk3 <- t(beta_chunk3)
beta_chunk3 <- as.data.frame(beta_chunk3)
beta_chunk3 <- beta_chunk3[, colnames(beta_chunk3) %in% probe_list]

gc()

# Chunk 4
beta_chunk4 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk4.rds")
beta_chunk4 <- t(beta_chunk4)
beta_chunk4 <- as.data.frame(beta_chunk4)
beta_chunk4 <- beta_chunk4[, colnames(beta_chunk4) %in% probe_list]

gc()

# Chunk 5
beta_chunk5 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk5.rds")
beta_chunk5 <- t(beta_chunk5)
beta_chunk5 <- as.data.frame(beta_chunk5)
beta_chunk5 <- beta_chunk5[, colnames(beta_chunk5) %in% probe_list]

gc()

# Chunk 6
beta_chunk6 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk6.rds")
beta_chunk6 <- t(beta_chunk6)
beta_chunk6 <- as.data.frame(beta_chunk6)
beta_chunk6 <- beta_chunk6[, colnames(beta_chunk6) %in% probe_list]

gc()

# Chunk 7
beta_chunk7 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk7.rds")
beta_chunk7 <- t(beta_chunk7)
beta_chunk7 <- as.data.frame(beta_chunk7)
beta_chunk7 <- beta_chunk7[, colnames(beta_chunk7) %in% probe_list]

gc()

# Chunk 8
beta_chunk8 <- readRDS("../CARDIA_Data/Raw/GRSet_fully_filtered_bmiq_chunk8.rds")
beta_chunk8 <- t(beta_chunk8)
beta_chunk8 <- as.data.frame(beta_chunk8)
beta_chunk8 <- beta_chunk8[, colnames(beta_chunk8) %in% probe_list]

gc()


### Combine chunks
all_betas <- rbind(beta_chunk1, beta_chunk2, beta_chunk3, beta_chunk4, 
                   beta_chunk5, beta_chunk6, beta_chunk7, beta_chunk8)

gc()

# Save as .csv
write.csv(all_betas, "../CARDIA_Data/CARDIA_SugdenProbes.csv")