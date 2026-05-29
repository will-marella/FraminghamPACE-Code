# Load probe list
probe_list <- read.csv("../../FittingPoA/Data/Initial/probelist_450k_epicv1_epicv2.csv")
probe_list <- probe_list$probe

### Load in chunks

# Chunk 1
beta_chunk1 <- readRDS("../HRS_Data/Initial/GRSet_fully_filtered_bmiq_chunk1.rds")
beta_chunk1 <- t(beta_chunk1)
beta_chunk1 <- as.data.frame(beta_chunk1)
beta_chunk1 <- beta_chunk1[, colnames(beta_chunk1) %in% probe_list]

# Chunk 2
beta_chunk2 <- readRDS("../HRS_Data/Initial/GRSet_fully_filtered_bmiq_chunk2.rds")
beta_chunk2 <- t(beta_chunk2)
beta_chunk2 <- as.data.frame(beta_chunk2)
beta_chunk2 <- beta_chunk2[, colnames(beta_chunk2) %in% probe_list]

# Chunk 3
beta_chunk3 <- readRDS("../HRS_Data/Initial/GRSet_fully_filtered_bmiq_chunk3.rds")
beta_chunk3 <- t(beta_chunk3)
beta_chunk3 <- as.data.frame(beta_chunk3)
beta_chunk3 <- beta_chunk3[, colnames(beta_chunk3) %in% probe_list]

# Chunk 4
beta_chunk4 <- readRDS("../HRS_Data/Initial/GRSet_fully_filtered_bmiq_chunk4.rds")
beta_chunk4 <- t(beta_chunk4)
beta_chunk4 <- as.data.frame(beta_chunk4)
beta_chunk4 <- beta_chunk4[, colnames(beta_chunk4) %in% probe_list]

# Chunk 5
beta_chunk5 <- readRDS("../HRS_Data/Initial/GRSet_fully_filtered_bmiq_chunk5.rds")
beta_chunk5 <- t(beta_chunk5)
beta_chunk5 <- as.data.frame(beta_chunk5)
beta_chunk5 <- beta_chunk5[, colnames(beta_chunk5) %in% probe_list]

# Find shared features across all chunks
shared_features <- Reduce(intersect, list(colnames(beta_chunk1), colnames(beta_chunk2), 
                                          colnames(beta_chunk3), colnames(beta_chunk4), 
                                          colnames(beta_chunk5)))

# Subset each chunk to only include shared features
beta_chunk1 <- beta_chunk1[, shared_features]
beta_chunk2 <- beta_chunk2[, shared_features]
beta_chunk3 <- beta_chunk3[, shared_features]
beta_chunk4 <- beta_chunk4[, shared_features]
beta_chunk5 <- beta_chunk5[, shared_features]

### Combine chunks
all_betas <- rbind(beta_chunk1, beta_chunk2, beta_chunk3, beta_chunk4, beta_chunk5)

# Save as .csv
write.csv(all_betas, "../HRS_Data/HRS_betas_450k_epicv1_epicv2.csv")


