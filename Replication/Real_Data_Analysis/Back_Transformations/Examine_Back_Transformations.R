################################################################
## Studying impact of transforming and back-transforming data ##
################################################################

# We ran this analysis on our university's high performance computer, which uses the Slurm
# workload manager.  For the exact specifications of the job we ran, see the 
# Real_Data_Analysis.sh file.  You may need to adjust some of the options depending on the
# computational resources you have access to.

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

library(dplyr)
library(Seurat)
library(dplyr)
library(LaplacesDemon)

############################
# Transformation Functions #
############################

transformX2H <- function(X, alpha0, alpha1, seed = 123456){
  set.seed(seed)
  n.zeros <- sum(colSums(X == 0))
  n.nonzeros <- sum(colSums(X > 0))
  H <- X  #leave any nonzero entries the same
  H[H > 0] <- (rgamma(n.nonzeros, shape = alpha1 + H[H > 0], rate = alpha1 + 1))
  H[H == 0] <- logit(rbeta(n.zeros, (alpha0*0.75), (0.25*alpha0+1))) #letting the prior be Beta(alpha-kappa, kappa), where kappa=0.25*alpha
  return(H)
}

backTransformH2X <- function(H, seed = 123456){
  set.seed(seed)
  X0 <- H
  X0[X0 <= 0] <- 0  
  binom_prob <- as.vector(exp(H[X0 > 0])/(1 + exp(H[X0 > 0])))
  binom_prob[!is.finite(binom_prob)] <- 1 #sometimes exponential of H is too large, but in practice this is (big / big +1) = 1
  
  X0[X0 > 0] <- rbinom(length(H[X0 > 0]), size = 1, prob = binom_prob)
  X <- X0 # leave all zero entries as zero
  X[X0 > 0] <- rpois(length(X[X0 > 0]), lambda = H[X0 > 0])
  return(X)
}

####################
# Load in the data #
####################

sample <- readRDS(paste0("Data/sample_", job.id, ".rds"))
spatial_count <- sample@assays$SCT@counts

############
# Settings #
############
alpha0 <- 2.0  # equivalent to alpha_0 = 1.5, kappa_0 = 0.5 in paper (beta(1.5,1.5))
alpha1 <- 0.1
seed_transform <- 123456 + 0:9
seed_back_transform <- 654321 + 0:9

################################
# Transform and back-transform #
################################
res <- data.frame(seed_transform = seed_transform,
                  seed_back_transform = seed_back_transform,
                   mean_diff = NA)

for (i in 1:nrow(res)){
  spatial_count_transform <- transformX2H(
    X = spatial_count,
    alpha0 = alpha0,
    alpha1 = alpha1,
    seed = res$seed_transform[i])
  
  spatial_count_back_transform <- backTransformH2X(
    H = spatial_count_transform,
    seed = res$seed_back_transform[i])
  
  res$mean_diff[i] <- mean(spatial_count - spatial_count_back_transform)
  
  gc()
  cat("Completed", i, "of", nrow(res), "transformations.\n")
}

write.csv(res, file = paste0("./mean_diff_transformation_and_back_sample", job.id, ".csv"), 
          row.names = FALSE, quote = FALSE)



