####################################################################
## Plot the results for Sparsity vs. RMSE in the OSCC simulations ##
####################################################################

# We do not have the simulation results saved to GitHub as they're too large. 
# The code used to prep the raw simulation results for plotting is included below.
# For the code used to solely plot the results, scroll down past the prep code.
# We have included the prepared (for plotting) simulation results, which may be found in
# Simulations/Simulation_Results/OSCC_Simulations_Sparsity_vs_RMSE_Reduction.txt
# All actual simulation results may be found on our OSF.

# Run this code from the parent directory where your results are stored
# We named the directory containing simulation results for sample ##sample_number##
# SPARSim_Results/Sample_##sample_number##/

###############################
# Prep the simulation results #
###############################

iseed <- c(56789, 123456 + 111111*0:98)
sample <- 2
ntotal <- 10
SC_replicate <- 1:10
ST_replicate <- 1:10
alpha0 <- c(-0.1, 0.0, 0.1)  #alpha_0 is actually xover.75minusx(invlogit(..)) of these four values
alpha1 <- c(0.1, 0.3, 0.5)

grid_temp <- expand.grid(alpha0, alpha1, SC_replicate, ST_replicate, sample)
colnames(grid_temp) <- c( "alpha0", "alpha1", "SC_replicate", "ST_replicate", "sample")
grid_temp$iseed <- rep(rep(iseed, rep(9, length(iseed))))
grid_temp$ntotal <- ntotal

lib_factor <- c(0.05, 0.1, 0.2, 0.4, 0.6)
Phi_factor <- c(50, 10, 4, 4, 4)
grid <- data.frame(alpha0 = rep(grid_temp$alpha0, 5), alpha1 = rep(grid_temp$alpha1, 5), SC_replicate = rep(grid_temp$SC_replicate, 5), ST_replicate =  rep(grid_temp$ST_replicate, 5),
                   lib_factor = rep(lib_factor, rep(900,5)), Phi_factor = rep(Phi_factor, rep(900, 5)), sample = 2, ntotal = 10)
files <- dir("SPARSim_Results/Sample_2/")
res_df <- data.frame()
for (j in 1:4500){
  file_name <- paste0("HGT_SPARSim_pseudo_OSCC_n", grid$ntotal[j], 
                      "_cells_library_factor_", grid$lib_factor[j],
                      "_Phi_factor_", grid$Phi_factor[j],
                      "_scRNAseq_rep_", grid$SC_replicate[j],
                      "_ST_replicate_", grid$ST_replicate[j],
                      "_alpha0_", grid$alpha0[j], 
                      "_alpha1_", grid$alpha1[j], ".rds")
  if (!(file_name %in% files)){next}
  res <- readRDS(paste0("SPARSim_Results/Sample_2/", file_name))
  tmp_df <- data.frame(sample = grid$sample[j], cell_count = grid$ntotal[j],
                       lib_factor = grid$lib_factor[j], Phi_factor = grid$Phi_factor[j],
                       SC_rep = grid$SC_replicate[j], ST_rep = grid$ST_replicate[j], 
                       alpha0 = grid$alpha0[j], alpha1 = grid$alpha1[j], 
                       CARD_rmse = res$rmse$CARD_rmse, Ave_rmse = res$rmse$Ave_rmse, 
                       j = j, WAIC = res$WAIC,
                       PercentMeanDiff = (1 - (res_df$Ave_rmse / res_df$CARD_rmse))*100)
  res_df <- rbind(res_df, tmp_df)
  if (j %% 10 == 0){
    gc()
    cat(j)
  }
}

# Determine sparsity levels for ease of plotting notation
res_df$sparsity_level <- NA
res_df$best_WAIC <- NA
for(j in 1:nrow(res_df)){
  if (res_df$lib_factor[j] == 0.05 & res_df$Phi_factor[j] ==  50){
    res_df$sparsity_level[j] <- 1
  }
  if (res_df$lib_factor[j] == 0.1 & res_df$Phi_factor[j] ==  10){
    res_df$sparsity_level[j] <- 2
  }
  if (res_df$lib_factor[j] == 0.2 & res_df$Phi_factor[j] ==  4){
    res_df$sparsity_level[j] <- 3
  }
  if (res_df$lib_factor[j] == 0.4 & res_df$Phi_factor[j] ==  4){
    res_df$sparsity_level[j] <- 4
  }
  if (res_df$lib_factor[j] == 0.6 & res_df$Phi_factor[j] ==  4){
    res_df$sparsity_level[j] <- 5
  }
  
  # Find which set of hyperparameters minimizes the WAIC
  res_df_WAIC <- res_df[res_df$lib_factor == res_df$lib_factor[j] & res_df$Phi_factor == res_df$Phi_factor[j] &
                        res_df$SC_rep == res_df$SC_rep[j] & res_df$ST_rep == res_df$ST_rep[j], ]
  min_WAIC <- min(res_df_WAIC$WAIC)
  res_df$best_WAIC[j] <- ifelse(res_df$WAIC[j] == min_WAIC, TRUE, FALSE)
}

res_df <- res_df[res_df$best_WAIC, ]
res_df$sparsity_level_f <- factor(res_df$sparsity_level)

# Save the dataframe now that it is fully prepped
write.table(res_df, "../Simulations/Simulation_Results/OSCC_Simulations_Sparsity_vs_RMSE_Reduction.txt", sep = "\t", row.names=F)

#################
# Make the plot #
#################

# Load in the data and plot #
res_df <- read.table("../Simulations/Simulation_Results/OSCC_Simulations_Sparsity_vs_RMSE_Reduction.txt", header = T)

library(ggplot2)
g <- ggplot(res_df_use, aes(x = sparsity_level_f, y = PercentMeanDiff)) + geom_boxplot(outlier.shape = NA, fill = "green4") + theme_bw() +
  labs(title = "RMSE Reduction vs. Sparsity", subtitle = "OSCC Sample 2", x = "Sparsity Level", y = "RMSE Reduction (%)") + coord_cartesian(ylim = c(-12, 20))
g
ggsave("../Figures/Figure3.png", g, "png", width = 8, units = "in")
