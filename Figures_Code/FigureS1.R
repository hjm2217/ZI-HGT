##########################################################################################
## Plot the results for number of different cell types vs. RMSE in the OSCC simulations ##
##########################################################################################

# We do not have the simulation results saved to GitHub as they're too large. 
# The code used to prep the raw simulation results for plotting is included below.
# For the code used to solely plot the results, scroll down past the prep code.
# We have included the prepared (for plotting) simulation results, which may be found in
# Simulations/Simulation_Results/OSCC_Simulations_num_cell_types_vs_RMSE_Reduction.txt
# All actual simulation results may be found on our OSF.

# If you do not wish to go through the preparation process

# Prepping the raw simulation results for plotting #

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
n_cell_types <- c(12,10,8,6)
alpha0 <- c(-0.1, 0.0, 0.1)  #alpha_0 is actually xover.75minusx(invlogit(..)) of these values
alpha1 <- c(0.1, 0.3, 0.5) # to match seq depth sims

grid <- expand.grid(alpha0, alpha1, SC_replicate, ST_replicate, n_cell_types, sample)
colnames(grid) <- c("alpha0", "alpha1", "SC_replicate", "ST_replicate","n_cell_types", "sample")
grid$iseed <- rep(rep(iseed, rep(9, length(iseed))))
grid$ntotal <- ntotal
grid$lib_factor <- 0.05
grid$Phi_factor <- 50


files <- dir("SPARSim_Results/Sample_2/")
res_df <- data.frame()
for (j in 1:3600){
  files <- dir(paste0("SPARSim_Results/Sample_", grid$sample[j], "/"))
  file_name <- paste0("./SPARSim_Results/Sample_", grid$sample[j],
                      "/HGT_SPARSim_pseudo_OSCC_n", grid$ntotal[j],
                      "_cells_n", grid$n_cell_types[j], "_cell_types",
                      "_library_factor_", grid$lib_factor[j],
                      "_Phi_factor_", grid$Phi_factor[j],
                      "_scRNAseq_rep_", grid$SC_replicate[j],
                      "_ST_replicate_", grid$ST_replicate[j],
                      "_alpha0_", grid$alpha0[j],
                      "_alpha1_", grid$alpha1[j], ".rds")
  res <- readRDS(file_name)
  tmp_df <- data.frame(sample = grid$sample[j], n_cell_types = grid$n_cell_types[j], 
                       cell_count = grid$ntotal[j],
                       lib_factor = grid$lib_factor[j], Phi_factor = grid$Phi_factor[j],
                       SC_rep = grid$SC_replicate[j], ST_rep = grid$ST_replicate[j], 
                       alpha0 = grid$alpha0[j], alpha1 = grid$alpha1[j], 
                       CARD_rmse = res$rmse$CARD_rmse, Ave_rmse = res$rmse$Ave_rmse, 
                       j = j, WAIC = res$WAIC,
                       PercentMeanDiff = 100*(1-(res$rmse$Ave_rmse / res$rmse$CARD_rmse)))
  res_df <- rbind(res_df, tmp_df)
  if (j %% 10 == 0){
    gc()
    cat(j)
  }
}

# Add all 14 cell types too
for (j in 1:900){
  files <- dir(paste0("SPARSim_Results/Sample_", grid$sample[j], "/"))
  file_name <- paste0("HGT_SPARSim_pseudo_OSCC_n", grid$ntotal[j], 
                      "_cells_library_factor_", grid$lib_factor[j],
                      "_Phi_factor_", grid$Phi_factor[j],
                      "_scRNAseq_rep_", grid$SC_replicate[j],
                      "_ST_replicate_", grid$ST_replicate[j],
                      "_alpha0_", grid$alpha0[j], 
                      "_alpha1_", grid$alpha1[j], ".rds")
  if (!(file_name %in% files)){next}
  res <- readRDS(paste0("SPARSim_Results/Sample_", grid$sample[j], "/", file_name))
  tmp_df <- data.frame(sample = grid$sample[j], n_cell_types = 14, 
                       cell_count = grid$ntotal[j],
                       lib_factor = grid$lib_factor[j], Phi_factor = grid$Phi_factor[j],
                       SC_rep = grid$SC_replicate[j], ST_rep = grid$ST_replicate[j], 
                       alpha0 = grid$alpha0[j], alpha1 = grid$alpha1[j], 
                       CARD_rmse = res$rmse$CARD_rmse, Ave_rmse = res$rmse$Ave_rmse, 
                       j = j, WAIC = res$WAIC,
                       PercentMeanDiff = 100*(1-(res$rmse$Ave_rmse / res$rmse$CARD_rmse)))
  res_df <- rbind(res_df, tmp_df)
  if (j %% 10 == 0){
    gc()
    cat(j)
  }
}

write.table(res_df, "../Simulations/Simulation_Results/OSCC_Simulations_num_cell_types_vs_RMSE_Reduction.txt", sep = "\t", row.names=F)


#############
# Make plot #
#############

res_df <- read.table("../Simulations/Simulation_Results/OSCC_Simulations_num_cell_types_vs_RMSE_Reduction.txt", header = T)

# Choose hyperparameters based on WAIC
res_df$best_WAIC <- NA
for(j in 1:nrow(res_df)){
  
  #Check if this row minimizes the WAIC
  res_df_WAIC <- res_df[res_df$n_cell_types == res_df$n_cell_types[j] &
                        res_df$SC_rep == res_df$SC_rep[j] & res_df$ST_rep == res_df$ST_rep[j], ]
  min_WAIC <- min(res_df_WAIC$WAIC)
  res_df$best_WAIC[j] <- ifelse(res_df$WAIC[j] == min_WAIC, TRUE, FALSE)
}
res_df <- res_df[res_df$best_WAIC, ]

# Now plot
res_df$nct_f <- factor(res_df$n_cell_types)
g <- ggplot(res_df, aes(x = nct_f, y = PercentMeanDiff)) + geom_boxplot(outlier.shape = NA, fill = "green4") + theme_bw() +
  labs(title = "RMSE Reduction vs. Number of Cell Types", subtitle = "OSCC Sample 2"
       , x = "Number of Cell Types", y = "RMSE Reduction (%)") + coord_cartesian(ylim = c(-15, 25))
g
ggsave("../Figures/Supplementary_Figures/FigureS1.png", g, "png", width = 8, units = "in")



