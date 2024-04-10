###################################################################################
## Plot the results for the other methods compared to the ZI-HGT + CARD and CARD ##
###################################################################################

# We do not have the simulation results saved to GitHub as they're too large. 
# The code used to prep the raw simulation results for plotting is included below.
# For the code used to solely plot the results, scroll down past the prep code.
# We have included the prepared (for plotting) simulation results, which may be found in
# Simulations/Simulation_Results/OSCC_Simulations_Method_Comp_incl_DT.txt
# All actual simulation results may be found on our OSF.

###############################
# Prep the simulation results #
###############################

# Run this code from the parent directory where your results are stored
# We named the directory containing simulation results for sample ##sample_number##
# SPARSim_Results/Sample_##sample_number##/

iseed <- c(56789, 123456 + 111111*0:98)
sample <- 2
ntotal <- 10
SC_replicate <- 1:10
ST_replicate <- 1:10
alpha0 <- c(-0.1, 0.0, 0.1)  #alpha_0 is actually xover.75minusx(invlogit(..)) of these values
alpha1 <- c(0.1, 0.3, 0.5) 

grid <- expand.grid(alpha0, alpha1, SC_replicate, ST_replicate, sample)
colnames(grid) <- c("alpha0", "alpha1", "SC_replicate", "ST_replicate", "sample")
grid$iseed <- rep(rep(iseed, rep(9, length(iseed))))
grid$ntotal <- ntotal
grid$lib_factor <- 0.05
grid$Phi_factor <- 50


# Get the ZI-HGT + CARD results
res_df <- data.frame()
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
                       Median_rmse = res$rmse$Median_rmse, j = j, WAIC = res$WAIC,
                       PercentMedianDiff = 100*(1-(res$rmse$Median_rmse / res$rmse$CARD_rmse)),
                       PercentMeanDiff = 100*(1-(res$rmse$Ave_rmse / res$rmse$CARD_rmse)))
  res_df <- rbind(res_df, tmp_df)
  if (j %% 10 == 0){
    gc()
    cat(j)
  }
}

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
rownames(res_df) <- 1:100

# Now get the results for the other methods
grid <- expand.grid(alpha0, alpha1, SC_replicate, ST_replicate, sample)
colnames(grid) <- c("alpha0", "alpha1", "SC_replicate", "ST_replicate", "sample")
grid$ntotal <- ntotal
grid$lib_factor <- 0.05
grid$Phi_factor <- 50

source("../Utilities/HGTfunctions3.R")

res_all_df <- data.frame()
for (j in 1:100){
  
  # Will need CARD results for comparison (alpha0 and alpha1 don't matter for CARD alone, which is why we aren't gridding on those here)
  file_name <- paste0("HGT_SPARSim_pseudo_OSCC_n", grid$ntotal[j], 
                      "_cells_library_factor_", grid$lib_factor[j],
                      "_Phi_factor_", grid$Phi_factor[j],
                      "_scRNAseq_rep_", grid$SC_replicate[j],
                      "_ST_replicate_", grid$ST_replicate[j],
                      "_alpha0_0", 
                      "_alpha1_0.3", ".rds")
  res <- readRDS(paste0("SPARSim_Results/Sample_", grid$sample[j], "/", file_name))
  true_prop <- res$cell_type_proportion_matrices$true_ct_prop_mat
  CARD_rmse <- res$rmse$CARD_rmse
  
  # SpatialDecon
  SpatialDecon_prop <- read.csv(paste0("SPARSim_Results/Sample_2/SpatialDecon/SpatialDecon_SPARSim_pseudo_OSCC_n",
                                         grid$ntotal[j], "_cells_library_factor_", grid$lib_factor[j],
                                         "_Phi_factor_", grid$Phi_factor[j], "_scRNAseq_rep_", grid$SC_replicate[j],
                                         "_ST_replicate_", grid$ST_replicate[j], ".rds"), header=T)
  rownames(SpatialDecon_prop) <- SpatialDecon_prop[ ,1]
  if (dim(SpatialDecon_prop)[2] != 15){next} # a couple of these proportions only have 12 cell types, we'll skip them
  SpatialDecon_prop <- SpatialDecon_prop[ ,2:15]
  SpatialDecon_rmse <- rmse(SpatialDecon_prop, true_prop)
  
  tmp_df1 <- data.frame(Method = "SpatialDecon", 
                        SC_rep = grid$SC_replicate[j], ST_rep = grid$ST_replicate[j], 
                        RMSE =  SpatialDecon_rmse,
                        CARD_RMSE = CARD_rmse,
                        PercentDiff = 100*(1-(SpatialDecon_rmse / CARD_rmse)))
  
  
  # SPOTlight
  SPOTlight_prop <- read.csv(paste0("SPARSim_Results/Sample_2/SPOTlight/SPOTlight_SPARSim_pseudo_OSCC_n",
                                       grid$ntotal[j], "_cells_library_factor_", grid$lib_factor[j],
                                       "_Phi_factor_", grid$Phi_factor[j], "_scRNAseq_rep_", grid$SC_replicate[j],
                                       "_ST_replicate_", grid$ST_replicate[j], ".rds"), header=T)
  rownames(SPOTlight_prop) <- SPOTlight_prop[ ,1]
  SPOTlight_prop <- SPOTlight_prop[ ,2:15]
  SPOTlight_prop <- SPOTlight_prop[, c("myofibroblast", "cancer.cell", "B.cell", "ecm.myCAF",
                                       "Intermediate.fibroblast", "detox.iCAF", "macrophage",
                                       "endothelial", "dendritic.", "mast",
                                       "conventional.CD4..T.helper.cells", "cytotoxic.CD8..T.",
                                       "Tregs", "cytotoxic.CD8..T.exhausted")]
  SPOTlight_rmse <- rmse(SPOTlight_prop, true_prop)
  tmp_df2 <- data.frame(Method = "SPOTlight", 
                        SC_rep = grid$SC_replicate[j], ST_rep = grid$ST_replicate[j], 
                        RMSE =  SPOTlight_rmse,
                        CARD_RMSE = CARD_rmse,
                        PercentDiff = 100*(1-(SPOTlight_rmse / CARD_rmse)))
  
  # STdeconvolve
  STdeconvolve_prop <- read.csv(paste0("SPARSim_Results/Sample_2/STdeconvolve/STdeconvolve_SPARSim_pseudo_OSCC_n",
                                    grid$ntotal[j], "_cells_library_factor_", grid$lib_factor[j],
                                    "_Phi_factor_", grid$Phi_factor[j], "_scRNAseq_rep_", grid$SC_replicate[j],
                                    "_ST_replicate_", grid$ST_replicate[j], ".rds"), header=T)
  rownames(STdeconvolve_prop) <- STdeconvolve_prop[ ,1]
  STdeconvolve_prop <- STdeconvolve_prop[ ,2:15]
  STdeconvolve_rmse <- rmse(STdeconvolve_prop, true_prop)
  tmp_df3 <- data.frame(Method = "STdeconvolve", 
                        SC_rep = grid$SC_replicate[j], ST_rep = grid$ST_replicate[j], 
                        RMSE =  STdeconvolve_rmse,
                        CARD_RMSE = CARD_rmse,
                        PercentDiff = 100*(1-(STdeconvolve_rmse / CARD_rmse)))
  res_all_df <- rbind(res_all_df, tmp_df1, tmp_df2, tmp_df3)
  gc()
  cat(j)
}

# Add CARD results with the deterministic transformation
iseed <- c(56789, 123456 + 111111*0:98)
sample <- 2
ntotal <- 10
SC_replicate <- 1:10
ST_replicate <- 1:10
grid <- expand.grid(SC_replicate, ST_replicate, sample)
colnames(grid) <- c( "SC_replicate", "ST_replicate", "sample")
grid$iseed <- iseed
grid$ntotal <- ntotal
grid$lib_factor <- 0.05
grid$Phi_factor <- 50
res_df_DT <- data.frame()
for (j in 1:100){
  file_name <- paste0("SPARSim_Results/Sample_2/DeterministicTransformation/DT_SPARSim_pseudo_OSCC_n", grid$ntotal[j],
                      "_cells_library_factor_", grid$lib_factor[j],
                      "_Phi_factor_", grid$Phi_factor[j],
                      "_scRNAseq_rep_", grid$SC_replicate[j],
                      "_ST_replicate_", grid$ST_replicate[j], ".rds")
  res <- readRDS(file_name)
  tmp_df <- data.frame(SC_rep = grid$SC_replicate[j], ST_rep = grid$ST_replicate[j], 
                       RMSE = res$rmse_T, CARD_RMSE = res$rmse_CARD, 
                       PercentDiff = 100*(1-(res$rmse_T / res$rmse_CARD)))
  res_df_DT <- rbind(res_df_DT, tmp_df)
  if (j %% 10 == 0){
    gc()
    cat(j)
  }
}
res_df_DT$Method <- "Det. Transformation + CARD"

# Make the HGT and DT results match the same formatting (and skip the problematic results for SpatialDecon)
res_df_HGT <- res_df[-(c(1:9,91:99)), c("SC_rep", "ST_rep", "Ave_rmse", "CARD_rmse", "PercentMeanDiff")]
colnames(res_df_HGT) <- c("SC_rep", "ST_rep", "RMSE", "CARD_RMSE", "PercentDiff")
res_df_HGT$Method <- "HGT + CARD"

res_df_DT <- res_df_DT[-c(1:9,91:99),]

res_all_df <- rbind(res_all_df, res_df_HGT, res_df_DT)
write.table(res_all_df, "../Simulations/Simulation_Results/OSCC_Simulations_Method_Comp_incl_DT.txt", sep = "\t", row.names=F)


#############
# Make plot #
#############

library(ggplot2)
res_all_df <- read.table("../Simulations/Simulation_Results/OSCC_Simulations_Method_Comp_incl_DT.txt", header = T)
res_all_df$Method[res_all_df$Method == "HGT + CARD"] <- "ZI-HGT + CARD"

res_all_df$Method_f <- factor(res_all_df$Method, levels = c("ZI-HGT + CARD", "Det. Transformation + CARD", "SPOTlight",
                                                            "SpatialDecon", "STdeconvolve"))

# Plot results
g <- ggplot(res_all_df, aes(x = Method_f, y = PercentDiff)) + geom_boxplot(outlier.shape = NA, fill = "green4") + theme_bw() +
  labs(title = "RMSE Reduction vs. CARD", subtitle = "OSCC Sample 2", x = "Method", y = "RMSE Reduction (%)") 
g
ggsave("../Figures/Supplementary_Figures/FigureS3.png", g, "png", width = 8, units = "in")



