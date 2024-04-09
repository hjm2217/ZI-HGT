####################################################################
## Plot the results for Scenario vs. RMSE in the CARD simulations ##
####################################################################

# We do not have the simulation results saved to GitHub as they're too large. 
# The code used to prep the raw simulation results for plotting is included below.
# For the code used to solely plot the results, scroll down past the prep code.
# We have included the prepared (for plotting) simulation results, which may be found in
# ../Simulations/Simulation_Results/HGTpCARD_SeqDepthSims_res.txt and 
# ../Simulations/Simulation_Results/SeqDepthSimsCARDres5.txt

# If you do not wish to go through the preparation process skip to 
#######################################
# Combine the data and make the plots #
#######################################

# Prepping the raw simulation results for plotting #

# Run this code from the parent directory where your results are stored - likely Simulation_Results
# We named the directory containing simulation results for the CARD simulations with 
# Mixnoise (percent of noisy locations, 0=0%, 1=20%, 2=40%, 3=60%) ##Mixnoise_number##
# and scenario (see the paper for a description) ##scenario_number##
# Results_CARD/Mixnoise##Mixnoise_number##/Sc##scenario_number##/

source("../HGTfunctions3.R")

###############################
# Prep the ZI-HGT + CARD Data #
###############################

mixnoise <- 0:3
replicate <- 1:100
scenarios <- 1:5
alpha0 <- c(-0.1, 0.0, 0.1)  #alpha_0 is actually xover.75minusx(invlogit(..)) of these four values
alpha1 <- c(0.1, 0.3, 0.5)
grid <- expand.grid(alpha0, alpha1, replicate, scenarios, mixnoise)
colnames(grid) <- c("alpha0", "alpha1", "replicate", "scenarios", "mixnoise")
grid$ntotal <- 10

# Collect all of the results
res_df <- data.frame()
for (j in 1:18000){
  file_name <- paste0("./Results_CARD/Mixnoise", grid$mixnoise[j],
                      "/Sc", grid$scenarios[j], 
                      "/HGTsim.pseudo.MOB.n10.cellType6.Mixnoise", grid$mixnoise[j],
                      ".rep", grid$replicate[j], 
                      ".alpha0_", grid$alpha0[j], 
                      ".alpha1_", grid$alpha1[j],  ".rds")
  res <- readRDS(file_name)
  tmp_df <- data.frame(cell_count = grid$ntotal[j], scenario = grid$scenarios[j],
                       mixnoise = grid$mixnoise[j], rep = grid$replicate[j], 
                       alpha0 = grid$alpha0[j], alpha1 = grid$alpha1[j], 
                       CARD_rmse = res$rmse$CARD_rmse, Ave_rmse = res$rmse$Ave_rmse, 
                       j = j, WAIC = res$WAIC,
                       PercentMeanDiff = 100*(1-(res$rmse$Ave_rmse / res$rmse$CARD_rmse)))
  res_df <- rbind(res_df, tmp_df)
  if (j %% 100 == 0){
    gc()
    cat(j)
  }
}

# We need to select the best hyperparameters based on the WAIC
mixnoise <- 0:3
scenario <- 1:5
replicate <- 1:100
grid <- expand.grid(replicate, scenario, mixnoise)
colnames(grid) <- c("replicate", "scenario", "mixnoise")

Seq.Depth.res <- data.frame("Mixnoise"=grid$mixnoise, "Scenario"=grid$scenario, "Replicate"=grid$replicate,
                            "CARD_RMSE"=NA, "Ave_RMSE"=NA, 
                            "Ave_Percent_CARD"=NA, 
                            "Ave_Percent_Reduction"=NA, 
                            "alpha0"=NA, "alpha1"=NA, "WAIC"=NA)

for (j in 1:nrow(Seq.Depth.res)){
  res_df_sub <- res_df[res_df$scenario == grid$scenario[j] & res_df$mixnoise == grid$mixnoise[j] & res_df$rep == grid$replicate[j],]
  res_df_WAIC <- res_df_sub[res_df_sub$WAIC == min(res_df_sub$WAIC), ]
  Seq.Depth.res$CARD_RMSE[j] <- res_df_WAIC$CARD_rmse
  Seq.Depth.res$Ave_RMSE[j] <- res_df_WAIC$Ave_rmse
  Seq.Depth.res$Ave_Percent_CARD[j] <- (res_df_WAIC$Ave_rmse / res_df_WAIC$CARD_rmse) * 100
  Seq.Depth.res$Ave_Percent_Reduction[j] <- (1 - (res_df_WAIC$Ave_rmse / res_df_WAIC$CARD_rmse)) * 100
  Seq.Depth.res$alpha0[j] <- res_df_WAIC$alpha0
  Seq.Depth.res$alpha1[j] <- res_df_WAIC$alpha1
  Seq.Depth.res$WAIC[j] <- res_df_WAIC$WAIC
  if(j %% 100 == 0){
    cat(j, " ")
    gc()
  }
}

write.table(Seq.Depth.res, "HGTpCARD_SeqDepthSims_res.txt", row.names=F)

#######################################################
# Prep the CARD with different sequencing depths data #
#######################################################
# We are assuming your CARD simulation results are in a directory Results_CARD_OG

res <- data.frame()
reslist <- dir("Results_CARD_OG/")

for (i in 1:length(reslist)){
  tmplist <- readRDS(paste0("Results_CARD_OG/", reslist[i]))
  tmp <- rbind(tmplist[[1]], tmplist[[2]], tmplist[[3]], tmplist[[4]], tmplist[[5]])
  res <- rbind(res, tmp)
}

#How much did things change for 9, 10, 11 cells for each scenario across the replicates
for (j in 1:nrow(res)){
  comp <- res[res$Ncells == 10 & res$Replicate == res$Replicate[j] & res$Scenario == res$Scenario[j] & res$Mixnoise == res$Mixnoise[j], ]
  res$RMSEcomp10[j] <- (1 - (res$RMSE[j] / comp$RMSE))*100
}

# We decided not to use 12 cells, so remove such data
res <- res[res$Ncells < 12,]
write.table(res, "SeqDepthSimsCARDres5.txt", row.names=F)


#######################################
# Combine the data and make the plots #
#######################################

#Boxplot with RMSE reduction vs Mixnoise for CARD at different seq depths and HGT

library(dplyr)
library(ggplot2)
library(latex2exp)

boxplotdf <- read.table("../Simulations/Simulation_Results/SeqDepthSimsCARDres5.txt", header=T)
# We decided not to use 12 cells, so remove such data
boxplotdf <- boxplotdf[boxplotdf$Ncells < 12,]
boxplotdf$`% Reduction` <- boxplotdf$RMSEcomp10
boxplotdf$SeqDepthRelativeToBaseline <- recode(boxplotdf$Ncells, "9" = "- 10%", "10" = "Baseline", "11" = "+ 10%")
boxplotdf <- boxplotdf[, c("Replicate","RMSE","Scenario","Mixnoise","% Reduction","SeqDepthRelativeToBaseline")]

HGTdf <- read.table("../Simulations/Simulation_Results/HGTpCARD_SeqDepthSims_res.txt", header=T)
HGTdf_mean <- HGTdf[, c("Replicate","Ave_RMSE","Scenario","Mixnoise","Ave_Percent_Reduction")]
HGTdf_mean$SeqDepthRelativeToBaseline <- "HGT (mean)"
colnames(HGTdf_mean)[c(2,5)] <- c("RMSE","% Reduction")

boxplotdf <- rbind(boxplotdf, HGTdf_mean)

boxplotdf$Noisiness <- recode(boxplotdf$Mixnoise, "0" = "0%", "1" = "20%", "2" = "40%", "3" = "60%")
boxplotdf$NoisinessF <- as.factor(boxplotdf$Noisiness)
boxplotdf$SeqDepthRelativeToBaseline <- recode(boxplotdf$SeqDepthRelativeToBaseline, "Baseline" = "CARD")
boxplotdf$SeqDepthRelativeToBaselineF <- factor(boxplotdf$SeqDepthRelativeToBaseline, levels = c("HGT (mean)", "- 10%", "CARD", "+ 10%"))
boxplotdf$ScenarioF <- as.factor(boxplotdf$Scenario)

boxplotdf <- boxplotdf[boxplotdf$SeqDepthRelativeToBaseline != "CARD",]
boxplotdf <- boxplotdf[boxplotdf$SeqDepthRelativeToBaseline != "CARD",]
boxplotdf$SeqDepthRelativeToBaseline[boxplotdf$SeqDepthRelativeToBaseline == "HGT (mean)"] <- "ZI-HGT + CARD"
boxplotdf$SeqDepthRelativeToBaselineF <- factor(boxplotdf$SeqDepthRelativeToBaseline, levels = c("ZI-HGT + CARD", "- 10%", "+ 10%"))

boxplotdf_0 <- boxplotdf[boxplotdf$Mixnoise == 0,]
boxplotdf_1 <- boxplotdf[boxplotdf$Mixnoise == 1,]
boxplotdf_2 <- boxplotdf[boxplotdf$Mixnoise == 2,]
boxplotdf_3 <- boxplotdf[boxplotdf$Mixnoise == 3,]


g0 <- ggplot(boxplotdf_0, aes(x = ScenarioF, y = `% Reduction`, fill = SeqDepthRelativeToBaselineF)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="bottom") +
  labs(x = "Simulation Scenario", title = "RMSE Reduction on CARD Simulations", subtitle = TeX("$p_n = 0\\%"), fill = "Method/Sequencing Depth", y = "% RMSE Reduction")
g0
g1 <- ggplot(boxplotdf_1, aes(x = ScenarioF, y = `% Reduction`, fill = SeqDepthRelativeToBaselineF)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="bottom") +
  labs(x = "Simulation Scenario", title = "RMSE Reduction on CARD Simulations", subtitle = TeX("$p_n = 20\\%"), fill = "Method/Sequencing Depth", y = "% RMSE Reduction")
g1
g2 <- ggplot(boxplotdf_2, aes(x = ScenarioF, y = `% Reduction`, fill = SeqDepthRelativeToBaselineF)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="bottom") +
  labs(x = "Simulation Scenario", title = "RMSE Reduction on CARD Simulations", subtitle = TeX("$p_n = 40\\%"), fill = "Method/Sequencing Depth", y = "% RMSE Reduction")
g2
g3 <- ggplot(boxplotdf_3, aes(x = ScenarioF, y = `% Reduction`, fill = SeqDepthRelativeToBaselineF)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.position="bottom") +
  labs(x = "Simulation Scenario", title = "RMSE Reduction on CARD Simulations", subtitle = TeX("$p_n = 60\\%"), fill = "Method/Sequencing Depth", y = "% RMSE Reduction")
g3

ggsave("../Figures/Figure4.png", g0, units = "in", width = 8)
ggsave("../Figures/Supplementary_Figures/FigureS4.png", g1, units = "in", width = 8)
ggsave("../Figures/Supplementary_Figures/FigureS5.png", g2, units = "in", width = 8)
ggsave("../Figures/Supplementary_Figures/FigureS6.png", g3, units = "in", width = 8)


