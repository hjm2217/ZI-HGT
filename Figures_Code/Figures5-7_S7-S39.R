#############################################
## Plotting output from real OSCC Analysis ##
#############################################

# This code generates cell-type proportions estimates for CARD and the ZI-HGT + CARD (plotted
# in a raw and binary fashion) and 95% pointwise BCI plots for the ZI-HGT + CARD for all 12 samples.
# Note that only sample 1 is present in the main manuscript (Figures 5-7), 
# and samples 2-11 are in the supplementary (Supplementary Figures 7-39).

library(CARD)
library(dplyr)
library(ggplot2)
library(cowplot)

source("../Utilities/HGTfunctions3.R")

##############################################################
# Choose hyperparameters based on WAIC for the ZI-HGT + CARD #
##############################################################

# Run this code to generate the min_WAIC dataframe from the parent directory of where you
# have stored the results for the ZI-HGT + CARD real data analysis with *all* hyperparameters.
# We will only include the results with the WAIC-chosen hyperparameters for sample 1 on Github for storage reasons (the remaining results are on our OSF).
# If you want to skip this step and just make the plots, scroll down to 
###########################################################################################
# Make plots of cell locations and proportions for CARD, HGT + CARD, and the HGT 95% BCIs #
###########################################################################################


dataset <- 1:12
alpha0 <- c(-0.1, 0.0, 0.1, 0.2)  
alpha1 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
grid <- expand.grid(alpha0, alpha1, dataset)
colnames(grid) <- c("alpha0", "alpha1", "dataset")
grid$WAIC <- NA

for (j in 1:nrow(grid)){
  tmp <- readRDS(paste0("../Real_Data_Analysis/Results/HGT_CARD_obj_", grid$dataset[j], "_alpha0_", grid[j, "alpha0"], "_alpha1_", grid[j, "alpha1"],  ".rds"))
  grid$WAIC[j] <- tmp$WAIC
  if (j %% 10 == 0){
    cat(j, "\t")
    gc()
  }
}

min_WAIC <- grid %>% group_by(dataset) %>% slice_min(WAIC) %>% as.data.frame()

write.table(min_WAIC, "../Real_Data_Analysis/Results/OSCC_min_WAIC_output.txt", row.names=F)

###########################################################################################
# Make plots of cell locations and proportions for CARD, HGT + CARD, and the HGT 95% BCIs #
###########################################################################################
min_WAIC <- read.table("../Real_Data_Analysis/Results/OSCC_min_WAIC_output.txt", header=T)

for (j in 1:nrow(min_WAIC)){
  CARD_output <- readRDS(paste0("../Real_Data_Analysis/Results/CARD_obj_", min_WAIC$dataset[j], ".rds"))
  HGT_CARD_output <- readRDS(paste0("../Real_Data_Analysis/Results/HGT_CARD_obj_", min_WAIC$dataset[j], "_alpha0_", min_WAIC[j, "alpha0"], "_alpha1_", min_WAIC[j, "alpha1"],  ".rds"))
  
  proportions <- CARD_output@Proportion_CARD
  HGT_proportions <- as.data.frame(HGT_CARD_output$cell_type_proportion_matrices$mean_ct_prop)
  HGT_prop_2.5 <- as.data.frame(HGT_CARD_output$cell_type_proportion_matrices$lower_95_prop)
  HGT_prop_97.5 <- as.data.frame(HGT_CARD_output$cell_type_proportion_matrices$upper_95_prop)
  cell_type_names <- colnames(proportions)
  
  # Make sure to save the figures to the correct locations/names
  if (min_WAIC$dataset[j] == 1){
    figure_names <- paste0("../Figures/Figure", c("5a", "6a", "7a", "7b"), ".png")
  } else {
    figure_names <- paste0("../Figures/Supplementary_Figures/FigureS", 
                           c( (min_WAIC$dataset[j]-1)*2 + 5, (min_WAIC$dataset[j]-1)*2 + 6, min_WAIC$dataset[j] + 27, min_WAIC$dataset[j] + 27),
                           c("a", "a", "a", "b"), ".png")
  } # the above is a rather confusing map, we apologize for that. The figures will be named correctly though, don't worry

  
  p2 <- HGT.CARD.visualize.prop(
    proportion = proportions,        
    spatial_location = CARD_output@spatial_location, 
    ct.visualize = cell_type_names,                 
    colors = c("mediumpurple4", "white", "darkgreen"),
    NumCols = 7)                             
  ggsave(figure_names[2], p2, device="png", width = 24, height = 12, units = "in")
  
  h2 <- HGT.CARD.visualize.prop(
    proportion = HGT_proportions,      
    spatial_location = CARD_output@spatial_location, 
    ct.visualize = cell_type_names,                 ### selected cell types to visualize
    colors = c("mediumpurple4", "white", "darkgreen"), ### if not provide, we will use the default colors
    NumCols = 7)                             ### point size in ggplot2 scatterplot  
  ggsave(figure_names[1], h2, device="png", width = 24, height = 12, units = "in")
  
  h3 <- HGT.CARD.visualize.prop(
    proportion = HGT_prop_2.5,      
    spatial_location = CARD_output@spatial_location, 
    ct.visualize = cell_type_names,                 ### selected cell types to visualize
    colors = c("mediumpurple4", "white", "darkgreen"), ### if not provide, we will use the default colors
    NumCols = 7)                             ### point size in ggplot2 scatterplot  
  ggsave(figure_names[3], h3, device="png", width = 24, height = 12, units = "in")
  
  h4 <- HGT.CARD.visualize.prop(
    proportion = HGT_prop_97.5,      
    spatial_location = CARD_output@spatial_location, 
    ct.visualize = cell_type_names,                 ### selected cell types to visualize
    colors = c("mediumpurple4", "white", "darkgreen"), ### if not provide, we will use the default colors
    NumCols = 7)                             ### point size in ggplot2 scatterplot  
  ggsave(figure_names[4], h4, device="png", width = 24, height = 12, units = "in")
  
  cat(j, "\n")
}

###############################################################################
# Plots of two colors, one if cell-type proportion > cutoff (0.05) one if not #
###############################################################################
for (j in 1:nrow(min_WAIC)){
  CARD_output <- readRDS(paste0("../Real_Data_Analysis/Results/CARD_obj_", min_WAIC$dataset[j], ".rds"))
  HGT_CARD_output <- readRDS(paste0("../Real_Data_Analysis/Results/HGT_CARD_obj_", min_WAIC$dataset[j], "_alpha0_", min_WAIC[j, "alpha0"], "_alpha1_", min_WAIC[j, "alpha1"],  ".rds"))
  
  proportions <- CARD_output@Proportion_CARD
  HGT_proportions <- as.data.frame(HGT_CARD_output$cell_type_proportion_matrices$mean_ct_prop)
  HGT_prop_2.5 <- as.data.frame(HGT_CARD_output$cell_type_proportion_matrices$lower_95_prop)
  HGT_prop_97.5 <- as.data.frame(HGT_CARD_output$cell_type_proportion_matrices$upper_95_prop)
  cell_type_names <- colnames(proportions)
  
  # Make sure to save the figures to the correct locations/names
  if (min_WAIC$dataset[j] == 1){
    figure_names <- paste0("../Figures/Figure", c("5b", "6b"), ".png")
  } else {
    figure_names <- paste0("../Figures/Supplementary_Figures/FigureS", 
                           c( (min_WAIC$dataset[j]-1)*2 + 5, (min_WAIC$dataset[j]-1)*2 + 6),
                           c("b", "b"), ".png")
  } # the above is a rather confusing map, we apologize for that. The figures will be named correctly though, don't worry
  
  
  p1 <- HGT.CARD.two.tone(
    proportion = proportions,        
    spatial_location = CARD_output@spatial_location, 
    ct.visualize = cell_type_names,                 
    colors = c("grey70","darkgreen"),
    cutoff = 0.05,
    NumCols = 7,
    pointSize = 3.0)                             
  ggsave(figure_names[2], p1, device="png", width = 24, height = 12, units = "in")
 
  h1 <- HGT.CARD.two.tone(
    proportion = HGT_proportions,      
    spatial_location = CARD_output@spatial_location, 
    ct.visualize = cell_type_names,      
    colors = c("grey70","darkgreen"),
    cutoff = 0.05,
    NumCols = 7,
    pointSize = 3.0)                              
  ggsave(figure_names[1], h1, device="png", width = 24, height = 12, units = "in")

  cat(j, "\n")
}
 
  