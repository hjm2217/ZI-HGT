#############################################
## Plotting output from real OSCC Analysis ##
#############################################

# This code generates cell-type proportions estimates for CARD and the ZI-HGT + CARD (plotted
# in a raw and binary fashion) and 95% pointwise BCI plots for the ZI-HGT + CARD for all 12 samples.
# Note that only sample 1's results for the ZI-HGT + CARD are present in the main manuscript (Figure 5), 
# and the rest of sample 1 plus samples 2-11 are in the supplementary (Supplementary Figures 7-41).

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
# If you follow our file structure, this is the Real_Data_Analysis directory.
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
    figure_names <- c("../Figures/Figure5a.png", "../Figures/Supplementary_Figures/FigureS8a.png", "../Figures/Figure6a.png", "../Figures/Figure6b.png")
    #figure_names <- paste0("../Figures/Figure", c("5a", "6a", "7a", "7b"), ".png")
  } else {
    figure_names <- paste0("../Figures/Supplementary_Figures/FigureS", 
                           c( (min_WAIC$dataset[j]-1)*2 + 7, (min_WAIC$dataset[j]-1)*2 + 8, min_WAIC$dataset[j] + 29, min_WAIC$dataset[j] + 29),
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
    figure_names <- paste0("../Figures/", c("Figure5b", "Supplementary_Figures/FigureS8b"), ".png")
  } else {
    figure_names <- paste0("../Figures/Supplementary_Figures/FigureS", 
                           c( (min_WAIC$dataset[j]-1)*2 + 7, (min_WAIC$dataset[j]-1)*2 + 8),
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
 

#################################################################################################
# Compare the mean cell-type proportions of the ZI-HGT + CARD & CARD to the scRNA-seq reference #
#################################################################################################

# Load the settings that produced the minimum WAIC
min_WAIC <- read.table("../Real_Data_Analysis/Results/OSCC_min_WAIC_output.txt", header=T)

# Load the single-cell reference dataset
load(file = "Data/puram_data.Robj")
sc_meta <- puram_data@meta.data

sc_cell_type_mean_props <- sc_meta %>% 
  group_by(cellype_fine) %>% 
  summarise(n = n()) %>% 
  mutate(mean_prop = n / sum(n)) 

# Set up comparisons for sample 1
CARD_output <- readRDS(paste0("./Results/CARD_obj_", min_WAIC$dataset[1], ".rds"))
CARD_prop <- CAD_output@Proportion_CARD
HGT_CARD_output <- readRDS(paste0("./Results/HGT_CARD_obj_", min_WAIC$dataset[1], "_alpha0_", min_WAIC[1, "alpha0"], "_alpha1_", min_WAIC[1, "alpha1"],  ".rds"))
HGT_prop <- HGT_CARD_output$cell_type_proportion_matrices$mean_ct_prop

# Get the mean proportions for each cell type
CARD_prop_mean <- colMeans(CARD_prop)
HGT_prop_mean <- colMeans(HGT_prop)

# Append to the sc_cell_type_mean_props df
cell_type_mean_props <- sc_cell_type_mean_props %>% 
  mutate(ZI_HGT = HGT_prop_mean, CARD = CARD_prop_mean) %>% 
  select(cellype_fine, mean_prop, ZI_HGT, CARD) %>% 
  arrange(desc(mean_prop))

  