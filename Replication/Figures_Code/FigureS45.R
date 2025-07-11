################################################################
## Plotting impact of transforming and back-transforming data ##
################################################################
library(dplyr)
library(Seurat)
library(dplyr)
library(LaplacesDemon)
library(ggplot2)

####################
# Load in the data #
####################
res <- read.csv("../Real_Data_Analysis/Back_Transformations/mean_diff_transformation_and_back_sample1.csv")
res$sample <- 1
for (i in 2:12){
  temp <- read.csv(paste0("../Real_Data_Analysis/Back_Transformations/mean_diff_transformation_and_back_sample", i, ".csv"))
  temp$sample <- i
  res <- rbind(res, temp)
}
res$sample <- as.factor(res$sample)

########
# Plot #
########
g <- ggplot(res, aes(x = sample, mean_diff)) +
  geom_boxplot(fill = "green4") +
  labs(title = "Reconstruction Error", subtitle = "OSCC Data vs. Back-Transformed ZI-HGT Data",
       x = "Sample",y = "Mean Error") +
  theme_minimal() + ylim(-1, 0)
g
ggsave("../Figures/Supplementary_Figures/FigureS45.png", g, device="png", width = 6, height = 6, units = "in")

