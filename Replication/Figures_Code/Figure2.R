########################################################################
## Plotting OSCC Sample 1 Read Counts before and after Transformation ##
########################################################################

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)

load("../Real_Data_Analysis/Data/sample_1.Robj")

sample_1_counts <- sample_1@assays$SCT@counts %>% as.matrix()
hist(sample_1_counts)

# Load all of the functions we may need
source("../Utilities/HGTfunctions3.R")

# Transform the counts (using WAIC-selected hyperparameters from real data analysis)
t_s1_counts <- transformX2H(sample_1_counts,
                            alpha0 = xover.75minusx(invlogit( 0 )),
                            alpha1 = 0.5,
                            seed = 123456)

# Convert the matrix to a data frame
df_t_counts <- data.frame(value = as.vector(t_s1_counts))
df_counts <- data.frame(value = as.vector(sample_1_counts))

# Specify custom limits so that we can actually see the shape of the data and aren't super zoomed out
custom_limits <- c(-10, 10)

# Create histograms
options(scipen = 999)
g <- ggplot(df_t_counts, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "green4", color = "black") +
  labs(title = "Distribution of Transformed Read Counts", subtitle = "OSCC Sample 1",
       x = "Transformed Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits)
g
ggsave("../Figures/Figure2b.png", g, device="png", width = 5, height = 5, units = "in")

h <- ggplot(df_counts, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "green4", color = "black") +
  labs(title = "Distribution of Raw Read Counts", subtitle = "OSCC Sample 1",
       x = "Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits)
h
ggsave("../Figures/Figure2a.png", h, device="png", width = 5, height = 5, units = "in")

