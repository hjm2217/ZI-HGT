################################################################
## Plotting OSCC Sample 1 Read Counts Post CARD Normalization ##
################################################################

# CARD normalizes by row, what do the distributions of counts and transformed counts
# look like afterwards?

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

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

# Transpose so that rows are spots and columns are genes
sample_1_counts <- t(sample_1_counts)
t_s1_counts <- t(t_s1_counts)

# Normalize by location (row-normalize)
sample_1_counts <- sample_1_counts / rowSums(sample_1_counts)
t_s1_counts <- t_s1_counts / rowSums(t_s1_counts)

# Convert the matrix to a data frame
df_t_counts <- data.frame(value = as.vector(t_s1_counts))
df_counts <- data.frame(value = as.vector(sample_1_counts))

# Specify custom limits so that we can actually see the shape of the data and aren't super zoomed out
custom_limits <- c(-0.01, 0.01)

# Create histograms
g <- ggplot(df_counts, aes(x = value)) +
  geom_histogram(binwidth = 0.0003, color = "black", fill = "green4") +
  labs(title = "CARD-Normalized Raw Read Counts", subtitle = "OSCC Sample 1",
       x = "CARD-Normalized Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits)
g
ggsave("../Figures/Supplementary_Figures/FigureS43a.png", g, device="png", width = 5, height = 5, units = "in")

h <- ggplot(df_t_counts, aes(x = value)) +
  geom_histogram(binwidth = 0.0003, color = "black", fill = "green4") +
  labs(title = "CARD-Normalized Transformed Read Counts", subtitle = "OSCC Sample 1",
       x = "CARD-Normalized Transformed Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits)
h
ggsave("../Figures/Supplementary_Figures/FigureS43b.png", h, device="png", width = 5, height = 5, units = "in")


