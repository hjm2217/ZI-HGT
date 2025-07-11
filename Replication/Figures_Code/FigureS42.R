########################################################
## Plotting OSCC Sample 1 Read Counts & where they go ##
########################################################

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

# Convert the matrix to a data frame
df_t_counts <- data.frame(value = as.vector(t_s1_counts))
df_counts <- data.frame(value = as.vector(sample_1_counts))

# Specify custom limits so that we can actually see the shape of the data and aren't super zoomed out
custom_limits <- c(-10, 10)

# Make a factor to color the different values in the ST data
df_counts$value_factor <- as.character(df_counts$value)
df_counts$value_factor[df_counts$value >= 5] <- ">=5"
df_counts$value_factor <- factor(df_counts$value_factor, levels = c("0", "1", "2", "3", "4", ">=5"))

df_t_counts$untransformed_value_factor <- df_counts$value_factor
df_t_counts_0 <- df_t_counts[df_t_counts$untransformed_value_factor == "0",]
df_t_counts_1 <- df_t_counts[df_t_counts$untransformed_value_factor == "1",]
df_t_counts_2 <- df_t_counts[df_t_counts$untransformed_value_factor == "2",]
df_t_counts_3 <- df_t_counts[df_t_counts$untransformed_value_factor == "3",]
df_t_counts_4 <- df_t_counts[df_t_counts$untransformed_value_factor == "4",]
df_t_counts_5_or_more <- df_t_counts[df_t_counts$untransformed_value_factor == ">=5",]

# Create histograms
g <- ggplot(df_counts, aes(x = value, fill = value_factor)) +
  geom_histogram(binwidth = 0.5, color = "black") +
  labs(title = "Raw Read Counts", subtitle = "OSCC Sample 1",
       x = "Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("0" = "#4DAF4A", "1" = "#377EB8", "2" = "#FF7F00", "3" = "#984EA3", "4" = "#E41A1C", ">=5" = "#FFFF33"))
g
ggsave("../Figures/Supplementary_Figures/FigureS42a.png", g, device="png", width = 5, height = 5, units = "in")

h0 <- ggplot(df_t_counts_0, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "#4DAF4A", color = "black") +
  labs(title = "Transformed Read Counts", subtitle = "Raw Count = 0",
       x = "Transformed Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("../Figures/Supplementary_Figures/FigureS42b.png", h0, device="png", width = 3, height = 3, units = "in")

h1 <- ggplot(df_t_counts_1, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "#377EB8", color = "black") +
  labs(title = "Transformed Read Counts", subtitle = "Raw Count = 1",
       x = "Transformed Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("../Figures/Supplementary_Figures/FigureS42c.png", h1, device="png", width = 3, height = 3, units = "in")

h2 <- ggplot(df_t_counts_2, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "#FF7F00", color = "black") +
  labs(title = "Transformed Read Counts", subtitle = "Raw Count = 2",
       x = "Transformed Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("../Figures/Supplementary_Figures/FigureS42d.png", h2, device="png", width = 3, height = 3, units = "in")

h3 <- ggplot(df_t_counts_3, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "#984EA3", color = "black") +
  labs(title = "Transformed Read Counts", subtitle = "Raw Count = 3",
       x = "Transformed Counts",y = "Frequency") +
  theme_minimal() +
  xlim(custom_limits) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("../Figures/Supplementary_Figures/FigureS42e.png", h3, device="png", width = 3, height = 3, units = "in")

h4 <- ggplot(df_t_counts_4, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "#E41A1C", color = "black") +
  labs(title = "Transformed Read Counts", subtitle = "Raw Count = 4") +
  theme_minimal() +
  xlim(custom_limits) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("../Figures/Supplementary_Figures/FigureS42f.png", h4, device="png", width = 3, height = 3, units = "in")

h5 <- ggplot(df_t_counts_5_or_more, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "#FFFF33", color = "black") +
  labs(title = "Transformed Read Counts", subtitle = "Raw Count >= 5",
       x = "Transformed Counts",y = "Frequency") +
  theme_minimal() +
  xlim(c(-10, 20)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave("../Figures/Supplementary_Figures/FigureS42g.png", h5, device="png", width = 3, height = 3, units = "in")


# These figures are then pasted together with Adobe Illustrator



