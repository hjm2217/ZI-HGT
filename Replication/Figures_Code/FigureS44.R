################################################################################################
# Make plots of variability from both sources (CARD & ZI-HGT) for all cell types and locations #
################################################################################################
min_WAIC <- read.table("../Real_Data_Analysis/Results/OSCC_min_WAIC_output.txt", header=T)

variability <- data.frame(dataset = NA, alpha0 = NA, alpha1 = NA, CARD_var = NA, HGT_var = NA, cell_type = NA)

for (j in 1:nrow(min_WAIC)){
  
  HGT_CARD_output <- readRDS(paste0("../Real_Data_Analysis/Results/HGT_CARD_obj_", min_WAIC$dataset[j], "_alpha0_", min_WAIC[j, "alpha0"], "_alpha1_", min_WAIC[j, "alpha1"],  ".rds"))
  
  CARD_var <- HGT_CARD_output$cell_type_proportion_matrices$CARD_var
  HGT_var <- HGT_CARD_output$cell_type_proportion_matrices$HGT_var
  rownames(HGT_var) <- rownames(CARD_var); colnames(HGT_var) <- colnames(CARD_var)
  cell_type_names <- colnames(CARD_var)
  
  tmp_df <- data.frame(dataset = min_WAIC$dataset[j],
                       alpha0 = min_WAIC[j, "alpha0"],
                       alpha1 = min_WAIC[j, "alpha1"],
                       CARD_var = as.numeric(CARD_var),
                       HGT_var = as.numeric(HGT_var),
                       cell_type = rep(cell_type_names, each = nrow(CARD_var)))
  
  variability <- rbind(variability, tmp_df)
  gc()
  cat(j, "\n")
}

variability <- variability[-1,]  #remove first row of NAs
variability$cell_type <- factor(variability$cell_type, levels = cell_type_names)
variability$dataset <- factor(variability$dataset, levels = 1:12)

library(tidyr)
variability_long <- variability %>%
  pivot_longer(cols = c("CARD_var", "HGT_var"), 
               names_to = "Variability_Source", 
               values_to = "Variability") %>%
  mutate(Variability_Source = recode(Variability_Source, 
                                     CARD_var = "CARD", 
                                     HGT_var = "ZI-HGT"))

variability_summary <- variability_long %>%
  group_by(dataset, Variability_Source) %>%
  summarise(mean_variability = mean(Variability, na.rm = TRUE)) %>%
  mutate(prop_variability = mean_variability / sum(mean_variability))

variability_double_summary <- variability_summary %>%
  group_by(Variability_Source) %>%
  summarise(mean_prop = mean(prop_variability)) 

# Make a boxplot of the variability for each sample/type of variability, facet_wrapped by cell type
p <- ggplot(variability_long, aes(x = dataset, y = Variability, color = Variability_Source)) +
  ggplot2::geom_boxplot(outliers = FALSE, fill = "white") +
  facet_wrap(~cell_type, scales = "free_y", ncol = 7) +
  scale_color_manual(values = c("mediumpurple", "green2")) +
  labs(x = "Sample", y = "Variability", fill = "Source of Variability", color = "Source of Variability") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
ggsave("../Figures/Supplementary_Figures/FigureS44.png", p, device="png", width = 24, height = 10, units = "in")
