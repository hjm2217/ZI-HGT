#################################################################
## Simulating scRNAseq data for use in the OSCC ST simulations ##
#################################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

# Install from source (downloaded from GitLab page)
#install.packages("sparsim-master.tar.gz", repos = NULL, type = "source")
library(SPARSim)

# We want to get parameters based on our actual single cell data, then change the library size
# to change the sparsity. 

###########################################################
# Define the SPARSim Input Parameters from the Puram Data #
###########################################################

puram_data <- readRDS("../Reference/puram_data.rds")
counts <- puram_data@assays$RNA@counts #23686 transcripts x 5902 cells

counts_norm <- scran_normalization(counts)

# We need to find the indices of all of the different cell types
ct.select <- c("myofibroblast", "cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast",
               "detox-iCAF", "macrophage", "endothelial", "dendritic ", "mast",
               "conventional CD4+ T-helper cells", "cytotoxic CD8+ T ", "Tregs", "cytotoxic CD8+ T exhausted")
cell.types <- list()
for (i in 1:14){
  cell.types[[i]] <- list()
}
names(cell.types) <- ct.select

for (i in 1:14){
  cell.types[[i]][[1]] <- subset(puram_data@meta.data, cellype_fine == names(cell.types)[i])
  cell.types[[i]][[2]] <- rownames(cell.types[[i]][[1]])
  cell.types[[i]][[3]] <- which(colnames(counts) %in% cell.types[[i]][[2]])
  cat(names(cell.types)[i], length(cell.types[[i]][[3]]), "\n") #print how many cells we have of each cell type
}

counts_matrix_cell_types_index <- list(
  "myofibroblast" = cell.types[[1]][[3]],
  "cancer cell" = cell.types[[2]][[3]],
  "B cell" = cell.types[[3]][[3]], 
  "ecm-myCAF" = cell.types[[4]][[3]],
  "Intermediate fibroblast" = cell.types[[5]][[3]],
  "detox-iCAF" = cell.types[[6]][[3]],
  "macrophage" = cell.types[[7]][[3]],
  "endothelial" = cell.types[[8]][[3]],
  "dendritic " = cell.types[[9]][[3]], 
  "mast" = cell.types[[10]][[3]],
  "conventional CD4+ T-helper cells" = cell.types[[11]][[3]],
  "cytotoxic CD8+ T " = cell.types[[12]][[3]], 
  "Tregs" = cell.types[[13]][[3]], 
  "cytotoxic CD8+ T exhausted" = cell.types[[14]][[3]]
)

# Make the simulation parameters using the different cell types all at once
SPARSim_sim_param <- SPARSim_estimate_parameter_from_data(raw_data = as.matrix(counts), 
                                                          norm_data = counts_norm, 
                                                          conditions = counts_matrix_cell_types_index)


###########################################################
# Now change the library sizes and biologic variabilities #
# to induce different amounts of sparsity #
###########################################################

SPARSim_sim_param_library_0.1_Phi_10 <- SPARSim_sim_param
SPARSim_sim_param_library_0.1_Phi_20 <- SPARSim_sim_param
SPARSim_sim_param_library_0.1_Phi_50 <- SPARSim_sim_param
SPARSim_sim_param_library_0.1_Phi_100 <- SPARSim_sim_param
SPARSim_sim_param_library_0.05_Phi_10 <- SPARSim_sim_param
SPARSim_sim_param_library_0.05_Phi_20 <- SPARSim_sim_param
SPARSim_sim_param_library_0.05_Phi_50 <- SPARSim_sim_param
SPARSim_sim_param_library_0.05_Phi_100 <- SPARSim_sim_param
SPARSim_sim_param_library_0.03_Phi_10 <- SPARSim_sim_param
SPARSim_sim_param_library_0.03_Phi_20 <- SPARSim_sim_param
SPARSim_sim_param_library_0.03_Phi_50 <- SPARSim_sim_param
SPARSim_sim_param_library_0.03_Phi_100 <- SPARSim_sim_param
SPARSim_sim_param_library_0.01_Phi_10 <- SPARSim_sim_param
SPARSim_sim_param_library_0.01_Phi_20 <- SPARSim_sim_param
SPARSim_sim_param_library_0.01_Phi_50 <- SPARSim_sim_param
SPARSim_sim_param_library_0.01_Phi_100 <- SPARSim_sim_param

for (i in 1:14){

  SPARSim_sim_param_library_0.1_Phi_10[[i]]$lib_size <- 0.1 * SPARSim_sim_param_library_0.1_Phi_10[[i]]$lib_size
  SPARSim_sim_param_library_0.1_Phi_10[[i]]$variability <- 10 * SPARSim_sim_param_library_0.1_Phi_10[[i]]$variability
  SPARSim_sim_param_library_0.1_Phi_20[[i]]$lib_size <- 0.1 * SPARSim_sim_param_library_0.1_Phi_20[[i]]$lib_size
  SPARSim_sim_param_library_0.1_Phi_20[[i]]$variability <- 20 * SPARSim_sim_param_library_0.1_Phi_20[[i]]$variability
  SPARSim_sim_param_library_0.1_Phi_50[[i]]$lib_size <- 0.1 * SPARSim_sim_param_library_0.1_Phi_50[[i]]$lib_size
  SPARSim_sim_param_library_0.1_Phi_50[[i]]$variability <- 50 * SPARSim_sim_param_library_0.1_Phi_50[[i]]$variability
  SPARSim_sim_param_library_0.1_Phi_100[[i]]$lib_size <- 0.1 * SPARSim_sim_param_library_0.1_Phi_100[[i]]$lib_size
  SPARSim_sim_param_library_0.1_Phi_100[[i]]$variability <- 1000 * SPARSim_sim_param_library_0.1_Phi_100[[i]]$variability

  SPARSim_sim_param_library_0.05_Phi_10[[i]]$lib_size <- 0.05 * SPARSim_sim_param_library_0.05_Phi_10[[i]]$lib_size
  SPARSim_sim_param_library_0.05_Phi_10[[i]]$variability <- 10 * SPARSim_sim_param_library_0.05_Phi_10[[i]]$variability
  SPARSim_sim_param_library_0.05_Phi_20[[i]]$lib_size <- 0.05 * SPARSim_sim_param_library_0.05_Phi_20[[i]]$lib_size
  SPARSim_sim_param_library_0.05_Phi_20[[i]]$variability <- 20 * SPARSim_sim_param_library_0.05_Phi_20[[i]]$variability
  SPARSim_sim_param_library_0.05_Phi_50[[i]]$lib_size <- 0.05 * SPARSim_sim_param_library_0.05_Phi_50[[i]]$lib_size
  SPARSim_sim_param_library_0.05_Phi_50[[i]]$variability <- 50 * SPARSim_sim_param_library_0.05_Phi_50[[i]]$variability
  SPARSim_sim_param_library_0.05_Phi_100[[i]]$lib_size <- 0.05 * SPARSim_sim_param_library_0.05_Phi_100[[i]]$lib_size
  SPARSim_sim_param_library_0.05_Phi_100[[i]]$variability <- 1000 * SPARSim_sim_param_library_0.05_Phi_100[[i]]$variability

  
}

j <- job.id
  seed <- (j-1)*10000
  set.seed(seed)
  
  sim_library_0.1_Phi_10 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.1_Phi_10)
  sim_library_0.1_Phi_20 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.1_Phi_20)
  sim_library_0.1_Phi_50 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.1_Phi_50)
  sim_library_0.1_Phi_100 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.1_Phi_100)
  sim_library_0.05_Phi_10 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.05_Phi_10)
  sim_library_0.05_Phi_20 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.05_Phi_20)
  sim_library_0.05_Phi_50 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.05_Phi_50)
  sim_library_0.05_Phi_100 <- SPARSim_simulation(dataset_parameter = SPARSim_sim_param_library_0.05_Phi_100)
  
  
  saveRDS(sim_library_0.1_Phi_10, paste0("./scRNAseq_sim_data/sim_library_0.1_Phi_10_rep_", j,".rds"))
  saveRDS(sim_library_0.1_Phi_20, paste0("./scRNAseq_sim_data/sim_library_0.1_Phi_20_rep_", j,".rds"))
  saveRDS(sim_library_0.1_Phi_50, paste0("./scRNAseq_sim_data/sim_library_0.1_Phi_50_rep_", j,".rds"))
  saveRDS(sim_library_0.1_Phi_100, paste0("./scRNAseq_sim_data/sim_library_0.1_Phi_100_rep_", j,".rds"))
  saveRDS(sim_library_0.05_Phi_10, paste0("./scRNAseq_sim_data/sim_library_0.05_Phi_10_rep_", j,".rds"))
  saveRDS(sim_library_0.05_Phi_20, paste0("./scRNAseq_sim_data/sim_library_0.05_Phi_20_rep_", j,".rds"))
  saveRDS(sim_library_0.05_Phi_50, paste0("./scRNAseq_sim_data/sim_library_0.05_Phi_50_rep_", j,".rds"))
  saveRDS(sim_library_0.05_Phi_100, paste0("./scRNAseq_sim_data/sim_library_0.05_Phi_100_rep_", j,".rds"))

  cat(j)






