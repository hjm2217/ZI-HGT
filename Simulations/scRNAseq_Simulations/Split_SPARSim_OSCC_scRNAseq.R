############################################
## Splitting the simulated  scRNAseq data ##
############################################

# The first split will be used to construct the simulated ST data
# The second split will be used to construct the reference basis matrix to run the deconvolution

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

library(Seurat)

# Set up the grid
lib_factor <- c(0.2, 0.01, 0.4, 0.6)
Phi_factor <- 4
SC_rep <- 1:10
grid <- expand.grid(lib_factor, Phi_factor, SC_rep)
colnames(grid) <- c("lib_factor", "Phi_factor", "SC_rep")
grid$seed <- rep(2000 * 25:28, 10) 



# Split 4 datasets on each array

for (j in ((job.id-1) * 4 + 1):(job.id*4)){

    # Load in the real and simulated scRNAseq data
    puram_data <- readRDS("../Reference/puram_data.rds")
    sim_sc_data <- readRDS(paste0("scRNAseq_sim_data/sim_library_", grid$lib_factor[j],
				  "_Phi_", grid$Phi_factor[j],
				  "_rep_", grid$SC_rep[j], ".rds"))

    # Get the simulated scRNAseq data counts matrix and add it as an assay to the real scRNAseq data (this is to make it easier to split and to get the meta.info)
    sc_matrix <- sim_sc_data$count_matrix
    sim_assay <- CreateAssayObject(counts = sc_matrix)
    puram_data[["simSC"]] <- sim_assay
    
    # Split the scRNAseq data into two pieces
    set.seed(grid$seed[j])
    puram_data@meta.data$split <- sample(c(1,2), nrow(puram_data@meta.data), replace=T)
    sim_splits <- SplitObject(puram_data, split.by = "split")
    names(sim_splits) <- c("Simulation", "Deconvolution")
    saveRDS(sim_splits, paste0("scRNAseq_sim_data/split_scRNAseq_sim_library_", grid$lib_factor[j],
    		                  "_Phi_", grid$Phi_factor[j],        	
                                  "_rep_", grid$SC_rep[j], ".rds"))

    gc()
    cat(j, "\t")
}