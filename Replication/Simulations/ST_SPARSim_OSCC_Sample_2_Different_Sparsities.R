################################################################
## Simulating the OSCC ST data based on SPARSim scRNAseq data ##
################################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

source("../Utilities/OSCC_Sims_functions.R")

library(CARD)
library(Seurat)
library(dplyr)
library(pbmcapply)

# Only simulate the datasets that we want to simulate to produce the correct sparsities
cell_counts <- c(10)
library_factor <- c(0.05,0.1)
Phi_factor <- c(50, 20, 10, 4)
SC_replicate <- 1:10
ST_replicate <- 1:10

# Set up the grid of parameters
grid <- expand.grid(cell_counts, library_factor, Phi_factor, SC_replicate, ST_replicate)
colnames(grid) <- c("cell_count", "library_factor", "Phi_factor", "SC_rep", "ST_rep")

library_factor_2 <- c(0.2, 0.3)
Phi_factor_2 <- c(4)
grid2 <- expand.grid(cell_counts, library_factor_2, Phi_factor_2, SC_replicate, ST_replicate)
colnames(grid2) <- c("cell_count", "library_factor", "Phi_factor", "SC_rep", "ST_rep")

library_factor_3 <- 0.4
Phi_factor_3 <- c(4)
grid3 <- expand.grid(cell_counts, library_factor_3, Phi_factor_3, SC_replicate, ST_replicate)
colnames(grid3) <- c("cell_count", "library_factor", "Phi_factor", "SC_rep", "ST_rep")

library_factor_4 <- 0.6
Phi_factor_4 <- c(4)
grid4 <- expand.grid(cell_counts, library_factor_4, Phi_factor_4, SC_replicate, ST_replicate)
colnames(grid4) <- c("cell_count", "library_factor", "Phi_factor", "SC_rep", "ST_rep")

grid <- rbind(grid, grid2, grid3, grid4)
grid$seed <- 1000*1:1200
grid$sample <- 2


# Start the loop, simulate 10 datasets per array
for (j in ((job.id - 1) * 10 + 1) : (job.id * 10)){
  

  ##################################################################################
  # Load in the spatial pattern and pathologist annotations the sample in question #
  ##################################################################################
  pattern_patho_label <- readRDS(paste0("./Reference/pattern_patho_label_", grid$sample[j], ".rds"))
  location <- readRDS(paste0("./Reference/location_", grid$sample[j], ".rds"))

  ct.varname <- "cellype_fine"
  ct.select <- c("myofibroblast", "cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast",
               "detox-iCAF", "macrophage", "endothelial", "dendritic ", "mast",
               "conventional CD4+ T-helper cells", "cytotoxic CD8+ T ", "Tregs", "cytotoxic CD8+ T exhausted")
  sample.varname <- "orig.ident"

  sample <- readRDS(paste0("../Real_Data_Analysis/Data/sample_", grid$sample[j], ".rds"))
  sample_transcripts <- rownames(sample@assays$SCT$counts)


  ################################################
  # Load in the real and simulated scRNAseq data #
  ################################################

  split_sim_scRNAseq <- readRDS(paste0("./scRNAseq_Simulations/scRNAseq_sim_data/split_scRNAseq_sim_library_",
  		                       grid$library_factor[j], "_Phi_", grid$Phi_factor[j], "_rep_", grid$SC_rep[j], ".rds"))

  eset.sub.split1 <- split_sim_scRNAseq$Simulation

  # Only include genes in the simulated ST data that are in the real ST data
  sim_sc_match <- eset.sub.split1@assays$simSC$counts[rownames(eset.sub.split1@assays$simSC$counts) %in% sample_transcripts, ]

  # Then add this new matrix to the Seurat object
  sim_sc_assay_match <- CreateAssayObject(counts = sim_sc_match)
  eset.sub.split1[["simScMatch"]] <- sim_sc_assay_match
  
  ##################################
  # Generate the simulated ST data #
  ##################################
  cat("#######################################################################################################################################\n")
  cat("Generating the", grid$ST_rep[j], "th replicate of simulated ST data with", 
      grid$cell_count[j], "cells per location, library factor", grid$library_factor[j], ", Phi factor", grid$Phi_factor[j],
      ", based on scRNAseq replicate", grid$SC_rep[j], "for sample", grid$sample[j], ".\n")
  cat("#######################################################################################################################################\n")    

  sim <- generateSpatial_norep_fixedProp_OSCC(seed = grid$seed[j],
                                              eset.sub.split1 = eset.sub.split1,
                                              ct.varname = "cellype_fine",
                                              sample.varname = "orig.ident",
                                              ct.select = ct.select,
                                              RNA.counts = eset.sub.split1@assays[["simScMatch"]]@counts,
                                              pattern_gp_label = pattern_patho_label,
                                              ntotal = grid$cell_count[j])
  saveRDS(sim, paste0("./SPARSim_ST_Data/Sample_", grid$sample[j], "/SPARSim_ST_n", grid$cell_count[j], 
                      "_cells_library_factor_", grid$library_factor[j], "_Phi_factor_", grid$Phi_factor[j],
		      "_scRNAseq_rep_", grid$SC_rep[j], "_ST_replicate_", grid$ST_rep[j], ".rds"))
  gc()
  
}












