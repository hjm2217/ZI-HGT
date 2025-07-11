############################################################################################################################
## Simulating the OSCC ST data based on SPARSim scRNAseq data with different numbers of cell types based on OSCC Sample 2 ##
############################################################################################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

source("../Utilities/OSCC_Sims_functions.R")

library(CARD)
library(Seurat)
library(dplyr)
library(pbmcapply)

n_cell_types <- c(12,10,8,6) #14 is already done in the main simulations
cell_counts <- 10
library_factor <- 0.05
Phi_factor <- 50
sample <- 2
SC_replicate <- 1:10
ST_replicate <- 1:10

# Set up the grid of parameters
grid <- expand.grid(n_cell_types, cell_counts, library_factor, Phi_factor, SC_replicate, ST_replicate, sample)
colnames(grid) <- c("n_cell_types", "cell_count", "library_factor", "Phi_factor", "SC_rep", "ST_rep", "sample")

grid$seed <- 1000*1101:1500

# Set up the cell types to use
# I randomly chose which cell types to remove
ct.select.list <- list()
ct.select.list[[1]] <- c("myofibroblast", "cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast",
               "detox-iCAF", "macrophage", "endothelial", "dendritic ",
               "conventional CD4+ T-helper cells", "cytotoxic CD8+ T ", "Tregs")
ct.select.list[[2]] <- c("cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast",
               "macrophage", "endothelial", "dendritic ",
               "conventional CD4+ T-helper cells", "cytotoxic CD8+ T ", "Tregs")
ct.select.list[[3]] <- c("cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast",
               "detox-iCAF", "macrophage", "dendritic ",
               "cytotoxic CD8+ T ")
ct.select.list[[4]] <- c("cancer cell", "ecm-myCAF", "Intermediate fibroblast",
               "macrophage", 
               "cytotoxic CD8+ T ", "Tregs")
	       



# Start the loop, simulate 40 datasets per array
for (j in ((job.id - 1) * 40 + 1) : (job.id * 40)){
  

  ##################################################################################
  # Load in the spatial pattern and pathologist annotations the sample in question #
  ##################################################################################
  pattern_patho_label <- readRDS(paste0("./Reference/pattern_patho_label_", grid$sample[j], ".rds"))
  location <- readRDS(paste0("./Reference/location_", grid$sample[j], ".rds"))

  ct.varname <- "cellype_fine"
  if (grid$n_cell_types[j] == 12){
     ct.select <- ct.select.list[[1]]
  } else if (grid$n_cell_types[j] == 10){
     ct.select <- ct.select.list[[2]]
  } else if (grid$n_cell_types[j] == 8){
     ct.select <- ct.select.list[[3]]
  } else if (grid$n_cell_types[j] == 6){
     ct.select <- ct.select.list[[4]]
  }

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
  saveRDS(sim, paste0("./SPARSim_ST_Data/Sample_", grid$sample[j], "/SPARSim_ST_n", grid$cell_count[j], "_cells_n", grid$n_cell_types[j], "_cell_types",
                      "_library_factor_", grid$library_factor[j], "_Phi_factor_", grid$Phi_factor[j],
		      "_scRNAseq_rep_", grid$SC_rep[j], "_ST_replicate_", grid$ST_rep[j], ".rds"))
  gc()
  
}












