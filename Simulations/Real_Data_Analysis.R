######################################################
## Running the ZI-HGT + CARD on the 12 OSCC Samples ##
######################################################

# We ran this analysis on our university's high performance computer, which uses the Slurm
# workload manager.  For the exact specifications of the job we ran, see the 
# Real_Data_Analysis.sh file.  You may need to adjust some of the options depending on the
# computational resources you have access to.

# Additionally, the data files are too large to comfortably store on GitHub.  You may find
# them on our OSF.


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

library(dplyr)
library(Seurat)
library(pbmcapply)
library(SingleCellExperiment)
library(CARD)

source("../Utilities/HGTfunctions3.R")  

#########################
# Load single-cell data #
#########################
load(file = "Data/puram_data.Robj")
sc_count <- puram_data@assays$RNA@counts
sc_meta <- puram_data@meta.data

###############
# Set up grid #
###############
dataset <- 1:12
alpha0 <- c(-0.1, 0.0, 0.1, 0.2)  #alpha_0 is actually xover.75minusx(invlogit(..)) of these four values
alpha1 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
grid <- expand.grid(alpha0, alpha1, dataset)
colnames(grid) <- c("alpha0", "alpha1", "dataset")
grid$seed <- 1000 * 1:240 + 1

###################################################################
# Run CARD with HGT across datasets and potential hyperparameters #
###################################################################

for (j in (job.id)){

    #########################
    # Load and prep ST data #
    #########################
    sample <- readRDS(paste0("Data/sample_", grid$dataset[j], ".rds"))
    sample@meta.data$sample <- NULL
    spatial_count <- sample@assays$SCT@counts
    
    spatial_location <- sample@images$tumor@coordinates %>% select(row, col)
    names(spatial_location) <- c("x","y")
    
    if (!identical(colnames(spatial_count), rownames(spatial_location))){
      #ensure that count matrix contains the same spots as the location matrix 
      spatial_count <- spatial_count[, colnames(spatial_count) %in% rownames(spatial_location)]
    }
    #match order
    spatial_location <- spatial_location[order(match(rownames(spatial_location), colnames(spatial_count))), , drop = FALSE]
    
    ##########################
    # Deconvolute using CARD #
    ##########################
    CARD_obj = createCARDObject(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = spatial_count,
      spatial_location = spatial_location,
      ct.varname = "cellype_fine",
      ct.select = unique(sc_meta$cellype_fine),
      sample.varname = "orig.ident",
      minCountGene = 100,  
      minCountSpot = 5) 
    
    CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
    
    proportions <- as.data.frame(CARD_obj@Proportion_CARD)
    proportions$dominant_celltype <- apply(proportions, 1, function(x) names(proportions)[which.max(x)])
    
    saveRDS(CARD_obj, paste0("Results/CARD_obj_", grid$dataset[j], ".rds"))
    write.table(proportions, paste0("Results/CARD_prop_", grid$dataset[j], ".rds"), row.names=F, sep="\t")

    cat("Done with original CARD.\n######################################################################################")

    #####################################
    # Deconvolute using the ZI-HGT+CARD #
    #####################################
    HGT_CARD_obj <- HGTplusCARD(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = spatial_count,
      spatial_location = spatial_location,
      ct.varname = "cellype_fine",
      ct.select = unique(sc_meta$cellype_fine),
      sample.varname = "orig.ident",
      minCountGene = 100,
      minCountSpot = 5,
      alpha0 = xover.75minusx(invlogit( grid[j, "alpha0"] )),
      alpha1 = grid[j, "alpha1"],
      n = 100,
      seed = grid$seed[j],
      doWAIC = TRUE)
    
    HGT_proportions <- as.data.frame(HGT_CARD_obj$cell_type_proportion_matrices$mean_ct_prop)
    HGT_proportions$dominant_celltype <- apply(HGT_proportions, 1, function(x) names(HGT_proportions)[which.max(x)])
    
    saveRDS(HGT_CARD_obj, paste0("Results/HGT_CARD_obj_", grid$dataset[j], "_alpha0_", grid[j, "alpha0"], "_alpha1_", grid[j, "alpha1"],  ".rds"))
    write.table(HGT_proportions, paste0("Results/HGT_CARD_prop_", grid$dataset[j], "_alpha0_", grid[j, "alpha0"], "_alpha1_", grid[j, "alpha1"],  ".rds"), row.names=F, sep="\t")
    
    cat("All done with dataset", grid$dataset[j], "alpha0=", grid[j, "alpha0"], "alpha1=", grid[j, "alpha1"], "\n")
}


  
    