#####################################################################################################################################################################
## Running the HGT + CARD on simulated data based on ALL OSCC samples with scRNAseq data simulated by SPARSim across our regular 9 combinations of hyperparameters ##
#####################################################################################################################################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

job.id <- job.id #+ 2000 # Adjust the job.id here if your HPC only goes to 1000 arrays like ours
cat(job.id)

library(pbmcapply)
library(SingleCellExperiment)
library(CARD)
library(data.table)

source("../HGTfunctions3.R")


###########################################
# Set up variables #
###########################################

ct.varname <- "cellype_fine"
ct.select <- c("myofibroblast", "cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast",
               "detox-iCAF", "macrophage", "endothelial", "dendritic ", "mast",
               "conventional CD4+ T-helper cells", "cytotoxic CD8+ T ", "Tregs", "cytotoxic CD8+ T exhausted")
sample.varname <- "orig.ident"


###############################
# Set up the grid of settings #
###############################

iseed <- c(56789, 123456 + 111111*0:98)
sample <- c(1,3:12)  # Already have output for sample 2 from main simulations
ntotal <- 10 
SC_replicate <- 1:10
ST_replicate <- 1:10
alpha0 <- c(-0.1, 0.0, 0.1)  #alpha_0 is actually xover.75minusx(invlogit(..)) of these values
alpha1 <- c(0.1, 0.3, 0.5) 

grid <- expand.grid(alpha0, alpha1, SC_replicate, ST_replicate, sample)
colnames(grid) <- c( "alpha0", "alpha1", "SC_replicate", "ST_replicate", "sample")
grid$iseed <- rep(rep(iseed, rep(9, length(iseed))))
grid$ntotal <- ntotal
grid$lib_factor <- 0.05
grid$Phi_factor <- 50 # Realistic sparsity settings


######################
# Run the HGT + CARD #
######################

# Need to only do 1 job at a time because some OSCC samples are much larger than others and take correspondingly longer
for (j in ((job.id - 1) * 1 + 1) : (job.id * 1)){

# Load the predefined layer label data and the location data
pattern_patho_label <- readRDS(paste0("./Reference/pattern_patho_label_", grid$sample[j], ".rds"))
location <- readRDS(paste0("./Reference/location_", grid$sample[j], ".rds"))

# Load the simulated scRNAseq data
split_sim_scRNAseq <- readRDS(paste0("./scRNAseq_Simulations/scRNAseq_sim_data/split_scRNAseq_sim_library_",
                                       grid$lib_factor[j], "_Phi_", grid$Phi_factor[j], "_rep_", grid$SC_replicate[j], ".rds"))
eset.sub.split1 <- split_sim_scRNAseq$Simulation ## was used to simulate the data
eset.sub.split2 <- split_sim_scRNAseq$Deconvolution ## use as the reference for downstream deconvolution analysis				       

# Load the simulated ST data
spatial.pseudo <- readRDS(paste0("./SPARSim_ST_Data/Sample_", grid$sample[j],
                                 "/SPARSim_ST_n", grid$ntotal[j], 
                                 "_cells_library_factor_", grid$lib_factor[j],
				                         "_Phi_factor_", grid$Phi_factor[j],
                                 "_scRNAseq_rep_", grid$SC_replicate[j],
                                 "_ST_replicate_", grid$ST_replicate[j], ".rds"))
                                 
spatial_location <- cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                     y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
rownames(spatial_location) <- colnames(spatial.pseudo$pseudo.data)

true_prop <- spatial.pseudo$true.p

file_name <- paste0("./SPARSim_Results/Sample_", grid$sample[j], 
                    "/HGT_SPARSim_pseudo_OSCC_n", grid$ntotal[j], 
                    "_cells_library_factor_", grid$lib_factor[j],
		                "_Phi_factor_", grid$Phi_factor[j],
                    "_scRNAseq_rep_", grid$SC_replicate[j],
                    "_ST_replicate_", grid$ST_replicate[j],
                    "_alpha0_", grid$alpha0[j], 
                    "_alpha1_", grid$alpha1[j], ".rds")
cat("Running", file_name, "\n")

sim <- SimOverCARD_OSCC(sc_count = eset.sub.split2@assays$simSC@counts,  
                        sc_meta = eset.sub.split2@meta.data,  
                        spatial_count = spatial.pseudo$pseudo.data,
                        spatial_location = spatial_location,  
                        ct.varname = "cellype_fine",
                        ct.select = c("myofibroblast", "cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast", "detox-iCAF", "macrophage", "endothelial", "dendritic ", "mast",
                                      "conventional CD4+ T-helper cells", "cytotoxic CD8+ T ", "Tregs", "cytotoxic CD8+ T exhausted"),
                        sample.varname = "orig.ident",
                        minCountGene = 100,
                        minCountSpot = 5,
                        truePropMat = spatial.pseudo$true.p,
                        scenario = 1,
                        alpha0 = xover.75minusx(invlogit( grid$alpha0[j] )),
                        alpha1 = grid$alpha1[j],
                        n = 100,
                        doWAIC = TRUE,
                        seed = grid$iseed[j])

saveRDS(sim, file_name)
gc()
cat("Done with row", j, "of grid!  This corresponds to sc_rep", grid$SC_replicate[j], "st_rep", grid$ST_replicate[j], "alpha0", grid$alpha0[j], "alpha1", grid$alpha1[j], "for lib_factor", grid$lib_factor[j], "and Phi_factor", grid$Phi_factor[j], ".\n")
}
