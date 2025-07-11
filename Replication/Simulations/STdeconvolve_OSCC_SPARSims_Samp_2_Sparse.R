#################################################################
## Run STdeconvolve on the simulated OSCC data as a comparison ##
#################################################################

# I'm adapting the code from the STdeconvolve vignette (https://bioconductor.org/packages/devel/bioc/vignettes/STdeconvolve/inst/doc/vignette.html#:~:text=STdeconvolve%20is%20an%20unsupervised%20machine,reliance%20on%20external%20single%2Dcell)
library(STdeconvolve)

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

# Set up the grid of settings
iseed <- c(56789, 123456 + 111111*0:98)
sample <- 2
ntotal <- 10
SC_replicate <- 1:10
ST_replicate <- 1:10

grid <- expand.grid(SC_replicate, ST_replicate, sample)
colnames(grid) <- c( "SC_replicate", "ST_replicate", "sample")
grid$iseed <- iseed
grid$ntotal <- ntotal
grid$lib_factor <- 0.05
grid$Phi_factor <- 50

j <- job.id

# Load in the pathologist labels and location data
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

# Set things up like STdeconvolve
pos <- location
rownames(pos) <- paste0(pos$x, "x", pos$y)
cd <- spatial.pseudo$pseudo.data
annot <- pattern_patho_label

## remove pixels with too few genes
counts <- cleanCounts(counts = cd,
                      min.lib.size = 100,
                      min.reads = 1,
                      min.detected = 1,
                      verbose = TRUE)

## feature select for genes
corpus <- restrictCorpus(counts,
                         removeAbove = 1.0,
                         removeBelow = 0.05,
                         alpha = 0.05,
                         plot = F,
                         verbose = TRUE)

## Note: the input corpus needs to be an integer count matrix of pixels x genes
# Note from me, Ks is the different numbers of cell types (as in we don't know it), so set Ks = 14
ldas <- fitLDA(t(as.matrix(corpus)), Ks = 14,
               perc.rare.thresh = 0.05,
               plot=FALSE,
               verbose=TRUE)

## select model with minimum perplexity
optLDA <- optimalModel(models = ldas, opt = "min")

## Extract pixel cell-type proportions (theta) and cell-type gene expression
## profiles (beta) for the given dataset.
## We can also remove cell-types from pixels that contribute less than 0% of the
## pixel proportion and scale the deconvolved transcriptional profiles by 1000 
results <- getBetaTheta(optLDA,
                        perc.filt = -0.01,
                        betaScale = 1000)

deconProp <- results$theta
colnames(deconProp)
deconGexp <- results$beta

# We need to identify which of the STdeconvolve cell types correspond to the cell types in the scRNAseq reference
scAnnot <- eset.sub.split2@meta.data$cellype_fine
names(scAnnot) <- rownames(eset.sub.split2@meta.data)
scProxyTheta <- model.matrix(~ 0 + scAnnot)

rownames(scProxyTheta) <- names(scAnnot)
# fix names
colnames(scProxyTheta) <- unlist(lapply(colnames(scProxyTheta), function(x) {
  unlist(strsplit(x, "scAnnot"))[2]
}))

sc_counts <- eset.sub.split2@assays$simSC$counts
sc_counts <- sc_counts[rownames(sc_counts) %in% colnames(deconGexp),]
scProxyGexp <- sc_counts %*% scProxyTheta

corMtx_beta <- getCorrMtx(m1 = as.matrix(deconGexp), # the deconvolved cell-type `beta` (celltypes x genes)
                          m2 = t(as.matrix(scProxyGexp)), # the reference `beta` (celltypes x genes) - comes from scRNAseq data for us
                          type = "b") # "b" = comparing beta matrices, "t" for thetas
rownames(corMtx_beta) <- paste0("decon_", seq(nrow(corMtx_beta)))

# Now we need to use the correlation matrix to assign cell types
assignment_df <- data.frame("Cell_type"= rep(NA, 14), "Decon" = 1:14)
corMtx_beta_update <- corMtx_beta
for (i in c(2,1,3:13)){
#for (i in c(2,1,3)){
  bestIdx <- which(corMtx_beta_update[,i] == max(corMtx_beta_update[,i]))[1]
  matchIdx <- which(rownames(corMtx_beta) == rownames(corMtx_beta_update)[bestIdx])
  assignment_df$Cell_type[matchIdx] <- colnames(corMtx_beta)[i]
  corMtx_beta_update <- corMtx_beta_update[-(bestIdx),]
  
  # Last one
  if (i ==13){
    assignment_df$Cell_type[is.na(assignment_df$Cell_type)] <- colnames(corMtx_beta)[i+1]
  }
}

colnames(deconProp) <- assignment_df$Cell_type
ct.select <- c("myofibroblast", "cancer cell", "B cell", "ecm-myCAF", "Intermediate fibroblast",
               "detox-iCAF", "macrophage", "endothelial", "dendritic ", "mast",
               "conventional CD4+ T-helper cells", "cytotoxic CD8+ T ", "Tregs", "cytotoxic CD8+ T exhausted")
deconProp_save <- deconProp[, ct.select]

out_dir = "SPARSim_Results/Sample_2/STdeconvolve"
dir.create(out_dir,recursive = TRUE, showWarnings = FALSE)
out_matrix_fp <- file.path(out_dir,sprintf(paste0("STdeconvolve_SPARSim_pseudo_OSCC_n", grid$ntotal[j],
                                               "_cells_library_factor_", grid$lib_factor[j],
                                               "_Phi_factor_", grid$Phi_factor[j],
                                               "_scRNAseq_rep_", grid$SC_replicate[j],
                                               "_ST_replicate_", grid$ST_replicate[j], ".rds")))

write.csv(deconProp_save,out_matrix_fp)