##############################################################
## Run SPOTlight on the simulated OSCC data as a comparison ##
##############################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

library(SPOTlight)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)

# I'm adapting the code from the SPOTlight documentation (https://bioconductor.org/packages/devel/bioc/vignettes/SPOTlight/inst/doc/SPOTlight_kidney.html#deconvolution)
# and the cell-type deconvolution comparison paper (https://github.com/leihouyeung/STdeconv_benchmark/blob/main/methods/SPOTlight/SPOTlight.R)


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

# Just need to set the spatial data up in genes on rows and locations on columns format
spe <- spatial.pseudo$pseudo.data

# Prep the single-cell data and transform the data (following SPOTlight)
sc_seurat <- CreateSeuratObject(counts = eset.sub.split2@assays$simSC@counts)
sce <- as.SingleCellExperiment(sc_seurat)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 3000)
colLabels(sce) <- eset.sub.split2@meta.data[, "cellype_fine"]
mgs <- scoreMarkers(sce)
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, average fold change > 2
  x <- x[x$mean.logFC.detected > 2, ] 
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.logFC.detected, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

# Downsample to at most 100 cells per cell type (following SPOTlight recommendation)
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$label)
# downsample to at most 100 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real life analysis
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

# Run SPOTlight, do downsampling as recommended
# Went from about 3k cells to 1.2k
# Make sure it runs with like 5 cells each then do 100 on the actual analysis
res <- SPOTlight(
  x = sce,
  y = spe,
  groups = sce$label,
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.logFC.detected",
  group_id = "cluster",
  gene_id = "gene")

# Extract deconvolution matrix 
decon_results <- res$mat

out_dir = "SPARSim_Results/Sample_2/SPOTlight"
#out_dir="./seqFISH_10000_Result"
dir.create(out_dir,recursive = TRUE, showWarnings = FALSE)
out_matrix_fp=file.path(out_dir,sprintf(paste0("SPOTlight_SPARSim_pseudo_OSCC_n", grid$ntotal[j],
                                                    "_cells_library_factor_", grid$lib_factor[j],
                                                    "_Phi_factor_", grid$Phi_factor[j],
                                                    "_scRNAseq_rep_", grid$SC_replicate[j],
                                                    "_ST_replicate_", grid$ST_replicate[j], ".rds")))

write.csv(decon_results,out_matrix_fp)
