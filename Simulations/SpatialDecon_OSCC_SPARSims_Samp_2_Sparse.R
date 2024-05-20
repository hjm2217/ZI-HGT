#################################################################
## Run SpatialDecon on the simulated OSCC data as a comparison ##
#################################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

library(SpatialDecon)
library(SeuratObject)
library(Matrix)

# I'm adapting the code from the SpatialDecon documentation (https://bioconductor.org/packages/release/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_ST.html)
# and the cell-type deconvolution comparison paper (https://github.com/leihouyeung/STdeconv_benchmark/blob/main/methods/SpatialDecon/SpatialDecon.R)

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

data <- list(
  st_counts = spatial.pseudo$pseudo.data,
  st_locations = spatial_location, 
  sc_counts = eset.sub.split2@assays$simSC@counts,
  sc_labels = eset.sub.split2@active.ident
)



preprocess=function(data){
  st_counts_norm = sweep(data$st_counts, 2, colSums(data$st_counts), "/") * mean(colSums(data$st_counts))
  st_object=CreateSeuratObject(counts=st_counts_norm,assay="Spatial")
  stopifnot(setequal(colnames(st_object),rownames(data$st_location)))
  st_object=AddMetaData(st_object,data$st_locations[colnames(st_object),1],col.name="x")
  st_object=AddMetaData(st_object,data$st_locations[colnames(st_object),2],col.name="y")
  
  stopifnot(all(colnames(data$sc_counts)==names(data$sc_labels)))
  
  sc_counts_matrix=as.matrix(data$sc_counts)
  sc_counts_matrix=Matrix::Matrix((sc_counts_matrix),sparse=TRUE)
  sc_labels_df=data.frame(cell_barcodes=names(data$sc_labels),sc_labels=data$sc_labels)
  sc_matrix <- create_profile_matrix(
    mtx = sc_counts_matrix,            # cell x gene count matrix
    cellAnnots = sc_labels_df,  # cell annotations with cell type and cell name as columns 
    cellTypeCol = "sc_labels",  # column containing cell type
    cellNameCol = "cell_barcodes",           # column containing cell ID/name
    matrixName = "custom_cell_type_matrix", # name of final profile matrix
    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
    normalize = TRUE,                
    minCellNum = 1,
    minGenes = 1
  ) 
  
  return(
    list(
      st_object=st_object,
      sc_matrix=sc_matrix
    )
  )
}

#data = load_seqFISH(10000)
out_dir = "SPARSim_Results/Sample_2/SpatialDecon"
#out_dir="./seqFISH_10000_Result"
dir.create(out_dir,recursive = TRUE, showWarnings = FALSE)
out_matrix_norm_fp=file.path(out_dir,sprintf(paste0("SpatialDecon_SPARSim_pseudo_OSCC_n", grid$ntotal[j],
                                                    "_cells_library_factor_", grid$lib_factor[j],
                                                    "_Phi_factor_", grid$Phi_factor[j],
                                                    "_scRNAseq_rep_", grid$SC_replicate[j],
                                                    "_ST_replicate_", grid$ST_replicate[j], ".rds")))

processed_data=preprocess(data)
start_time <- Sys.time()
res = runspatialdecon(object = processed_data$st_object,
                      bg = 0.01,
                      X = processed_data$sc_matrix,
                      align_genes = TRUE)
end_time <- Sys.time()

weights=t(res$beta)
norm_weights=sweep(weights, 1, rowSums(weights), "/")
norm_weights <- norm_weights[, colnames(true_prop)]
write.csv(as.matrix(norm_weights),out_matrix_norm_fp)