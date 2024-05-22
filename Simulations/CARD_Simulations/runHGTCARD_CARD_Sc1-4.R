###########################################################################################################
## Running CARD with HGT on simulated MOB data (ncells = 10) for scenarios 1-4 with all noisiness levels ##
###########################################################################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)
job.id <- job.id # + 1000 # Adjust the job.id here if your HPC only goes to 1000 arrays like ours
cat(job.id, "\n")

library(pbmcapply)
library(SingleCellExperiment)
library(CARD)

source("../HGTfunctions3.R")  


#############
# Load Data #
#############

load("./Reference/split.scRNAseq.forSimu.RData")
eset.sub.split1 = split$eset.sub.split1 ## simulate the data
eset.sub.split2 = split$eset.sub.split2 ## use as the reference for downstream deconvolution analysis
ct.varname = "cellType"
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
sample.varname = "sampleID"

#sc_eset.scenario5 = readRDS("./GSE109447.sceset.scenario5.RDS") ## single cell reference for scenario 5

## load the predefined layer label data and ST MOB location data
load("./Reference/pattern_gp_label.RData")  #loads in pattern_gp_label, which is just a vector named with each spatial location and the value refers to what layer of tissue was sampled there (1,2,3)
load("./Reference/sim_MOB_location.RData")  #loads in location dataframe, which contains the x and y coordinates of each spatial location for simulation

#####################
# Run CARD with HGT #
#####################

iseed <- c(56789, 123456 + 111111*0:98)
mixnoise <- 0:3
replicate <- 1:100
alpha0 <- c(-0.1, 0.0, 0.1)  #alpha_0 is actually xover.75minusx(invlogit(..)) of these four values
alpha1 <- c(0.1, 0.3, 0.5)


grid <- expand.grid(alpha0, alpha1, replicate, mixnoise)
colnames(grid) <- c("alpha0", "alpha1", "replicate", "mixnoise")
grid$iseed <- rep(iseed, rep(9, length(iseed)))

scenarios <- 1:4 # Only do first 4 scenarios for this job, 5 needs its own


for (j in ((job.id - 1) * 2  + 1):(job.id*2)){

    #load in data
    spatial.pseudo <- readRDS(paste0("./Simulated_ST_Data/Mixnoise", grid$mixnoise[j],
                                     "/sim.pseudo.MOB.n10.cellType6.Mixnoise", grid$mixnoise[j],
                                     ".rep", grid$replicate[j], ".rds"))
    eset.sub.split2 = split$eset.sub.split2
    ct.varname = "cellType"
    ct.select <- list()
    ct.select[[1]] <- c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
    ct.select[[2]] <- ct.select[[1]][ct.select[[1]] != "Neurons"]
    ct.select[[3]] <-  c(ct.select[[1]], "Blood")

    ct.select[[4]] <- ct.select[[1]]
    Combinations = combn(ct.select[[1]],2)
    #### for example, when we merge the cell types Astrocytes & Ependymal into one cell type
    icom = 5
    imerge = paste(Combinations[,icom],collapse = "_")
    #### create the new single cell RNAseq with the merged cell types
    eset.sub.split2.merged = eset.sub.split2
    ct.select.merged = c(ct.select[[4]][!(ct.select[[4]] %in% Combinations[,icom])],"Merged")
    pData(eset.sub.split2.merged)[,ct.varname] <- as.character(pData(eset.sub.split2.merged)[,ct.varname])
    pData(eset.sub.split2.merged)[,ct.varname][pData(eset.sub.split2.merged)[,ct.varname] %in% Combinations[,icom]] <- "Merged"

    #ct.select[[5]] <- c("Astrocytes","Neuron","Oligodendrocytes","Endothelial","Microglia","Ependymal")
    #sc_eset.scenario5 = sc_eset.scenario5[,pData(sc_eset.scenario5)[,ct.varname] %in% ct.select[[5]]]

    sample.varname = "sampleID"
    spatial_location = cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                      y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
    rownames(spatial_location) = colnames(spatial.pseudo$pseudo.data)
    truePropMat = spatial.pseudo$true.p


    for (sc in scenarios){
    	 filename <- paste0("./Results_CARD/Mixnoise", grid$mixnoise[j],
    	                    "/Sc", sc, 
    	                    "/HGTsim.pseudo.MOB.n10.cellType6.Mixnoise", grid$mixnoise[j],
    	                    ".rep", grid$replicate[j],
    	                    ".alpha0_", grid$alpha0[j], 
    	                    ".alpha1_", grid$alpha1[j], ".rds")
	     if (file.exists(filename)){
	       cat(paste0("Scenario ", sc, "already done. Next.\n"))
	       next
	     }
    	 if (sc %in% c(1,2,3)){
	       sim <- SimOverCARD(sc_count = exprs(eset.sub.split2),  
                     sc_meta = pData(eset.sub.split2),  #meta info for the 10k cells (cellname, celltype, sampleID)
                     spatial_count = spatial.pseudo$pseudo.data,
                     spatial_location = spatial_location,  #spatial locations
                     ct.varname = "cellType",
                     ct.select = ct.select[[sc]],  
                     sample.varname = "sampleID",
                     minCountGene = 100,
                     minCountSpot = 5,
                     truePropMat = spatial.pseudo$true.p,
                     scenario = sc,
                     alpha0 = xover.75minusx(invlogit( grid[j, "alpha0"] )),
                     alpha1 = grid[j, "alpha1"],  
                     n = 100,
		                 doWAIC = TRUE,
                     seed = grid[j, "iseed"])
	     }

	     if (sc == 4){
         sim <- SimOverCARD(sc_count = exprs(eset.sub.split2.merged),
               	     sc_meta = pData(eset.sub.split2.merged),
             	       spatial_count = spatial.pseudo$pseudo.data,
             	       spatial_location = spatial_location,
             	       ct.varname = "cellType",
             	       ct.select = ct.select.merged,
             	       sample.varname = "sampleID",
             	       minCountGene = 100,
	                   minCountSpot = 5,
                     truePropMat = spatial.pseudo$true.p,
                     scenario = sc,
		                 alpha0 = xover.75minusx(invlogit( grid[j, "alpha0"] )),
                     alpha1 = grid[j, "alpha1"],
                     n = 100,
		                 doWAIC = TRUE,
                     seed = grid[j, "iseed"])
       }

	     if (sc == 5){
         sim <- SimOverCARD(sc_count = exprs(sc_eset.scenario5),
                     sc_meta = pData(sc_eset.scenario5),
             	       spatial_count = spatial.pseudo$pseudo.data,
             	       spatial_location = spatial_location,
             	       ct.varname = "cellType",
             	       ct.select = ct.select[[sc]],
             	       sample.varname = "sampleID",
             	       minCountGene = 100,
             	       minCountSpot = 5,
	                   truePropMat = spatial.pseudo$true.p,
                     scenario = sc,
                     alpha0 = xover.75minusx(invlogit( grid[j, "alpha0"] )),
                     alpha1 = grid[j, "alpha1"],
                     n = 100,
		                 doWAIC = TRUE,
                     seed = grid[j, "iseed"])
       }
		

	     saveRDS(sim, paste0("./Results_CARD/Mixnoise", grid$mixnoise[j],
	                         "/Sc", sc, 
	                         "/HGTsim.pseudo.MOB.n10.cellType6.Mixnoise", grid$mixnoise[j],
	                         ".rep", grid$replicate[j],
	                         ".alpha0_", grid$alpha0[j], 
	                         ".alpha1_", grid$alpha1[j], ".rds"))
     	gc()
    	cat("#############################\nDONE WITH SCENARIO", sc, "/ 5! \n#############################")
    }
    cat("#############################\nDONE WITH DATASET", j, "\n#############################")

}


  
    