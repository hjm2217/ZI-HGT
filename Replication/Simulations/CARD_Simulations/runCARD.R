###############################################################################################################
## Running CARD on simulated MOB data for all noisiness levels and scenarios and different sequencing depths ##
###############################################################################################################

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

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

sc_eset.scenario5 = readRDS("./Reference/GSE109447.sceset.scenario5.RDS") ## single cell reference for scenario 5

## load the predefined layer label data and ST MOB location data
load("./Reference/pattern_gp_label.RData")  #loads in pattern_gp_label, which is just a vector named with each spatial location and the value refers to what layer of tissue was sampled there (1,2,3)
load("./Reference/sim_MOB_location.RData")  #loads in location dataframe, which contains the x and y coordinates of each spatial location for simulation


############
# Run CARD #
############

iseed <- c(56789, 123456 + 111111*0:98)
ntotal <- 9:12
mixnoise <- 0:3
replicate <- 1:100

grid <- expand.grid(replicate, mixnoise, ntotal)
colnames(grid) <- c("replicate", "mixnoise", "ntotal")
grid$iseed <- iseed

scenarios <- 1:5
results <- list()

results[[1]] <- data.frame()
results[[2]] <-	data.frame()
results[[3]] <-	data.frame()
results[[4]] <- data.frame()
results[[5]] <- data.frame()

#Use 100 arrays to each run 16 sets of data with all 5 scenarios

for (j in ((job.id - 1) * 16 + 1):(job.id*16)){

    #load in data
    spatial.pseudo <- readRDS(paste0("./Simulated_ST_Data/Mixnoise", grid$mixnoise[j],
                                     "/sim.pseudo.MOB.n", grid$ntotal[j],
                                     ".cellType6.Mixnoise", grid$mixnoise[j],
                                     ".rep", grid$replicate[j], ".rds"))
    eset.sub.split2 = split$eset.sub.split2
    ct.varname = "cellType"
    ct.select <- list()
    ct.select[[1]] <- c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
    ct.select[[2]] <- ct.select[[1]][ct.select[[1]] != "Neurons"]
    ct.select[[3]] <-  c(ct.select[[1]], "Blood")
    
    ct.select[[4]] <- ct.select[[1]]
    Combinations = combn(ct.select[[1]],2)
    #### for example, when we merge the cell types of Astrocytes Ependymal into one cell type
    icom = 5
    imerge = paste(Combinations[,icom],collapse = "_")
    #### create the new single cell RNAseq with the merged cell types
    eset.sub.split2.merged = eset.sub.split2
    ct.select.merged = c(ct.select[[4]][!(ct.select[[4]] %in% Combinations[,icom])],"Merged")
    pData(eset.sub.split2.merged)[,ct.varname] <- as.character(pData(eset.sub.split2.merged)[,ct.varname])
    pData(eset.sub.split2.merged)[,ct.varname][pData(eset.sub.split2.merged)[,ct.varname] %in% Combinations[,icom]] <- "Merged" 

    ct.select[[5]] <- c("Astrocytes","Neuron","Oligodendrocytes","Endothelial","Microglia","Ependymal")
    sc_eset.scenario5 = sc_eset.scenario5[,pData(sc_eset.scenario5)[,ct.varname] %in% ct.select[[5]]]
    
    sample.varname = "sampleID"
    spatial_location = cbind.data.frame(x=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",1)),
                                      y=as.numeric(sapply(strsplit(as.character(colnames(spatial.pseudo$pseudo.data)),split="x"),"[",2)))
    rownames(spatial_location) = colnames(spatial.pseudo$pseudo.data)
    truePropMat = spatial.pseudo$true.p

    tmplist <- list()
    tmplist[[1]] <- data.frame(Ncells = NA, SeqDepthTot = NA, SeqDepthPerLoc = NA, Replicate = NA, Seed = NA, RMSE = NA, RMSEcomp10 = NA, Scenario = NA, Mixnoise = NA)
    tmplist[[2]] <- data.frame(Ncells = NA, SeqDepthTot = NA, SeqDepthPerLoc = NA, Replicate = NA, Seed = NA, RMSE = NA, RMSEcomp10 = NA, Scenario = NA, Mixnoise = NA)
    tmplist[[3]] <- data.frame(Ncells = NA, SeqDepthTot = NA, SeqDepthPerLoc = NA, Replicate = NA, Seed = NA, RMSE = NA, RMSEcomp10 = NA, Scenario = NA, Mixnoise = NA)
    tmplist[[4]] <- data.frame(Ncells = NA, SeqDepthTot = NA, SeqDepthPerLoc = NA, Replicate = NA, Seed = NA, RMSE = NA, RMSEcomp10 = NA, Scenario = NA, Mixnoise = NA)
    tmplist[[5]] <- data.frame(Ncells = NA, SeqDepthTot = NA, SeqDepthPerLoc = NA, Replicate = NA, Seed = NA, RMSE = NA, RMSEcomp10 = NA, Scenario = NA, Mixnoise = NA)

    for (sc in scenarios){
    	 
	    tmplist[[sc]][1, c("Ncells", "Replicate", "Seed", "Mixnoise")] <- grid[j, c("ntotal", "replicate", "iseed", "mixnoise")]

	    if (sc %in% c(1,2,3)){
	      CARD_obj = createCARDObject(
	         sc_count = exprs(eset.sub.split2),  #18263 genes x 10208 cells dataset
    	     sc_meta = pData(eset.sub.split2),  #meta info for the 10k cells (cellname, celltype, sampleID)
    	     spatial_count = spatial.pseudo$pseudo.data,  #18263 genes x 260 locations
    	     spatial_location = spatial_location,  #spatial locations
   	       ct.varname = "cellType",
    	     ct.select = ct.select[[sc]],  #select all cells for appropriate scenario
    	     sample.varname = "sampleID",
    	     minCountGene = 100,
    	     minCountSpot = 5)
	    }

	    if (sc == 4){
	      CARD_obj = createCARDObject(
    	     sc_count = exprs(eset.sub.split2.merged),
    	     sc_meta = pData(eset.sub.split2.merged),
	         spatial_count = spatial.pseudo$pseudo.data,
    	     spatial_location = spatial_location,
    	     ct.varname = "cellType",
    	     ct.select = ct.select.merged,
    	     sample.varname = "sampleID",
    	     minCountGene = 100,
    	     minCountSpot = 5) 
	    }

	    if (sc == 5){
        CARD_obj = createCARDObject(
           sc_count = exprs(sc_eset.scenario5),
    	     sc_meta = pData(sc_eset.scenario5),
           spatial_count = spatial.pseudo$pseudo.data,
           spatial_location = spatial_location,
           ct.varname = "cellType",
           ct.select = ct.select[[sc]],
           sample.varname = "sampleID",
           minCountGene = 100,
           minCountSpot = 5)
       }
	 
  	   CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  	   saveRDS(CARD_obj, paste0("./Mixnoise", grid$mixnoise[j],
  	                            "/Sc", sc,
  	                            "/CARDresults/Outputsim.pseudo.MOB.n", grid$ntotal[j],
  	                            ".cellType6.Mixnoise",  grid$mixnoise[j],
  	                            ".rep", grid$replicate[j], ".rds"))
  	   saveRDS(sim, paste0("./Results_CARD/Mixnoise", grid$mixnoise[j],
  	                       "/Sc", sc, 
  	                       "/HGTsim.pseudo.MOB.n10.cellType6.Mixnoise", grid$mixnoise[j],
  	                       ".rep", grid$replicate[j],
  	                       ".alpha0_", grid$alpha0[j], 
  	                       ".alpha1_", grid$alpha1[j], ".rds"))
  	   
	     tmplist[[sc]]$SeqDepthTot[1] <- sum(spatial.pseudo$pseudo.data)
  	   tmplist[[sc]]$SeqDepthPerLoc[1] <- sum(spatial.pseudo$pseudo.data) / ncol(spatial.pseudo$pseudo.data)
  	   tmplist[[sc]]$RMSE[1] <- rmse(CARD_obj@Proportion_CARD, truePropMat, scenario = sc)
	     tmplist[[sc]]$Scenario[1] <- sc
	     tmplist[[sc]]$Mixnoise[1] <- grid$mixnoise[j]
	     tmplist[[sc]]$Ncells[1] <- grid$ntotal[j]
	     tmplist[[sc]]$Replicate[1] <- grid$replicate[j]
	     tmplist[[sc]]$Seed[1] <- grid$iseed[j]

	     results[[sc]] <- rbind(results[[sc]], tmplist[[sc]])
	     gc()
	
    }

    cat("Done with dataset", j, ".\n")

}

saveRDS(results, paste0("Results_CARD/CARD.sim.MOB.res.", job.id, ".rds"))

  
    