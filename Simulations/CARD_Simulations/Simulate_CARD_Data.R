#########################################################################
## Simulating CARD data with different Sequencing Depths and Mixnoises ##
#########################################################################

# This simulated data is based on the Mouse Olfactory Bulb data as shown in the CARD paper
# https://www.nature.com/articles/s41587-022-01273-7
# We make two small changes to their simulations.  1) We analyze 100 simulated 
# datasets rather than 5.  2) We adjust the overall sequencing depth  by increasing/
# decreasing the number of cells present at each location.  We then compare CARD with 
# different sequencing depths to the ZI-HGT + CARD.

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
job.id <- as.numeric(slurm_arrayid)

library(pbmcapply)
library(SingleCellExperiment)
library(CARD)

load("./Reference/split.scRNAseq.forSimu.RData")
eset.sub.split1 = split$eset.sub.split1 ## simulate the data
eset.sub.split2 = split$eset.sub.split2 ## use as the reference for downstream deconvolution analysis
ct.varname = "cellType"
ct.select = c("Astrocytes","Neurons","Oligos","Vascular","Immune","Ependymal")
sample.varname = "sampleID"

## load the predefined layer label data and ST MOB location data
load("./Reference/pattern_gp_label.RData")  #loads in pattern_gp_label, which is just a vector named with each spatial location and the value refers to what layer of tissue was sampled there (1,2,3)
load("./Reference/sim_MOB_location.RData")  #loads in location dataframe, which contains the x and y coordinates of each spatial location for simulation


#############################
# CARD Simulation Functions #
#############################

#### generate random numbers from Dirichelet distribution
generateMultiN <- function(pattern_gp_label,ipt,ntotal,mix,ct.select){
  library(MCMCpack)
  message(paste0("Generating cell type proportions for pattern",ipt))  #the different ipt's are patterns (1,2, or 3) of cell-types for layers 1, 2, 3 of tissue.  See lines 55-64
  nobs = sum(pattern_gp_label == ipt)  #pattern_gp (group)_label assigns group label to each observation, these are all 1 for each location
  sample = matrix(0,nrow = nobs,ncol = length(ct.select))  #set up the sample matrix for whichever pattern (layer, ipt) we're looking at
  colnames(sample) = ct.select  #name the columns
  sampleNames = names(pattern_gp_label)[pattern_gp_label == ipt]  #gets the names (locations) for which we've observed the layer (ipt)
  prop = matrix(0,nrow = nobs,ncol = length(ct.select))  #set up the cell-type proportion matrix
  colnames(prop) = ct.select  #name the columns
  ## sample total number of cell types in each layer: main cell type + colocalized cell types
  numCT = sample(1:length(ct.select),1) 
  if(ipt == 1){
    main = "Neurons" ### defined one dominant cell type in layer 1
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 2){ 
    main = "Astrocytes" ### defined one dominant cell type in layer 2
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if(ipt == 3){
    main = "Oligos" ### defined one dominant cell type in layer 3
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }
  propSample = rdirichlet(nrow(prop),concen)  #generates samples from a Dirichlet distribution with parameters = 1, note that the Dirichlet dist is a multivariate generalization of the beta dist, and the support is X = (X1, X2,..., Xn), where  sum(Xi) = 1 and Xi in [0,1] for all i
  ct.select.sub.sample = sample(ct.select[ct.select != main],numCT-1)  #all cell types included other than the dominant cell type, in a random order
  ct.select.sub = c(main,ct.select.sub.sample)  #put the main cell type first
  Fix_Dirichlet = sample(1:nrow(sample),round(nobs * mix[1])) #set the locations that aren't noisy - i.e. have the dominant cell type receiving the highest proportion
  mix_Dirichlet = c(1:nrow(sample))[!(c(1:nrow(sample)) %in% Fix_Dirichlet)] #set the locations that are noisy, the dominant cell type doesn't necessarily get the highest proportion, set these according to mix, which comes from pn:proportion of noisy locations
  if(length(mix_Dirichlet) > 0){
    propSample[mix_Dirichlet,] = rdirichlet(length(mix_Dirichlet),rep(1,numCT))  #randomly set proportions for noisy locations
  }
  print(ct.select.sub)
  if(length(Fix_Dirichlet) > 0){  #set proportions for non-noisy locations
    propSample[Fix_Dirichlet,] = t(sapply(Fix_Dirichlet,function(i){
      propSample[i,][order(propSample[i,],decreasing = T)] ### for non-noisy positions, order the proportion to assign the largest proportion to the dominant cell type
    }))
  }
  colnames(propSample) = ct.select.sub  #set the cell-type proportion matrix column names (note that this is for whatever layer we're working on)
  prop[,ct.select.sub] = propSample  #insert the cell-type proportions into the matrix that isn't organized by dominant cell type first (so it'll match the format for the other layers)
  sample = round(ntotal * prop,digits = 0)  #number of cells of each type at each location, multiply total cells at each location (10) by cell-type proportion generated from simulation
  ##### to avoid no numbers sample
  index = which(rowSums(sample) == 0)  #if for some reason we don't have any cells in a location, fix that by using a multinomial distribution
  if(length(index) > 0){
    sample[index,] = t(sapply(index,function(i){rmultinom(1, ntotal, prob = prop[i,])}))
  }
  ####
  return(list(sample = sample))
}

#made it through the first function, move onto the second one later

generateSpatial_norep_fixedProp <- function(seed, eset.sub.split1, ct.varname, sample.varname,ct.select,sample.withRep = F,pattern_gp_label, ntotal,mix1,mix2,mix3){
  #### using split 1 to sample the single cell RNAseq data
  phenoData <- eset.sub.split1@phenoData@data  #10k by 3 dataset with cell name, cell type, and sample ID
  # number of cell type of interest
  k <- length(unique(ct.select))
  message(paste('Using',k,'cell types to generate pseudo spatial dataset'))
  # select donors for each pseudo bulk sample.varname
  Sample_random = matrix(data = 0,ncol = length(ct.select), nrow = length(pattern_gp_label))
  rownames(Sample_random) = names(pattern_gp_label)
  colnames(Sample_random) = ct.select
  ##### Total number of cells in the subset
  ct.id <- droplevels(as.factor(eset.sub.split1@phenoData@data[,ct.varname]))
  library(Hmisc)
  ### random generate the number of cells for each spatial location in each layer
  pattern1 = generateMultiN(pattern_gp_label,1,ntotal,mix1,ct.select)  #the 1,2,3 are the ipt's, which are just the different tissue layers
  pattern2 = generateMultiN(pattern_gp_label,2,ntotal,mix2,ct.select)  #these are matrices with locations on rows and cell types on columns, with cell counts as the entries
  pattern3 = generateMultiN(pattern_gp_label,3,ntotal,mix3,ct.select)
  
  Sample_random[pattern_gp_label == 1,] = pattern1$sample  #stick the cell counts into the sample matrix
  Sample_random[pattern_gp_label == 2,] = pattern2$sample
  Sample_random[pattern_gp_label == 3,] = pattern3$sample
  
  message(paste0("Generating pseudo spatial dataset for ",length(unique(pattern_gp_label))," patterns"))
  ##### calculate the number of cells we need to sample in each cell type
  set.seed(seed)
  temp.exprs <- exprs(eset.sub.split1)  #expression from first split of scRNA-seq data, genes on rows cells on cols
  temp.nct <- Sample_random  #temporary number of each type of cell
  true.p = sweep(temp.nct,1,rowSums(temp.nct),"/")  #divide the cell counts at each location by the number of cells at that location to get the true cell type proportions
  ##### Randomness in pbmclapply
  temp.pseudo = pbmclapply(1:nrow(temp.nct),function(isample){  #This creates the spatially resolved gene expression matrix
    temp.sample = temp.nct[isample,]  #take the row isample from the cell counts matrix temp.nct
    temp.sample.pseudo = sapply(ct.select,function(ict){
      temp.vec <- temp.exprs[,ct.id %in% ict]  #get the expression for cells for which their types are in a specific layer
      if(temp.nct[isample,ict] > ncol(temp.vec)){  #if the number of  cells in row isample column names ict (cell types that live in layer ict) is larger than the columns of the temporary expression vector, we'll need to sample with replacement
        sample.withRep = T
      }
      temp.id <- sample(1:ncol(temp.vec), temp.nct[isample,ict], replace = sample.withRep)  #sample from 1:n, where n is the number of cells expressed in the current layer, sample size is temp.nct[isample,ict] (this is choosing which cells from say astrocytes is being expressed)
      if(length(temp.id) > 1){
        rowSums(temp.vec[,temp.id])
      }else if(length(temp.id) == 1){
        temp.vec[,temp.id]
      }else if(length(temp.id) == 0){
        rep(0,nrow(temp.vec))
      }
    })
    rowSums(temp.sample.pseudo)
  },mc.cores = 2,mc.set.seed = F) ## Set the number of cores that can be used in the system
  temp.pseudo = do.call("cbind",temp.pseudo)
  colnames(temp.pseudo) = rownames(Sample_random)
  rownames(temp.pseudo) = rownames(temp.exprs)
  return(list(pseudo.data = temp.pseudo, true.p = true.p, ntotal = ntotal,Sample_random = Sample_random))
}  



#####################
# Simulate the data #
#####################

#We want 100 replicates, 4 levels of noisiness (mixnoise = 0,1,2,3), and 4 numbers of cells (ntotal = 9, 10 (baseline), 11, 12)
iseed <- c(56789, 123456 + 111111*0:98)
ntotal <- 9:12
mixnoise <- 0:3
replicate <- 1:100

grid <- expand.grid(replicate, mixnoise, ntotal)
colnames(grid) <- c("replicate", "mixnoise", "ntotal")
grid$iseed <- iseed

#Use 100 arrays to each create 16 sets of data

for (j in ((job.id - 1) * 16 + 1):(job.id*16)){

    imix <- grid$mixnoise[j]
    mix1 = mix2 = mix3 = c(1 - (0.2 * imix),0.2*imix)

    set.seed(grid$iseed[j])

    spatial.pseudo = generateSpatial_norep_fixedProp(
    seed = grid$iseed[j],
    eset.sub.split1 = eset.sub.split1,
    ct.varname = ct.varname,
    sample.varname = sample.varname,
    ct.select = ct.select,
    sample.withRep = F,
    pattern_gp_label = pattern_gp_label,
    ntotal = grid$ntotal[j],
    mix1 = mix1,
    mix2 = mix2,
    mix3 = mix3)

    saveRDS(spatial.pseudo, paste0("./Simulated_ST_Data/Mixnoise", grid$mixnoise[j],
                                   "/sim.pseudo.MOB.n", grid$ntotal[j],
                                   ".cellType6.Mixnoise", grid$mixnoise[j],
                                   ".rep", grid$replicate[j], ".rds"))

    gc()
    cat("Simulation", j, "done.\n")
}
