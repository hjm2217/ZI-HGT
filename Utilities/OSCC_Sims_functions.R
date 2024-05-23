library(pbmcapply)
library(SingleCellExperiment)
library(CARD)
library(Seurat)
library(dplyr)

#############################
# CARD Simulation Functions #
#############################

library(pbmcapply)
library(SingleCellExperiment)
library(CARD)
library(Seurat)
library(dplyr)

###############################################
# CARD Simulation Functions for the OSCC data #
###############################################

#### generate random numbers from Dirichelet distribution
generateMultiN_OSCC <- function(pattern_gp_label,ipt,ntotal,mix,ct.select){
  library(MCMCpack)
  message(paste0("Generating cell type proportions for pattern",ipt))  #the different ipt's are patterns (1,2, or 3) of cell-types for layers 1, 2, 3 of tissue.  See lines 55-64
  nobs = sum(pattern_gp_label == ipt)  #pattern_gp (group)_label assigns group label to each observation, these are all 1 for each location
  if (nobs == 0){
    message(paste0("Pattern", ipt, " is not in this sample. Skipping..."))
    return(NA)
  }

  sample = matrix(0,nrow = nobs,ncol = length(ct.select))  #set up the sample matrix for whichever pattern (layer, ipt) we're looking at
  colnames(sample) = ct.select  #name the columns
  sampleNames = names(pattern_gp_label)[pattern_gp_label == ipt]  #gets the names (locations) for which we've observed the layer (ipt)
  prop = matrix(0,nrow = nobs,ncol = length(ct.select))  #set up the cell-type proportion matrix
  colnames(prop) = ct.select  #name the columns
  ## sample total number of cell types in each layer: main cell type + at least one colocalized cell type
  numCT = sample(2:length(ct.select),1)
  if(ipt == 1){
    main = "endothelial" ### defined one dominant cell type in layer 1
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 2){
    main = "Intermediate fibroblast" ### defined one dominant cell type in layer 2
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 3){
    main = "ecm-myCAF" ### defined one dominant cell type in layer 3
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 4){
    main = "cytotoxic CD8+ T " ### defined one dominant cell type in layer 4
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 5){
    main = "myofibroblast" ### defined one dominant cell type in layer 5
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 6){
    main = "macrophage" ### defined one dominant cell type in layer 6
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 7){
    main = "cancer cell" ### defined one dominant cell type in layer 7
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

generateSpatial_norep_fixedProp_OSCC <- function(seed, eset.sub.split1, ct.varname, sample.varname, ct.select, RNA.counts, sample.withRep = F, pattern_gp_label, ntotal, mix = 1:0){
  #### using split 1 to sample the single cell RNAseq data
  # # phenoData <- eset.sub.split1@phenoData@data  #10k by 3 dataset with cell name, cell type, and sample ID
  phenoData <- data.frame(cellname = rownames(eset.sub.split1@meta.data),
                          cellype_fine = eset.sub.split1@meta.data$cellype_fine,
                          orig.ident = eset.sub.split1@meta.data$orig.ident)
  rownames(phenoData) <- phenoData$cellname

  # number of cell type of interest
  k <- length(unique(ct.select))
  message(paste('Using',k,'cell types to generate pseudo spatial dataset'))
  # select donors for each pseudo bulk sample.varname
  Sample_random = matrix(data = 0,ncol = length(ct.select), nrow = length(pattern_gp_label))
  rownames(Sample_random) = names(pattern_gp_label)
  colnames(Sample_random) = ct.select
  ##### Total number of cells in the subset
  # # ct.id <- droplevels(as.factor(eset.sub.split1@phenoData@data[,ct.varname]))
  ct.id <- droplevels(as.factor(phenoData[,ct.varname]))

  library(Hmisc)
  ### random generate the number of cells for each spatial location in each layer - even the layers not present in the data
  pattern <- list()
  for (i in 1:7){
    pattern[[i]] <- generateMultiN_OSCC(pattern_gp_label,i,ntotal,mix,ct.select)
    if (!is.na(pattern[[i]])){
      Sample_random[pattern_gp_label == i,] = pattern[[i]]$sample  #stick the cell counts into the sample matrix
    }
  }
  # # pattern1 = generateMultiN_OSCC(pattern_gp_label,1,ntotal,mix,ct.select)  #the 1,2,3 are the ipt's, which are just the different tissue layers
  # # pattern2 = generateMultiN_OSCC(pattern_gp_label,2,ntotal,mix,ct.select)  #these are matrices with locations on rows and cell types on columns, with cell counts as the entries
  # # pattern3 = generateMultiN_OSCC(pattern_gp_label,3,ntotal,mix,ct.select)
  # # pattern4 = generateMultiN_OSCC(pattern_gp_label,4,ntotal,mix,ct.select)
  # # pattern5 = generateMultiN_OSCC(pattern_gp_label,5,ntotal,mix,ct.select)
  # # pattern6 = generateMultiN_OSCC(pattern_gp_label,6,ntotal,mix,ct.select)
  # # pattern7 = generateMultiN_OSCC(pattern_gp_label,7,ntotal,mix,ct.select)

  # # Sample_random[pattern_gp_label == 1,] = pattern1$sample  #stick the cell counts into the sample matrix
  # # Sample_random[pattern_gp_label == 2,] = pattern2$sample
  # # Sample_random[pattern_gp_label == 3,] = pattern3$sample
  # # Sample_random[pattern_gp_label == 4,] = pattern4$sample
  # # Sample_random[pattern_gp_label == 5,] = pattern5$sample
  # # Sample_random[pattern_gp_label == 6,] = pattern6$sample
  # # Sample_random[pattern_gp_label == 7,] = pattern7$sample

  message(paste0("Generating pseudo spatial dataset for ",length(unique(pattern_gp_label))," patterns"))
  ##### calculate the number of cells we need to sample in each cell type
  set.seed(seed)
  # # temp.exprs <- exprs(eset.sub.split1)  #expression from first split of scRNA-seq data, genes on rows cells on cols
  temp.exprs <- RNA.counts
  temp.nct <- Sample_random  #temporary number of each type of cell
  true.p = sweep(temp.nct,1,rowSums(temp.nct),"/")  #divide the cell counts at each location by the number of cells at that location to get the true cell type proportions
  ##### Randomness in pbmclapply
  temp.pseudo = pbmclapply(1:nrow(temp.nct),function(isample){  #This creates the spatially resolved gene expression matrix
    temp.sample = temp.nct[isample,]  #take the row isample from the cell counts matrix temp.nct
    temp.sample.pseudo = sapply(ct.select,function(ict){
      temp.vec <- temp.exprs[,ct.id %in% ict]  #get the expression for cells for which their types are in a specific layer
      if(temp.nct[isample,ict] > ncol(temp.vec)){  #if the number of  cells in row isample column names ict (cell types that live in layer ict) is larger than the columns of the temporary expression vector, we'll need to sample with replacement
        sample.withRep = T
      } else {
        sample.withRep = F
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
  temp.pseudo <- round(temp.pseudo, digits = 0)
  return(list(pseudo.data = temp.pseudo, true.p = true.p, ntotal = ntotal,Sample_random = Sample_random))
}



