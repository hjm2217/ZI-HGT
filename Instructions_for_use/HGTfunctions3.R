library(CARD)
library(MuSiC)
library(LaplacesDemon)
library(loo)
library(rdist)
library(reshape2)
library(SummarizedExperiment)
library(ALRA)

###################################################
## Calculate the root mean square error for CARD ##
###################################################

rmse <- function(simProp, trueProp, scenario = 1){
  
  if (scenario == 2){
    simProp <- cbind(simProp, rep(0,nrow(simProp)))
    colnames(simProp)[6] <- "Neurons"
    simProp <- simProp[,c(1,6,2:5)]
  }
  
  if (scenario == 3){
    trueProp <- cbind(trueProp, rep(0, nrow(trueProp)))
    colnames(trueProp)[7] <- "Blood"
  }

  if (scenario == 4){
    trueProp <- cbind(trueProp, trueProp[, "Astrocytes"] + trueProp[, "Ependymal"])
    colnames(trueProp)[7] <- "Merged"
    trueProp <- trueProp[, c("Neurons", "Oligos", "Vascular", "Immune", "Merged")]
  }

  N <- nrow(simProp)
  K <- ncol(simProp)
  
  diffMat2 <- (trueProp - simProp)^2
  out <- 1/(N*K) * sqrt(sum(diffMat2))
  return(out)  
}

##############################################
## Logit, 1/(1-x), inverse logit, 1/(.75-x) ##
##############################################

logit <- function(x){
  out <- log( (x) / (1 - x) )
  return(out)
}

invlogit <- function(x){
  out <- 1 / (1 + exp(-x))
  return(out)
}

xover1minusx <- function(x){
  out <- x / (1-x)
  return(out)
}

xover.75minusx <- function(x){
  out <- x / (.75 - x)
  return(out)
}

###########################################################
## Transform the spatial gene expression data with a ZiP ##
###########################################################

transformX2H <- function(X, alpha0, alpha1, seed = 123456){
  set.seed(seed)
  n.zeros <- sum(colSums(X == 0))
  n.nonzeros <- sum(colSums(X > 0))
  H <- X  #leave any nonzero entries the same
  H[H > 0] <- (rgamma(n.nonzeros, shape = alpha1 + H[H > 0], rate = alpha1 + 1))
  H[H == 0] <- logit(rbeta(n.zeros, (alpha0*0.75), (0.25*alpha0+1))) #letting the prior be Beta(alpha-kappa, kappa), where kappa=0.25*alpha
  return(H)
}

##########################################################
## Get CARD to stop removing negative values in ST data ##
##########################################################

CARD_deconvolution_HGT <- function(CARD_object,markerCutoff=0.001){
  ct.select = CARD_object@info_parameters$ct.select
  ct.varname = CARD_object@info_parameters$ct.varname
  sample.varname = CARD_object@info_parameters$sample.varname
  cat(paste0("## create reference matrix from scRNASeq...\n"))
  sc_eset = CARD_object@sc_eset
  Basis_ref = createscRef(sc_eset, ct.select, ct.varname, sample.varname)
  Basis = Basis_ref$basis
  Basis = Basis[,colnames(Basis) %in% ct.select]
  Basis = Basis[,match(ct.select,colnames(Basis))]
  spatial_count = CARD_object@spatial_countMat
  commonGene = intersect(rownames(spatial_count),rownames(Basis))
  #### remove mitochondrial and ribosomal genes
  commonGene  = commonGene[!(commonGene %in% commonGene[grep("mt-",commonGene)])]
  cat(paste0("## Select Informative Genes! ...\n"))
  common = selectInfo_HGT(Basis,sc_eset,commonGene,ct.select,ct.varname,markerCutoff)
  Xinput = spatial_count
  rm(spatial_count)
  B = Basis
  rm(Basis)
  ##### match the common gene names
  Xinput = Xinput[order(rownames(Xinput)),]
  B = B[order(rownames(B)),]
  B = B[rownames(B) %in% common,]
  Xinput = Xinput[rownames(Xinput) %in% common,]
  ##### #DON'T# filter out non expressed genes or cells again
  ## Xinput = Xinput[rowSums(Xinput) > 0,]
  ## Xinput = Xinput[,colSums(Xinput) > 0]
  ##### normalize count data
  colsumvec = colSums(Xinput)
  Xinput_norm = sweep(Xinput,2,colsumvec,"/")
  B = B[rownames(B) %in% rownames(Xinput_norm),]    
  B = B[match(rownames(Xinput_norm),rownames(B)),]
  #### spatial location
  spatial_location = CARD_object@spatial_location
  spatial_location = spatial_location[rownames(spatial_location) %in% colnames(Xinput_norm),]
  spatial_location = spatial_location[match(colnames(Xinput_norm),rownames(spatial_location)),]
  
  ##### normalize the coordinates without changing the shape and relative position
  norm_cords = spatial_location[ ,c("x","y")]
  norm_cords$x = norm_cords$x - min(norm_cords$x)
  norm_cords$y = norm_cords$y - min(norm_cords$y)
  scaleFactor = max(norm_cords$x,norm_cords$y)
  norm_cords$x = norm_cords$x / scaleFactor
  norm_cords$y = norm_cords$y / scaleFactor
  ##### initialize the proportion matrix
  ED <- rdist(as.matrix(norm_cords))##Euclidean distance matrix
  ED <- as.matrix(ED)
  cat(paste0("## Deconvolution Starts! ...\n"))
  set.seed(20200107)
  Vint1 = as.matrix(rdirichlet(ncol(Xinput_norm), rep(10,ncol(B))))
  colnames(Vint1) = colnames(B)
  rownames(Vint1) = colnames(Xinput_norm)
  b = rep(0,length(ct.select))
  ###### parameters that need to be set
  isigma = 0.1 ####construct Gaussian kernel with the default scale /length parameter to be 0.1
  epsilon = 1e-04  #### convergence epsion 
  phi = c(0.01,0.1,0.3,0.5,0.7,0.9,0.99) #### grided values for phi
  kernel_mat <- exp(-ED^2 / (2 * isigma^2))
  diag(kernel_mat) <- 0
  rm(ED)
  rm(Xinput)
  rm(norm_cords)
  gc()
  ###### scale the Xinput_norm and B to speed up the convergence. 
  mean_X = mean(Xinput_norm)
  mean_B = mean(B)
  Xinput_norm = Xinput_norm * 1e-01 / mean_X
  B = B * 1e-01 / mean_B
  gc()
  ResList = list()
  Obj = c()
  for(iphi in 1:length(phi)){
    res = CARDref(
      XinputIn = as.matrix(Xinput_norm),
      UIn = as.matrix(B),
      WIn = kernel_mat, 
      phiIn = phi[iphi],
      max_iterIn =1000,
      epsilonIn = epsilon,
      initV = Vint1,
      initb = rep(0,ncol(B)),
      initSigma_e2 = 0.1, 
      initLambda = rep(10,length(ct.select)))
    rownames(res$V) = colnames(Xinput_norm)
    colnames(res$V) = colnames(B)
    ResList[[iphi]] = res
    Obj = c(Obj,res$Obj)
  }
  Optimal = which(Obj == max(Obj))
  Optimal = Optimal[length(Optimal)] #### just in case if there are two equal objective function values
  OptimalPhi = phi[Optimal]
  OptimalRes = ResList[[Optimal]]
  cat(paste0("## Deconvolution Finish! ...\n"))
  CARD_object@info_parameters$phi = OptimalPhi
  CARD_object@Proportion_CARD = sweep(OptimalRes$V,1,rowSums(OptimalRes$V),"/")
  CARD_object@algorithm_matrix = list(B = B * mean_B / 1e-01, Xinput_norm = Xinput_norm * mean_X / 1e-01, Res = OptimalRes)
  CARD_object@spatial_location = spatial_location
  return(CARD_object)
}


##################################################################################################################
## Require marker genes to have average expression of at least some value in one cell type in the scRNAseq data ##
##################################################################################################################

selectInfo_HGT <- function(Basis,sc_eset,commonGene,ct.select,ct.varname,markerCutoff){
#### log2 mean fold change >0.5
gene1 = lapply(ct.select,function(ict){
rest = rowMeans(Basis[,colnames(Basis) != ict])
FC = log((Basis[,ict] + 1e-06)) - log((rest + 1e-06))
rownames(Basis)[FC > 1.25 & Basis[,ict] > markerCutoff]
})
gene1 = unique(unlist(gene1))
gene1 = intersect(gene1,commonGene)
counts = assays(sc_eset)$counts
counts = counts[rownames(counts) %in% gene1,]
##### only check the cell type that contains at least 2 cells
ct.select = names(table(colData(sc_eset)[,ct.varname]))[table(colData(sc_eset)[,ct.varname]) > 1]
sd_within = sapply(ct.select,function(ict){
  temp = counts[,colData(sc_eset)[,ct.varname] == ict]
  apply(temp,1,var) / apply(temp,1,mean)
  })
##### remove the outliers that have high dispersion across cell types
gene2 = rownames(sd_within)[apply(sd_within,1,mean,na.rm = T) < quantile(apply(sd_within,1,mean,na.rm = T),prob = 0.99,na.rm = T)]
return(gene2)
}


##########################################
## Run the HGT over CARD in Simulations ##
##########################################

#We need to loop 50ish times using the CARD method and different versions of H
SimOverCARD <- function(sc_count,  
                        sc_meta,  
                        spatial_count,
                        spatial_location,
                        ct.varname,
                        ct.select, 
                        sample.varname,
                        minCountGene,
                        minCountSpot,
                        truePropMat,
                        scenario = 1,
                        alpha0,
                        alpha1,
                        n = 100,
                        doWAIC = TRUE,
                        seed = 1){
  time.start <- proc.time()
  #things to spit out: original spatial gene expression matrix, original cell-type proportion matrix, list of all H matrices, list of cell-type proportion matrices, true cell type proportion matrix, average cell-type proportions across runs, RMSE for each of 50 tries, RMSE for average cell-type proportion, RMSE for original CARD, alpha, number of runs, original seed 
  #run CARD originally
  CARDobjOG <- createCARDObject(sc_count, sc_meta, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene, minCountSpot)
  CARDobjOG <- CARD_deconvolution(CARDobjOG)
  
  #and record the original spatial gene expression matrix, the true and CARD cell type proportion matrices, and the CARD RMSE
  #spatialCountList <- list(original_spatial_count = spatial_count), can't afford to save
  cellTypePropMatList <- list(true_ct_prop_mat = truePropMat, CARD_ct_prop_mat = CARDobjOG@Proportion_CARD)
  rmseOG <- rmse(CARDobjOG@Proportion_CARD, truePropMat, scenario = scenario)
  rmseList <- list(CARD_rmse = rmseOG)
  cat("done with original CARD\n")
  
  #Initialize the list of transformations and RMSE for each run of MCMC
  TcellTypePropMatList <- list()
  TrmseList <- list()
  rm(CARDobjOG)
  
  logliklist <- list()
  
  for(i in 1:n){
    Tspatial_count <- transformX2H(spatial_count, alpha0 = alpha0, alpha1 = alpha1, seed = (seed + (i-1)*5000000))  #transform spatial counts matrix

    
    CARDobjtemp <- createCARDObject(sc_count, sc_meta, Tspatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene, minCountSpot)  #create CARD object with transformed spatial counts matrix
    CARDobjtemp <- CARD_deconvolution_HGT(CARDobjtemp)  #run CARD with transformed spatial counts
    rmsetemp <- rmse(CARDobjtemp@Proportion_CARD, truePropMat, scenario = scenario)  #find RMSE for this run
    
    #TspatialCountList[[i]] <- Tspatial_count
    TcellTypePropMatList[[i]] <- CARDobjtemp@Proportion_CARD
    TrmseList[[i]] <- rmsetemp

   # We want to record the WAIC to help choose hyperparameters.  WAIC = -2 (lpd - p_WAIC), and we can calculate this in R with the WAIC function, which requires the loglikelihood in an N x S matrix, where N is total data points and S the number of MCMC samples
  if (doWAIC){

     # To find the WAIC, we'll need 1) the original data (restricted to cell-type informative genes and genes/locations that pass QC)
     Xinput_norm <- CARDobjtemp@algorithm_matrix$Xinput_norm
     spatial_count_restricted <- spatial_count[rownames(spatial_count) %in% rownames(Xinput_norm), colnames(spatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
     spatial_count_restricted <- spatial_count_restricted[ match(rownames(Xinput_norm), rownames(spatial_count_restricted)), match(colnames(Xinput_norm), colnames(spatial_count_restricted))]

     # 2) the transformed data under the same restrictions
     Tspatial_count_restricted <- Tspatial_count[rownames(Tspatial_count) %in% rownames(Xinput_norm), colnames(Tspatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
     Tspatial_count_restricted <- Tspatial_count_restricted[ match(rownames(Xinput_norm), rownames(Tspatial_count_restricted)), match(colnames(Xinput_norm), colnames(Tspatial_count_restricted))]

     # 3) the things we need from the CARD model, B, V, and sigma_e^2
     B <- CARDobjtemp@algorithm_matrix$B
     V <- CARDobjtemp@algorithm_matrix$Res$V
     sigma_e2 <- CARDobjtemp@algorithm_matrix$Res$sigma_e2

     # 4) the latent mean of the CARD model, X^(2) = BV'
     X2 <- B %*% t(V)

     # 5) the log-likelihood at each data point
     #logliklist <- list()
     logliktemp <- matrix(0, dim(spatial_count_restricted)[1], dim(spatial_count_restricted)[2])
     Idx0 <- which(spatial_count_restricted == 0)
     Idxg0 <- which(spatial_count_restricted > 0)
     rownames(logliktemp) <- rownames(spatial_count_restricted)
     colnames(logliktemp) <- colnames(spatial_count_restricted)

     #for the zero data, we have a point mass at 0 * point mass at X^(B)_ij * dbern(0, p = invlogit( X^(B)_ij ) * dnorm( X^(1)_ij, mean = X^(2)_ij, var = sigma_e2) and X^(B) = X^(1)
     logliktemp[Idx0] <- dbern(0, prob = invlogit(Tspatial_count_restricted[Idx0]), log=TRUE) + dnorm(Tspatial_count_restricted[Idx0], mean = X2[Idx0], sd = sqrt(sigma_e2), log=TRUE)
     logliktemp[Idxg0] <- dpois(spatial_count_restricted[Idxg0], lambda = Tspatial_count_restricted[Idxg0], log=TRUE) + dnorm(Tspatial_count_restricted[Idxg0], mean = X2[Idxg0], sd = sqrt(sigma_e2), log=TRUE)
     logliklist[[i]] <- logliktemp

     rm(logliktemp)
  }

    rm(CARDobjtemp)
    rm(Tspatial_count)
    rm(rmsetemp)
    gc()
    cat("done with transformed run", i, "\n")
  }
  
  #now work with average cell-type proportions from our different runs
  #later provide 95% credible intervals on proportions, but I won't worry about it now
  cellTypePropMatAve <- Reduce("+", TcellTypePropMatList) / n
  rmseAve <- rmse(cellTypePropMatAve, truePropMat, scenario = scenario)
  
  cellTypePropMatList <- c(cellTypePropMatList, Ave_ct_prop_mat = cellTypePropMatAve, transformed_ct_prop_mat_list = TcellTypePropMatList)
  rmseList <- c(rmseList, Ave_rmse = rmseAve, transformed_rmse_list = TrmseList)
  
  #Get the WAIC
  if (doWAIC){

     # We can only look at locations and genes that each iteration has in common (luckily it's rare that any iteration loses locations or genes)
     common_rows <- rownames(logliklist[[1]])
     common_cols <- colnames(logliklist[[1]])
     for (i in 1:n){
	  common_rows <- intersect(common_rows, rownames(logliklist[[i]]))
	  common_cols <- intersect(common_cols, colnames(logliklist[[i]]))
     }

     totaldatapoints <- length(common_rows)*length(common_cols)
     allLogLikelihood <- matrix(0, nrow = totaldatapoints, ncol = n)
     cat(length(common_rows), "common genes and ", length(common_cols), "common locations\n")
     cat(dim(logliklist[[1]])[1], "genes in iteration 1 and ", dim(logliklist[[1]])[2], "common locations\n\n")  

     for (i in 1:n){
     	tmploglik <- logliklist[[i]]
	tmploglik <- tmploglik[ rownames(tmploglik) %in% common_rows, colnames(tmploglik) %in% common_cols ]
	allLogLikelihood[, i] <- tmploglik
	rm(tmploglik)
	gc()
     }

    WAICout <- waic(allLogLikelihood)$estimates["waic", "Estimate"]  #3 seconds-ish for 10 runs
  } else {
    WAICout <- NULL
  }
  
  time.end <- proc.time()
  time_elapsed = time.end[3] - time.start[3]
  out <- list(cell_type_proportion_matrices = cellTypePropMatList,
              rmse = rmseList,
              WAIC = WAICout,
              alpha0 = alpha0,
              alpha1 = alpha1,
              n = n,
              original_seed = seed,
              scenario = scenario,
              time = time_elapsed)
  return(out)
  cat("All done!\n")
}

###########################
## Run the HGT over CARD ##
###########################

#We need to loop 50ish times using the CARD method and different versions of H
HGTplusCARD <- function(sc_count,
                        sc_meta,
                        spatial_count,
                        spatial_location,
                        ct.varname,
                        ct.select,
                        sample.varname,
                        minCountGene,
                        minCountSpot,
                        alpha0,
                        alpha1,
                        n = 100,
			                  markerCutoff = 0.05,
                        doWAIC = TRUE,
                        seed = 1){
  time.start <- proc.time()
  #things to spit out: original spatial gene expression matrix, list of all X^(1) matrices, list of cell-type proportion matrices, alpha0, alpha1, WAIC, average cell-type proportions
  
  #Initialize the list of transformations, log-likelihoods, and Fisher information matrices
  TcellTypePropMatList <- list()
  logliklist <- list()
  inverseFIMatList <- list() 

  # Will need Gaussian kernel for FI calculation later
  norm_cords = spatial_location[ ,c("x","y")] # CARD uses normalized coordinates
  norm_cords$x = norm_cords$x - min(norm_cords$x)
  norm_cords$y = norm_cords$y - min(norm_cords$y)
  scaleFactor = max(norm_cords$x,norm_cords$y)
  norm_cords$x = norm_cords$x / scaleFactor
  norm_cords$y = norm_cords$y / scaleFactor
  distance_matrix <- rdist(as.matrix(norm_cords)) # Euclidean distance matrix
  sigma_W = 0.1 # Construct Gaussian kernel with the default scale /length parameter to be 0.1
  W_mat <- as.matrix(exp(-distance_matrix^2 / (2 * sigma_W^2)))
  diag(W_mat) <- 0
  D_mat <- diag(rowSums(W_mat)) # D = diag(W_1+, ..., W_n+)

  pb <- txtProgressBar(min = 0, max = n, style = 3) # progress bar
  for (i in 1:n){
    Tspatial_count <- transformX2H(spatial_count, alpha0 = alpha0, alpha1 = alpha1, seed = (seed + (i-1)*5000000))  #transform spatial counts matrix
    invisible(capture.output(CARDobjtemp <- createCARDObject_HGT(sc_count, sc_meta, Tspatial_count, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene, minCountSpot)))  #create CARD object with transformed spatial counts matrix
    invisible(capture.output(CARDobjtemp <- CARD_deconvolution_HGT(CARDobjtemp,markerCutoff=markerCutoff)))  #run CARD with transformed spatial counts
    
    TcellTypePropMatList[[i]] <- CARDobjtemp@Proportion_CARD
    Vtemp <- CARDobjtemp@algorithm_matrix$Res$V
    # We can sometimes get very small negative numbers instead of 0s in the cell-type proportions, so we'll need to zero those out and rescale the rows
    TcellTypePropMatList[[i]][TcellTypePropMatList[[i]] < 0] <- 0
    TcellTypePropMatList[[i]] <- sweep(TcellTypePropMatList[[i]], 1, rowSums(TcellTypePropMatList[[i]]), FUN = "/")
    Vtemp[Vtemp < 0] <- 0
    
   # We want to record the WAIC to help choose hyperparameters.  WAIC = -2 (lpd - p_WAIC), and we can calculate this in R with the WAIC function, which requires the loglikelihood in an N x S matrix, where N is total data points and S the number of MCMC samples
    if (doWAIC){

      # To find the WAIC, we'll need 1) the original data (restricted to cell-type informative genes and genes/locations that pass QC)
      Xinput_norm <- CARDobjtemp@algorithm_matrix$Xinput_norm
      spatial_count_restricted <- spatial_count[rownames(spatial_count) %in% rownames(Xinput_norm), colnames(spatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
      spatial_count_restricted <- spatial_count_restricted[ match(rownames(Xinput_norm), rownames(spatial_count_restricted)), match(colnames(Xinput_norm), colnames(spatial_count_restricted))]

      # 2) the transformed data under the same restrictions
      Tspatial_count_restricted <- Tspatial_count[rownames(Tspatial_count) %in% rownames(Xinput_norm), colnames(Tspatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
      Tspatial_count_restricted <- Tspatial_count_restricted[ match(rownames(Xinput_norm), rownames(Tspatial_count_restricted)), match(colnames(Xinput_norm), colnames(Tspatial_count_restricted))]

      # 3) the things we need from the CARD model, B, V, and sigma_e^2
      B <- CARDobjtemp@algorithm_matrix$B
      V <- Vtemp
      sigma_e2 <- CARDobjtemp@algorithm_matrix$Res$sigma_e2

      # 4) the latent mean of the CARD model, X^(2) = BV'
      X2 <- B %*% t(V)

      # 5) the log-likelihood at each data point
      #logliklist <- list()
      logliktemp <- matrix(0, dim(spatial_count_restricted)[1], dim(spatial_count_restricted)[2])
      Idx0 <- which(spatial_count_restricted == 0)
      Idxg0 <- which(spatial_count_restricted > 0)
      rownames(logliktemp) <- rownames(spatial_count_restricted)
      colnames(logliktemp) <- colnames(spatial_count_restricted)

      #for the zero data, we have a point mass at 0 * point mass at X^(B)_ij * dbern(0, p = invlogit( X^(B)_ij ) * dnorm( X^(1)_ij, mean = X^(2)_ij, var = sigma_e2) and X^(B) = X^(1)
      logliktemp[Idx0] <- dbern(0, prob = invlogit(Tspatial_count_restricted[Idx0]), log=TRUE) + dnorm(Tspatial_count_restricted[Idx0], mean = X2[Idx0], sd = sqrt(sigma_e2), log=TRUE)
      logliktemp[Idxg0] <- dpois(spatial_count_restricted[Idxg0], lambda = Tspatial_count_restricted[Idxg0], log=TRUE) + dnorm(Tspatial_count_restricted[Idxg0], mean = X2[Idxg0], sd = sqrt(sigma_e2), log=TRUE)
      logliklist[[i]] <- logliktemp

      rm(logliktemp)
    }

  # We need to construct the Fisher information matrix to be able to (approximately) determine the pointwise Bayesian credible intervals for V_j,k (loc j cell type k)
  # The FI is 1/(sigma_e^2) * B_k'B_k + 1/(lambda_k) * L_j,j, where B_k is the kth column of B, L = D - phi*W, D = diag(W_1+, ..., W_n+), and W is the gaussian kernel matrix, W_i,j = exp(distance/2*0.1^2)
  sigma_e2 <- CARDobjtemp@algorithm_matrix$Res$sigma_e2
  B <- CARDobjtemp@algorithm_matrix$B *1e-01 / mean(CARDobjtemp@algorithm_matrix$B) # rescale B to be on the same scale as in CARD's deconvolution
  lambda <- unlist(lapply(CARDobjtemp@algorithm_matrix$Res$lambda, function(x) min(x, 25))) # lambda parameter can get very big, but isn't super meaningful past a certain point, so we'll cap it
  # We have the distance matrix, W_mat, and the diagonal matrix D_mat already
  L_mat <- D_mat - CARDobjtemp@info_parameters$phi * W_mat

  V <- CARDobjtemp@Proportion_CARD
  ((t(V[,2] - mean(V[,2])) %*% L_mat %*% (V[,2] - mean(V[,2])) ) / 2 + 500) / (2 + 500)

  FI <- matrix(0, nrow = nrow(CARDobjtemp@Proportion_CARD), ncol = ncol(CARDobjtemp@Proportion_CARD))
  rownames(FI) <- rownames(CARDobjtemp@Proportion_CARD); colnames(FI) <- colnames(CARDobjtemp@Proportion_CARD)
  for (j in 1:nrow(FI)){
    for (k in 1:ncol(FI)){
      FI[j, k] <- 1/(sigma_e2) * (t(B[, k]) %*% B[, k]) + 1/lambda[k] * L_mat[j, j]
    }
  }
  inverseFI <- 1/FI # Take the pointwise inverse of the Fisher information matrix
  inverseFIMatList[[i]] <- inverseFI
  
    rm(CARDobjtemp)
    rm(Tspatial_count)
    gc()
    setTxtProgressBar(pb, i) # increment progress bar
  }
  close(pb)

  # now work with the cell-type proportions from our different runs
  # We'll provide the mean, median, and 95% Bayesian credible intervals
  # Mean
  cellTypePropMatAve <- Reduce("+", TcellTypePropMatList) / n

  # Median
  medianArray <- array(unlist(TcellTypePropMatList), c(dim(TcellTypePropMatList[[1]]), length(TcellTypePropMatList)))
  medianMat <- apply(medianArray, 1:2, median)
  cellTypePropMatMedian <- sweep(medianMat, 1, rowSums(medianMat), FUN = "/")
  rownames(cellTypePropMatMedian) <- rownames(cellTypePropMatAve)
  colnames(cellTypePropMatMedian) <- colnames(cellTypePropMatAve)

  # Find the variance of each cell type proportion across the runs (Var(E[V_j,k]))
  var_expected_V <- apply(medianArray, 1:2, var)
  mean_inverseFI <- Reduce("+", inverseFIMatList) / n  
  xi_squared <- mean_inverseFI + var_expected_V # xi^2 = CARD var estimate + HGT var estimate

  # 95% Pointwise Bayesian Credible Intervals
  # For each element in the cellTypePropMatAve, we add 1.96sqrt(xi_squared) to get the upper bound and subtract it to get the lower bound, but bound at 0 and 1
  lowerMat <- cellTypePropMatAve - 1.96 * sqrt(xi_squared)
  lowerMat[lowerMat < 0] <- 0
  upperMat <- cellTypePropMatAve + 1.96 * sqrt(xi_squared)
  upperMat[upperMat > 1] <- 1


  cellTypePropMatList <- list(mean_ct_prop = cellTypePropMatAve, median_ct_prop = cellTypePropMatMedian, lower_95_prop = lowerMat, upper_95_prop = upperMat, transformed_ct_prop_mat_list = TcellTypePropMatList, CARD_var = mean_inverseFI, HGT_var = var_expected_V, total_var = xi_squared)

   #Get the WAIC
  if (doWAIC){

     # We can only look at locations and genes that each iteration has in common (luckily it's rare that any iteration loses locations or genes)
     common_rows <- rownames(logliklist[[1]])
     common_cols <- colnames(logliklist[[1]])
     for (i in 1:n){
          common_rows <- intersect(common_rows, rownames(logliklist[[i]]))
          common_cols <- intersect(common_cols, colnames(logliklist[[i]]))
     }

     totaldatapoints <- length(common_rows)*length(common_cols)
     allLogLikelihood <- matrix(0, nrow = totaldatapoints, ncol = n)
     cat(length(common_rows), "common genes and ", length(common_cols), "common locations\n")
     cat(dim(logliklist[[1]])[1], "genes in iteration 1 and ", dim(logliklist[[1]])[2], "common locations\n\n")

     for (i in 1:n){
	    tmploglik <- logliklist[[i]]
        tmploglik <- tmploglik[ rownames(tmploglik) %in% common_rows, colnames(tmploglik) %in% common_cols ]
        allLogLikelihood[, i] <- tmploglik
        rm(tmploglik)
        gc()
     }

    WAICout <- waic(allLogLikelihood)$estimates["waic", "Estimate"]  #3 seconds-ish for 10 runs
  } else {
    WAICout <- NULL
  }

  time.end <- proc.time()
  time_elapsed = time.end[3] - time.start[3]
  out <- list(cell_type_proportion_matrices = cellTypePropMatList,
              WAIC = WAICout,
              alpha0 = alpha0,
              alpha1 = alpha1,
              n = n,
              original_seed = seed,
              time = time_elapsed)
  return(out)
  cat("All done!\n")
}

###############################
# Relabel the OSCC cell types #
###############################
OSCC_cell_types <- list(
  "myofibroblast" = "Myofibroblasts",
  "cancer cell" = "Cancer Cells",
  "B cell" = "B Cells",
  "ecm-myCAF" = "ECM-myCAFs",                       
  "Intermediate fibroblast" = "Intermediate Fibroblasts",
  "detox-iCAF" = "Detox-iCAFs",
  "macrophage" = "Macrophages",
  "endothelial" = "Endothelial Cells",
  "dendritic " = "Dendritic Cells",
  "mast" = "Mast Cells",                            
  "conventional CD4+ T-helper cells" = "CD4+ Th Cells",
  "cytotoxic CD8+ T " = "Cytotoxic CD8+ T Cells",               
  "Tregs" = "Tregs",
  "cytotoxic CD8+ T exhausted" = "Cytotoxic CD8+ Exhausted T Cells"  
)
OSCC_labeller <- function(value){
  return(OSCC_cell_types[[value]])
}

#######################################################
# Plot the absolute proportions of all the cell types #
#######################################################

HGT.CARD.visualize.prop <- function(proportion,spatial_location,ct.visualize = ct.visualize,colors = c("lightblue","lightyellow","red"),NumCols, pointSize = 3.0){
  if(is.null(colors)){
    colors = c("lightblue","lightyellow","red")
  }else{
    colors = colors
  }
  res_CARD = as.data.frame(proportion)
  res_CARD = res_CARD[,order(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  ct.select = ct.visualize
  res_CARD = res_CARD[,ct.select]
  #if(!is.null(ncol(res_CARD))){
  #  res_CARD_scale = as.data.frame(apply(res_CARD,2,function(x){
  #    (x - min(x)) / (max(x) - min(x))
  #  } ))}else{
  #    res_CARD_scale = as.data.frame((res_CARD - min(res_CARD)) / (max(res_CARD) - min(res_CARD)))
  #    colnames(res_CARD_scale) = ct.visualize
  #  }
  res_CARD$x = as.numeric(location$x)
  res_CARD$y = as.numeric(location$y)
  mData = melt(res_CARD,id.vars = c("x","y"))
  colnames(mData)[3] <- "Cell_Type"
  b = c(0,1)
  p = suppressMessages(ggplot(mData, aes(x, y)) + 
                         geom_point(aes(colour = value),size = pointSize) +
                         scale_color_gradientn(colours = colors, limits = c(0,1)) + 
                         #scale_color_viridis_c(option = 2)+
                         scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
                         facet_wrap(~Cell_Type,ncol = NumCols, labeller = labeller(Cell_Type = label_wrap_gen(20)))+ 
                         coord_fixed()+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               #legend.position=c(0.14,0.76),
                               panel.background = element_blank(),
                               plot.background = element_blank(), #element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 14,face="bold"),
                               legend.text=element_text(size = 11),
                               strip.text = element_text(size = 12,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm')))
  return(p)
}

########################################################
# Plot two-tone plots (cell type prop > cutoff vs. not #
########################################################

HGT.CARD.two.tone <- function(proportion, spatial_location, ct.visualize = ct.visualize, colors = c("grey70","darkgreen"), NumCols, pointSize = 3.0, cutoff = 0.05){
  if(is.null(colors)){
    colors = c("grey70","darkgreen")
  }else{
    colors = colors
  }
  props = as.data.frame(proportion)
  res_CARD = apply(props, 2,function(x) ifelse(x >= cutoff, 1, 0))
  res_CARD = as.data.frame(res_CARD)
  res_CARD = res_CARD[,order(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  ct.select = ct.visualize
  res_CARD = res_CARD[,ct.select]
  #if(!is.null(ncol(res_CARD))){
  res_CARD$x = as.numeric(location$x)
  res_CARD$y = as.numeric(location$y)
  mData = melt(res_CARD,id.vars = c("x","y"))
  colnames(mData)[3] <- "Cell_Type"
  mData$value <- factor(mData$value)
  b = c(0,1)
  p = suppressMessages(ggplot(mData, aes(x, y)) +
                         geom_point(aes(colour = value),size = pointSize) +
                         scale_color_manual(values = colors) +
                         #scale_color_viridis_c(option = 2)+
                         scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+
                         facet_wrap(~Cell_Type,ncol = NumCols, labeller = labeller(Cell_Type = label_wrap_gen(20))) +
                         coord_fixed() + labs(colour = paste0("Proportion > ", cutoff)) + 
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               #legend.position=c(0.14,0.76),
                               panel.background = element_blank(),
                               plot.background = element_blank(), #element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 14,face="bold"),
                               legend.text=element_text(size = 11),
                               strip.text = element_text(size = 12,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm')))
  return(p)
}


################################
# Plot the Heterogeniety Score #
################################

HGT.visualize.heterogeneity <- function(het_df,spatial_location,colors = c("white","darkgreen"), pointSize = 3.0){

  res_het = as.data.frame(het_df)
  rownames(res_het) <- res_het$location
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_het)==rownames(location))!= nrow(res_het)){
    stop("The rownames of heterogeneity data does not match with the rownames of spatial location data")
  }
  res_het = res_het[,c("x","y","het_score")]
  res_het = res_het[rownames(location), ]
  p = suppressMessages(ggplot(res_het, aes(x, y)) +
                         geom_point(aes(colour = het_score),size = pointSize) +
                         scale_color_gradientn(colours = colors) +
                         #scale_color_viridis_c(option = 2)+
                         scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+
                         coord_fixed()+ labs(color = "Heterogeneity\nScore") + 
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               #legend.position=c(0.14,0.76),
                               panel.background = element_blank(),
                               plot.background = element_blank(), #element_blank(),
                               panel.border = element_blank(),  #element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 14,face="bold"),
                               legend.text=element_text(size = 11),
                               strip.text = element_text(size = 12,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm')))
  return(p)
}


####################################
# Plot the Pathologist Annotations #
####################################

HGT.plot.patho.annot <- function(patho_annot_df, colors, pointSize = 3.0){

#  res_het = as.data.frame(het_df)
#  rownames(res_het) <- res_het$location
#  location = as.data.frame(spatial_location)
#  if(sum(rownames(res_het)==rownames(location))!= nrow(res_het)){
#    stop("The rownames of heterogeneity data does not match with the rownames of spatial location data")
#  }
#  res_het = res_het[,c("x","y","het_score")]
#  res_het = res_het[rownames(location), ]
  p = suppressMessages(ggplot(patho_annot_df, aes(x, y)) +
                         geom_point(aes(colour = patho_annot_f),size = pointSize) +
                         scale_colour_manual(values = colors) +
                         scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+
                         coord_fixed()+ labs(color = "Pathologist Annotation") +
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               #legend.position=c(0.14,0.76),
                               panel.background = element_blank(),
                               plot.background = element_blank(), #element_blank(),
                               panel.border = element_blank(),  #element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 14,face="bold"),
                               legend.text=element_text(size = 11),
                               strip.text = element_text(size = 12,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm')))
  return(p)
}


###################################################
## Run the HGT over CARD in the OSCC Simulations ##
###################################################

#We need to loop 50ish times using the CARD method and different versions of H
SimOverCARD_OSCC <- function(sc_count,  
                             sc_meta,  
                             spatial_count,
                             spatial_location,
                             ct.varname,
                             ct.select, 
                             sample.varname,
                             minCountGene,
                             minCountSpot,
                             truePropMat,
                             scenario = 1,
                             alpha0,
                             alpha1,
                             n = 100,
			     markerCutoff = 0.05,
                             doWAIC = TRUE,
                             seed = 1){
  time.start <- proc.time()
  #things to spit out: original spatial gene expression matrix, original cell-type proportion matrix, list of all H matrices, list of cell-type proportion matrices, true cell type proportion matrix, average cell-type proportions across runs, RMSE for each of 50 tries, RMSE for average cell-type proportion, RMSE for original CARD, alpha, number of runs, original seed 
  #run CARD originally
  CARDobjOG <- createCARDObject(sc_count, sc_meta, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene = 100, minCountSpot = 5)
  CARDobjOG <- CARD_deconvolution(CARDobjOG)
  
  #and record the original spatial gene expression matrix, the true and CARD cell type proportion matrices, and the CARD RMSE
  #spatialCountList <- list(original_spatial_count = spatial_count), can't afford to save
  rmseOG <- rmse(CARDobjOG@Proportion_CARD, truePropMat, scenario = scenario)
  #rmseList <- list(CARD_rmse = rmseOG)
  cat("done with original CARD\n")
  
  #Initialize the list of transformations and RMSE for each run of MCMC
  TcellTypePropMatList <- list()
  TrmseList <- list()
    
  logliklist <- list()
  
  for(i in 1:n){
    Tspatial_count <- transformX2H(spatial_count, alpha0 = alpha0, alpha1 = alpha1, seed = (seed + (i-1)*5000000))  #transform spatial counts matrix
    
    CARDobjtemp <- createCARDObject_HGT(sc_count, sc_meta, Tspatial_count, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene, minCountSpot)  #create CARD object with transformed spatial counts matrix
    CARDobjtemp <- CARD_deconvolution_HGT(CARDobjtemp,markerCutoff=markerCutoff)  #run CARD with transformed spatial counts
    
    TcellTypePropMatList[[i]] <- CARDobjtemp@Proportion_CARD
    Vtemp <- CARDobjtemp@algorithm_matrix$Res$V
    # We can sometimes get very small negative numbers instead of 0s in the cell-type proportions, so we'll need to zero those out and rescale the rows
    TcellTypePropMatList[[i]][TcellTypePropMatList[[i]] < 0] <- 0
    TcellTypePropMatList[[i]] <- sweep(TcellTypePropMatList[[i]], 1, rowSums(TcellTypePropMatList[[i]]), FUN = "/")
    Vtemp[Vtemp < 0] <- 0
    
    rmsetemp <- rmse(TcellTypePropMatList[[i]], truePropMat, scenario = scenario)  #find RMSE for this run
    TrmseList[[i]] <- rmsetemp
    
    # We want to record the WAIC to help choose hyperparameters.  WAIC = -2 (lpd - p_WAIC), and we can calculate this in R with the WAIC function, which requires the loglikelihood in an N x S matrix, where N is total data points and S the number of MCMC samples
    if (doWAIC){
      
      # To find the WAIC, we'll need 1) the original data (restricted to cell-type informative genes and genes/locations that pass QC)
      Xinput_norm <- CARDobjtemp@algorithm_matrix$Xinput_norm
      spatial_count_restricted <- spatial_count[rownames(spatial_count) %in% rownames(Xinput_norm), colnames(spatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
      spatial_count_restricted <- spatial_count_restricted[ match(rownames(Xinput_norm), rownames(spatial_count_restricted)), match(colnames(Xinput_norm), colnames(spatial_count_restricted))]
      
      # 2) the transformed data under the same restrictions
      Tspatial_count_restricted <- Tspatial_count[rownames(Tspatial_count) %in% rownames(Xinput_norm), colnames(Tspatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
      Tspatial_count_restricted <- Tspatial_count_restricted[ match(rownames(Xinput_norm), rownames(Tspatial_count_restricted)), match(colnames(Xinput_norm), colnames(Tspatial_count_restricted))]
      
      # 3) the things we need from the CARD model, B, V, and sigma_e^2
      B <- CARDobjtemp@algorithm_matrix$B
      V <- Vtemp
      sigma_e2 <- CARDobjtemp@algorithm_matrix$Res$sigma_e2
      
      # 4) the latent mean of the CARD model, X^(2) = BV'
      X2 <- B %*% t(V)
      
      # 5) the log-likelihood at each data point
      #logliklist <- list()
      logliktemp <- matrix(0, dim(spatial_count_restricted)[1], dim(spatial_count_restricted)[2])
      Idx0 <- which(spatial_count_restricted == 0)
      Idxg0 <- which(spatial_count_restricted > 0)
      rownames(logliktemp) <- rownames(spatial_count_restricted)
      colnames(logliktemp) <- colnames(spatial_count_restricted)
      
      #for the zero data, we have a point mass at 0 * point mass at X^(B)_ij * dbern(0, p = invlogit( X^(B)_ij ) * dnorm( X^(1)_ij, mean = X^(2)_ij, var = sigma_e2) and X^(B) = X^(1)
      logliktemp[Idx0] <- dbern(0, prob = invlogit(Tspatial_count_restricted[Idx0]), log=TRUE) + dnorm(Tspatial_count_restricted[Idx0], mean = X2[Idx0], sd = sqrt(sigma_e2), log=TRUE)
      logliktemp[Idxg0] <- dpois(spatial_count_restricted[Idxg0], lambda = Tspatial_count_restricted[Idxg0], log=TRUE) + dnorm(Tspatial_count_restricted[Idxg0], mean = X2[Idxg0], sd = sqrt(sigma_e2), log=TRUE)
      logliklist[[i]] <- logliktemp
      
      rm(logliktemp)
    }
    
    rm(CARDobjtemp)
    rm(Tspatial_count)
    rm(rmsetemp)
    gc()
    cat("done with transformed run", i, "\n")
  }
  
 
  #now work with average cell-type proportions from our different runs
  #later provide 95% credible intervals on proportions, but I won't worry about it now
  cellTypePropMatAve <- Reduce("+", TcellTypePropMatList) / n
  rmseAve <- rmse(cellTypePropMatAve, truePropMat, scenario = scenario)
  
  medianArray <- array(unlist(TcellTypePropMatList), c(dim(TcellTypePropMatList[[1]]), length(TcellTypePropMatList)))
  medianMat <- apply(medianArray, 1:2, median)
  cellTypePropMatMedian <- sweep(medianMat, 1, rowSums(medianMat), FUN = "/")
  rownames(cellTypePropMatMedian) <- rownames(cellTypePropMatAve)
  colnames(cellTypePropMatMedian) <- colnames(cellTypePropMatAve)
  rmseMedian <- rmse(cellTypePropMatMedian, truePropMat, scenario = scenario)
  
  cellTypePropMatList <- list(true_ct_prop_mat = truePropMat, CARD_ct_prop_mat = CARDobjOG@Proportion_CARD, mean_ct_prop = cellTypePropMatAve, median_ct_prop = cellTypePropMatMedian, transformed_ct_prop_mat_list = TcellTypePropMatList)
  rmseList <- list(CARD_rmse = rmseOG, Ave_rmse = rmseAve, Median_rmse = rmseMedian, transformed_rmse_list = TrmseList)

  #Get the WAIC
  if (doWAIC){
    
    # We can only look at locations and genes that each iteration has in common (luckily it's rare that any iteration loses locations or genes)
    common_rows <- rownames(logliklist[[1]])
    common_cols <- colnames(logliklist[[1]])
    for (i in 1:n){
      common_rows <- intersect(common_rows, rownames(logliklist[[i]]))
      common_cols <- intersect(common_cols, colnames(logliklist[[i]]))
    }
    
    totaldatapoints <- length(common_rows)*length(common_cols)
    allLogLikelihood <- matrix(0, nrow = totaldatapoints, ncol = n)
    cat(length(common_rows), "common genes and ", length(common_cols), "common locations\n")
    cat(dim(logliklist[[1]])[1], "genes in iteration 1 and ", dim(logliklist[[1]])[2], "common locations\n\n")  
    
    for (i in 1:n){
      tmploglik <- logliklist[[i]]
      tmploglik <- tmploglik[ rownames(tmploglik) %in% common_rows, colnames(tmploglik) %in% common_cols ]
      allLogLikelihood[, i] <- tmploglik
      rm(tmploglik)
      gc()
    }
    
    WAICout <- waic(allLogLikelihood)$estimates["waic", "Estimate"]  #3 seconds-ish for 10 runs
  } else {
    WAICout <- NULL
  }
  
  time.end <- proc.time()
  time_elapsed = time.end[3] - time.start[3]
  out <- list(cell_type_proportion_matrices = cellTypePropMatList,
              rmse = rmseList,
              WAIC = WAICout,
              alpha0 = alpha0,
              alpha1 = alpha1,
              n = n,
              original_seed = seed,
              scenario = scenario,
              time = time_elapsed)
  return(out)
  cat("All done!\n")
}

createCARDObject_HGT <- function(sc_count,sc_meta,spatial_count,spatial_count_OG,spatial_location,ct.varname,ct.select,sample.varname,minCountGene = 100,minCountSpot =5){  

#### QC on scRNASeq dataset
cat(paste0("## QC on scRNASeq dataset! ...\n"))
if(is(sc_count,"matrix")){
	sc_countMat  <- as(as.matrix(sc_count), "sparseMatrix")
	}else if(is(sc_count,"vector")){
		sc_countMat  <- as(t(as.matrix(sc_count)), "sparseMatrix")
		}else if(is(sc_count,"sparseMatrix")){
			sc_countMat  <- sc_count
			}else{
				stop("scRNASeq counts has to be of following forms: vector,matrix or sparseMatrix")
}
if (missing(x = sc_countMat)) {
	stop("Please provide scRNASeq count data")
	} else if (is.null(sample.varname) || missing(sample.varname)) {
		sample.varname = "Sample"
		sc_meta = as.data.frame(sc_meta)
		sc_meta$sampleID = "Sample"
		} else if (any(rownames(x = sc_countMat) == '')) {
			stop("Feature names of sc_count matrix cannot be empty", call. = FALSE)
			} else if(sum(rownames(sc_meta) == colnames(sc_countMat)) != ncol(sc_countMat)){
				stop("Cell name in scRNAseq count data does not match with the rownames of metaData")
				} else if(ncol(sc_countMat)!=nrow(sc_meta)){
					stop("The number of cells in scRNA-seq counts and sc_meta should be consistent! (sc_count -- p x c; sc_meta -- c x 2)")
}
if (is.null(ct.varname)){
	stop("Please provide the column name indicating the cell type information in the meta data of scRNA-seq")
	}else if (is.null(ct.select)){
		cat(paste0("No cell types selected, we will use all the cell types in the scRNA-seq data\n"))
		ct.select <- unique(sc_meta[, ct.varname])
	}
ct.select <- as.character(ct.select[!is.na(ct.select)])
sc_eset = sc_QC(sc_countMat,sc_meta,ct.varname,ct.select,sample.varname)
#### Check the spatial count dataset
#### QC on spatial dataset
cat(paste0("## QC on spatially-resolved dataset! ...\n"))

if(is(spatial_count,"matrix")){
	spatial_countMat  <- as(as.matrix(spatial_count), "sparseMatrix")
	}else if(is(spatial_count,"vector")){
		spatial_countMat  <- as(t(as.matrix(spatial_count)), "sparseMatrix")
		}else if(is(spatial_count,"sparseMatrix")){
			spatial_countMat <- spatial_count
			}else{
				stop("spatial resolved transcriptomic counts has to be of following forms: vector,matrix or sparseMatrix")
}
if (any(rownames(x = spatial_countMat) == '')) {
			stop("Gene names of spatial count matrix cannot be empty", call. = FALSE)
}
commonGene = intersect(rownames(spatial_countMat),rownames(assays(sc_eset)$counts))
if (length(commonGene) == 0) {
			stop("There are no common gene names in spatial count data and single cell RNAseq count data", call. = FALSE)
}
if(is.null(spatial_location)){
	stop("Please provide the matched spatial location data frame")
}
if(ncol(spatial_countMat)!=nrow(spatial_location)){
	stop("The number of spatial locations in spatial_count and spatial_location should be consistent! (spatial_count -- p x n; spatial_location -- n x 2)")
	}# end fi
## check data order should consistent
if(!identical(colnames(spatial_countMat), rownames(spatial_location))){
	stop("The column names of spatial_count and row names of spatial_location should be should be matched each other! (spatial_count -- p x n; spatial_location -- n x 2)")
	}# end fi
#### QC on spatial dataset

########################################################################################################################################################################
# Remove genes/locations in untransformed spatial count matrix and then remove the corresponding rows and columns in spatial_countMat, which is the transformed counts #
########################################################################################################################################################################
spatial_count_OG = spatial_count_OG[rowSums(spatial_count_OG > 0) > minCountSpot,]
spatial_count_OG = spatial_count_OG[,(colSums(spatial_count_OG) >= minCountGene & colSums(spatial_count_OG) <= 1e6)]
spatial_countMat <- spatial_countMat[rownames(spatial_countMat) %in% rownames(spatial_count_OG), ]
spatial_countMat <- spatial_countMat[, colnames(spatial_countMat) %in% colnames(spatial_count_OG)]
#spatial_countMat = spatial_countMat[rowSums(spatial_countMat > 0) > minCountSpot,]
#spatial_countMat = spatial_countMat[,(colSums(spatial_countMat) >= minCountGene & colSums(spatial_countMat) <= 1e6)]
spatial_location = spatial_location[rownames(spatial_location) %in% colnames(spatial_countMat),]
spatial_location = spatial_location[match(colnames(spatial_countMat),rownames(spatial_location)),]

object <- new(
		Class = "CARD",
		sc_eset = sc_eset,
		spatial_countMat = spatial_countMat,
		spatial_location = spatial_location,
		project = "Deconvolution",
		info_parameters = list(ct.varname = ct.varname,ct.select = ct.select,sample.varname = sample.varname)
		)
return(object)
}

#############################################
## Run CARD with a log(1+x) transformation ##
#############################################

SimOverCARD_OSCC_DT <- function(sc_count,  
                               sc_meta,  
                               spatial_count,
                               spatial_location,
                               ct.varname,
                               ct.select, 
                               sample.varname,
                               truePropMat,
                               scenario = 1){
  time.start <- proc.time()
  #things to spit out: original spatial gene expression matrix, original cell-type proportion matrix, list of all H matrices, list of cell-type proportion matrices, true cell type proportion matrix, average cell-type proportions across runs, RMSE for each of 50 tries, RMSE for average cell-type proportion, RMSE for original CARD, alpha, number of runs, original seed 
  #run CARD originally
  CARDobjOG <- createCARDObject(sc_count, sc_meta, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene = 100, minCountSpot = 5)
  CARDobjOG <- CARD_deconvolution(CARDobjOG)
  
  rmseOG <- rmse(CARDobjOG@Proportion_CARD, truePropMat, scenario = scenario)
  
  #Transform the data with the log(1+x) transformation
  Tspatial_count <- log(1.05+spatial_count)
  
  #Apply CARD
  CARDobjT <- createCARDObject(sc_count, sc_meta, Tspatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene = 100, minCountSpot = 5)  #create CARD object with transformed spatial counts matrix
  CARDobjT <- CARD_deconvolution(CARDobjT)  #run CARD with transformed spatial counts
  rmseT <- rmse(CARDobjT@Proportion_CARD, truePropMat, scenario = scenario)
  
  time.end <- proc.time()
  time_elapsed = time.end[3] - time.start[3]
  out <- list(CARD_prop = CARDobjOG@Proportion_CARD,
              CARD_prop_T = CARDobjT@Proportion_CARD,
              rmse_CARD = rmseOG,
              rmse_T = rmseT,
              scenario = scenario,
              time = time_elapsed)
  return(out)
  cat("All done!\n")
}

#####################################
## Run CARD with imputation (ALRA) ##
#####################################

SimOverCARD_OSCC_ALRA <- function(sc_count,  
                                  sc_meta,  
                                  spatial_count,
                                  spatial_location,
                                  ct.varname,
                                  ct.select, 
                                  sample.varname,
                                  truePropMat,
                                  scenario = 1,
				  seed = 1){
  time.start <- proc.time()
  #things to spit out: original spatial gene expression matrix, original cell-type proportion matrix, list of all H matrices, list of cell-type proportion matrices, true cell type proportion matrix, average cell-type proportions across runs, RMSE for each of 50 tries, RMSE for average cell-type proportion, RMSE for original CARD, alpha, number of runs, original seed 
  #run CARD originally
  CARDobjOG <- createCARDObject(sc_count, sc_meta, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene = 100, minCountSpot = 5)
  CARDobjOG <- CARD_deconvolution(CARDobjOG)
  
  rmseOG <- rmse(CARDobjOG@Proportion_CARD, truePropMat, scenario = scenario)
  
  #Impute the data with ALRA
  set.seed(seed)
  spatial_count_alra <- alra(t(spatial_count))
  Tspatial_count <- t(spatial_count_alra[[1]])
  colnames(Tspatial_count) <- colnames(spatial_count)
  
  #Apply CARD
  CARDobjT <- createCARDObject(sc_count, sc_meta, Tspatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene = 100, minCountSpot = 5)  #create CARD object with transformed spatial counts matrix
  CARDobjT <- CARD_deconvolution(CARDobjT)  #run CARD with transformed spatial counts
  rmseT <- rmse(CARDobjT@Proportion_CARD, truePropMat, scenario = scenario)
  
  time.end <- proc.time()
  time_elapsed = time.end[3] - time.start[3]
  out <- list(CARD_prop = CARDobjOG@Proportion_CARD,
              CARD_prop_T = CARDobjT@Proportion_CARD,
              rmse_CARD = rmseOG,
              rmse_T = rmseT,
              scenario = scenario,
              time = time_elapsed)
  return(out)
  cat("All done!\n")
}

####################################
## Run CARD with denoising (MIST) ##
####################################

SimOverCARD_OSCC_MIST <- function(sc_count,  
                                  sc_meta,  
                                  spatial_count,
                                  spatial_location,
                                  ct.varname,
                                  ct.select, 
                                  sample.varname,
                                  truePropMat,
                                  scenario = 1,
                                  seed = 1,
                                  file_name_for_MIST,
                                  file_name_from_MIST){
  time.start <- proc.time()
  #things to spit out: original spatial gene expression matrix, original cell-type proportion matrix, list of all H matrices, list of cell-type proportion matrices, true cell type proportion matrix, average cell-type proportions across runs, RMSE for each of 50 tries, RMSE for average cell-type proportion, RMSE for original CARD, alpha, number of runs, original seed 
  #run CARD originally
  CARDobjOG <- createCARDObject(sc_count, sc_meta, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene = 100, minCountSpot = 5)
  CARDobjOG <- CARD_deconvolution(CARDobjOG)
  
  rmseOG <- rmse(CARDobjOG@Proportion_CARD, truePropMat, scenario = scenario)
  
  # Denoise the data with MIST
  # 1) Save the original spatial counts matrix as a .csv with column names as genes, row names as *xloc*x*yloc*
  spatial_count_transpose <- t(spatial_count)
  write.csv(spatial_count_transpose, file = file_name_for_MIST)
  
  # 2) Run MIST
  input_csv <- file_name_for_MIST
  output_csv <- file_name_from_MIST
 # system("bash -c 'source activate ReST'")
  system(paste("python run_MIST.py", input_csv, output_csv))
  
  # 3) Load in the denoised data
  spatial_count_mist_transpose <- read.csv(file_name_from_MIST, row.names = 1)
 
  Tspatial_count <- t(spatial_count_mist_transpose)
  colnames(Tspatial_count) <- colnames(spatial_count)
  
  #Apply CARD
  CARDobjT <- createCARDObject_HGT(sc_count = sc_count,
                                   sc_meta = sc_meta,
                                   spatial_count = Tspatial_count,
                                   spatial_count_OG = spatial_count,
                                   spatial_location = spatial_location,
                                   ct.varname = ct.varname,
                                   ct.select = ct.select,
                                   sample.varname = sample.varname,
                                   minCountGene = 100,
                                   minCountSpot = 5)  #create CARD object with transformed spatial counts matrix
  CARDobjT <- CARD_deconvolution(CARDobjT)  #run CARD with transformed spatial counts
  rmseT <- rmse(CARDobjT@Proportion_CARD, truePropMat, scenario = scenario)
  
  time.end <- proc.time()
  time_elapsed = time.end[3] - time.start[3]
  out <- list(CARD_prop = CARDobjOG@Proportion_CARD,
              CARD_prop_T = CARDobjT@Proportion_CARD,
              rmse_CARD = rmseOG,
              rmse_T = rmseT,
              scenario = scenario,
              time = time_elapsed)
  return(out)
  cat("All done!\n")
}



transformX2H_noZI <- function(X, alpha1, seed = 123456){
  set.seed(seed)
  n.zeros <- sum(colSums(X == 0))
  n.nonzeros <- sum(colSums(X > 0))
  H <- X  #leave any nonzero entries the same
  H[H > 0] <- (rgamma(n.nonzeros, shape = alpha1 + H[H > 0], rate = alpha1 + 1))
  # H[H == 0] <- logit(rbeta(n.zeros, (alpha0*0.75), (0.25*alpha0+1))) #don't mess with the zeros for this one
  return(H)
}

SimOverCARD_OSCC_noZI <- function(sc_count,  
                             sc_meta,  
                             spatial_count,
                             spatial_location,
                             ct.varname,
                             ct.select, 
                             sample.varname,
                             minCountGene,
                             minCountSpot,
                             truePropMat,
                             scenario = 1,
                             alpha1,
                             n = 100,
                             doWAIC = TRUE,
                             seed = 1){
  time.start <- proc.time()
  #things to spit out: original spatial gene expression matrix, original cell-type proportion matrix, list of all H matrices, list of cell-type proportion matrices, true cell type proportion matrix, average cell-type proportions across runs, RMSE for each of 50 tries, RMSE for average cell-type proportion, RMSE for original CARD, alpha, number of runs, original seed 
  #run CARD originally
  CARDobjOG <- createCARDObject(sc_count, sc_meta, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene = 100, minCountSpot = 5)
  CARDobjOG <- CARD_deconvolution(CARDobjOG)
  
  #and record the original spatial gene expression matrix, the true and CARD cell type proportion matrices, and the CARD RMSE
  #spatialCountList <- list(original_spatial_count = spatial_count), can't afford to save
  rmseOG <- rmse(CARDobjOG@Proportion_CARD, truePropMat, scenario = scenario)
  #rmseList <- list(CARD_rmse = rmseOG)
  cat("done with original CARD\n")
  
  #Initialize the list of transformations and RMSE for each run of MCMC
  TcellTypePropMatList <- list()
  TrmseList <- list()
    
  logliklist <- list()
  
  for(i in 1:n){
    Tspatial_count <- transformX2H_noZI(spatial_count, alpha1 = alpha1, seed = (seed + (i-1)*5000000))  #transform spatial counts matrix
    
    CARDobjtemp <- createCARDObject_HGT(sc_count, sc_meta, Tspatial_count, spatial_count, spatial_location, ct.varname, ct.select, sample.varname, minCountGene, minCountSpot)  #create CARD object with transformed spatial counts matrix
    CARDobjtemp <- CARD_deconvolution_HGT(CARDobjtemp)  #run CARD with transformed spatial counts
    
    TcellTypePropMatList[[i]] <- CARDobjtemp@Proportion_CARD
    Vtemp <- CARDobjtemp@algorithm_matrix$Res$V
    # We can sometimes get very small negative numbers instead of 0s in the cell-type proportions, so we'll need to zero those out and rescale the rows
    TcellTypePropMatList[[i]][TcellTypePropMatList[[i]] < 0] <- 0
    TcellTypePropMatList[[i]] <- sweep(TcellTypePropMatList[[i]], 1, rowSums(TcellTypePropMatList[[i]]), FUN = "/")
    Vtemp[Vtemp < 0] <- 0
    
    rmsetemp <- rmse(TcellTypePropMatList[[i]], truePropMat, scenario = scenario)  #find RMSE for this run
    TrmseList[[i]] <- rmsetemp
    
    # We want to record the WAIC to help choose hyperparameters.  WAIC = -2 (lpd - p_WAIC), and we can calculate this in R with the WAIC function, which requires the loglikelihood in an N x S matrix, where N is total data points and S the number of MCMC samples
    if (doWAIC){
      
      # To find the WAIC, we'll need 1) the original data (restricted to cell-type informative genes and genes/locations that pass QC)
      Xinput_norm <- CARDobjtemp@algorithm_matrix$Xinput_norm
      spatial_count_restricted <- spatial_count[rownames(spatial_count) %in% rownames(Xinput_norm), colnames(spatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
      spatial_count_restricted <- spatial_count_restricted[ match(rownames(Xinput_norm), rownames(spatial_count_restricted)), match(colnames(Xinput_norm), colnames(spatial_count_restricted))]
      
      # 2) the transformed data under the same restrictions
      Tspatial_count_restricted <- Tspatial_count[rownames(Tspatial_count) %in% rownames(Xinput_norm), colnames(Tspatial_count) %in% colnames(Xinput_norm)] #get the genes and locations that pass QC
      Tspatial_count_restricted <- Tspatial_count_restricted[ match(rownames(Xinput_norm), rownames(Tspatial_count_restricted)), match(colnames(Xinput_norm), colnames(Tspatial_count_restricted))]
      
      # 3) the things we need from the CARD model, B, V, and sigma_e^2
      B <- CARDobjtemp@algorithm_matrix$B
      V <- Vtemp
      sigma_e2 <- CARDobjtemp@algorithm_matrix$Res$sigma_e2
      
      # 4) the latent mean of the CARD model, X^(2) = BV'
      X2 <- B %*% t(V)
      
      # 5) the log-likelihood at each data point
      #logliklist <- list()
      logliktemp <- matrix(0, dim(spatial_count_restricted)[1], dim(spatial_count_restricted)[2])
      Idx0 <- which(spatial_count_restricted == 0)
      Idxg0 <- which(spatial_count_restricted > 0)
      rownames(logliktemp) <- rownames(spatial_count_restricted)
      colnames(logliktemp) <- colnames(spatial_count_restricted)
      
      #for the zero data, we have a point mass at 0 * point mass at X^(B)_ij * dbern(0, p = invlogit( X^(B)_ij ) * dnorm( X^(1)_ij, mean = X^(2)_ij, var = sigma_e2) and X^(B) = X^(1)
      logliktemp[Idx0] <- dbern(0, prob = invlogit(Tspatial_count_restricted[Idx0]), log=TRUE) + dnorm(Tspatial_count_restricted[Idx0], mean = X2[Idx0], sd = sqrt(sigma_e2), log=TRUE)
      logliktemp[Idxg0] <- dpois(spatial_count_restricted[Idxg0], lambda = Tspatial_count_restricted[Idxg0], log=TRUE) + dnorm(Tspatial_count_restricted[Idxg0], mean = X2[Idxg0], sd = sqrt(sigma_e2), log=TRUE)
      logliklist[[i]] <- logliktemp
      
      rm(logliktemp)
    }
    
    rm(CARDobjtemp)
    rm(Tspatial_count)
    rm(rmsetemp)
    gc()
    cat("done with transformed run", i, "\n")
  }
  
 
  # Now work with average cell-type proportions from our different runs
  cellTypePropMatAve <- Reduce("+", TcellTypePropMatList) / n
  rmseAve <- rmse(cellTypePropMatAve, truePropMat, scenario = scenario)
  
  medianArray <- array(unlist(TcellTypePropMatList), c(dim(TcellTypePropMatList[[1]]), length(TcellTypePropMatList)))
  medianMat <- apply(medianArray, 1:2, median)
  cellTypePropMatMedian <- sweep(medianMat, 1, rowSums(medianMat), FUN = "/")
  rownames(cellTypePropMatMedian) <- rownames(cellTypePropMatAve)
  colnames(cellTypePropMatMedian) <- colnames(cellTypePropMatAve)
  rmseMedian <- rmse(cellTypePropMatMedian, truePropMat, scenario = scenario)
  
  cellTypePropMatList <- list(true_ct_prop_mat = truePropMat, CARD_ct_prop_mat = CARDobjOG@Proportion_CARD, mean_ct_prop = cellTypePropMatAve, median_ct_prop = cellTypePropMatMedian, transformed_ct_prop_mat_list = TcellTypePropMatList)
  rmseList <- list(CARD_rmse = rmseOG, Ave_rmse = rmseAve, Median_rmse = rmseMedian, transformed_rmse_list = TrmseList)

  #Get the WAIC
  if (doWAIC){
    
    # We can only look at locations and genes that each iteration has in common (luckily it's rare that any iteration loses locations or genes)
    common_rows <- rownames(logliklist[[1]])
    common_cols <- colnames(logliklist[[1]])
    for (i in 1:n){
      common_rows <- intersect(common_rows, rownames(logliklist[[i]]))
      common_cols <- intersect(common_cols, colnames(logliklist[[i]]))
    }
    
    totaldatapoints <- length(common_rows)*length(common_cols)
    allLogLikelihood <- matrix(0, nrow = totaldatapoints, ncol = n)
    cat(length(common_rows), "common genes and ", length(common_cols), "common locations\n")
    cat(dim(logliklist[[1]])[1], "genes in iteration 1 and ", dim(logliklist[[1]])[2], "common locations\n\n")  
    
    for (i in 1:n){
      tmploglik <- logliklist[[i]]
      tmploglik <- tmploglik[ rownames(tmploglik) %in% common_rows, colnames(tmploglik) %in% common_cols ]
      allLogLikelihood[, i] <- tmploglik
      rm(tmploglik)
      gc()
    }
    
    WAICout <- waic(allLogLikelihood)$estimates["waic", "Estimate"]  #3 seconds-ish for 10 runs
  } else {
    WAICout <- NULL
  }
  
  time.end <- proc.time()
  time_elapsed = time.end[3] - time.start[3]
  out <- list(cell_type_proportion_matrices = cellTypePropMatList,
              rmse = rmseList,
              WAIC = WAICout,
              alpha1 = alpha1,
              n = n,
              original_seed = seed,
              scenario = scenario,
              time = time_elapsed)
  return(out)
  cat("All done!\n")
}