# Using the ZI-HGT

## Installation

Until the ZI-HGT becomes a package (at which point this README will look a bit different), please download the HGTfunctions3.R script and source it, as shown in the tutorial.

## If you're interested in using the ZI-HGT, we anticipate that it is for one of two reasons.

### You want to use the ZI-HGT + CARD to conduct more accurate cell-type deconvolution with built in uncertainty quantification.

In this case, please follow along with the [tutorial](https://hjm2217.github.io/ZI-HGT/Tutorial.html).

You may need to do some preprocessing of your reference scRNAseq data however you typically do (we use Seurat, and you can see an example in _Replication/Real
_Data_Analysis/ProcessReference.R_). Be sure you've determined the cell types in your reference scRNAseq data and added that information to the meta data.

Once your reference scRNAseq data is prepared, you can run the analysis as in _Tutorial.html_.  You will need:
* an object containing the spatial transcriptomics counts
* an object containing the spatial locations
* an object containing the scRNAseq counts
* an object containing the scRNAseq meta data
* to source the functions in _HGTfunctions3.R_

If you need help or run into any issues, please let us know!

### You want to use the ZI-HGT as a wrapper function for your another method.

This will be a little bit trickier, as you will need to adapt and edit our code for the ZI-HGT to work with your method.  On the bright side, assuming your method works on the spatial transcriptomics counts, and you would like a more normal version of those counts and built in UQ, it should not be too difficult.

Essentially, you will need to adapt the structure of the code found in the _HGTplusCARD()_ function in _HGTfunctions3.R_to work with your method.  The steps you need to accomplish are:
1) Get all of the data, reference data, locations, and anything else you may need to run your analysis together.
2) Following the _HGTplusCARD()_ function, set up a loop from 1 to n (we recommend starting with n = 100).  Within each iteration, transform the spatial count data with the ZI-HGT using the hyperparameters and a reproducible seed.  Apply your method and record the results.
3) If you wish to calculate the WAIC to decide between choices of hyperparameters _alpha0_ and _alpha1_, you'll have to determine the log likelihood (as we did in the if (doWAIC) {...} section of the the _HGTplusCARD()_ function).
4) After the loop, collect summaries of the results - the mean, median, variance, 95% BCIs, etc. as we did.  If you wish to calculate the WAIC, you'll do that here too.
5) If you've decided to use the WAIC, after you've run all of the analysis across all hyperparameters and samples, choose your hyperparameters for each sample by minimizing the WAIC.

Those are the general steps to implementing the ZI-HGT agnostically around an unknown method.  If you run into any trouble, or you'd like to discuss how to use the ZI-HGT with your method, please reach out and let us know!  We are happy to help.
