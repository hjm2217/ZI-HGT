If you're interested in using the ZI-HGT, we anticipate that it is for one of two reasons.

### You want to use the ZI-HGT + CARD to conduct more accurate cell-type deconvolution with built in uncertainty quantification.

In this case, you should navigate to the _Real_Data_Analysis_ directory.  You'll want to preprocess your reference scRNAseq data however you typically do (we use Seurat, and you can see an example in _ProcessReference.R_).  Be sure you've determined the cell types in your referene scRNAseq data and added that information to the meta data.

Once your reference scRNAseq data is prepared, you can run the analysis as in _Real_Data_Analysis.R_.  You will need
* an object containing the scRNAseq counts
* an object containing the scRNAseq meta data
* to source the functions in /Utilities/HGTfunctions3.R
* Assuming you are running the ZI-HGT + CARD in parallel across multiple arrays on an HPC, to determine the array id and job.id as in _Real_Data_Analysis.R_.
* to build a dataframe similar to our _grid_ that contains the hyperparameters, sample (called dataset in _Real_Data_Analysis.R_ - you can just set grid$dataset <- 1 if you only have one sample), and seed.

Once you have these set up, you can run the ZI-HGT + CARD.  Within your loop over the job.id, be sure to load in the correct sample and determine the spatial_count and spatial_location dataframes.  Then plug these, along with the scRNAseq reference data and metadata, the cell type variable name, the vector of cell types, etc. into the _HGTplusCARD()_ function as in _Real_Data_Analysis.R_.  We recommend using the same values for the hyperparameters $\alpha_0$ and $\alpha_1$ that we did, letting minCountGene and minCountSpot be 100 and 5 respectively (as in CARD), letting n, the number of iterations of the ZI-HGT + CARD, be 100, and letting doWAIC = TRUE.  The last of these ensures that the WAIC will be calculated for each set of hyperparameters and the best set can be chosen by minimizing it.

After running the analysis, you'll likely want to build cell-type proportions plots as we did.  For this, we recommend following along with _/Figures_Code/Figure5_S7-S41.R_.  You will select the best set of hyperparameters using the WAIC then construct cell-type proportions plots for all samples and all cell-types.  Simply comment out the lines regarding plots _p2_ and _p1_ (as these construct comparison plots using CARD only).

If you run into any trouble, please reach out and let us know.  We would be happy to help!

### You want to use the ZI-HGT as a wrapper function for your another method.

This will be a little bit trickier, as you will need to adapt and edit our code to work with your method.  On the bright side, assuming your method works on the counts, and you would like a more normal version of those counts and built in UQ, it should not be too difficult.

Essentially, you will need to adapt the structure of the code found in _Real_Data_Analysis.R_ and the _HGTplusCARD()_ function in _/Utilities/HGTfunctions3.R_to work with your method.  The steps you need to accomplish are:
1) Get all of the data, reference data, locations, and anything else you may need to run your analysis together.
2) Set up a dataframe similar to our _grid_ to run the analysis over.  This should have the ZI-HGT hyperparameters (which we would recommend leaving as what we have for a starting point), the sample, and the seed for each run.  We recommend that each row of the dataset correspond to a separate array on your local HPC to minimize computational burden.
3) Following the _HGTplusCARD()_ function, set up a loop from 1 to n (we recommend starting with n = 100).  Within each iteration, transform the spatial count data with the ZI-HGT using the hyperparameters and seed.  Apply your method and record the results (we used a list to do so).
4) If you wish to calculate the WAIC, you'll have to determine the log likelihood (as we did in the if (doWAIC) {...} section of the the _HGTplusCARD()_ function).
5) After the loop, collect summaries of the results - the mean, median, 95% BCIs, etc. as we did.  If you wish to calculate the WAIC, you'll do that here too.
6) If you've decided to use the WAIC, after you've run all of the analysis across all hyperparameters and samples, choose your hyperparameters for each sample by minimizing the WAIC.

Those are the steps to implementing the ZI-HGT agnostically around an unknown method.  If you run into any trouble, please reach out and let us know.  We are happy to help!
