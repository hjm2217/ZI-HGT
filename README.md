The ZI-HGT
================

This GitHub repository contains instructions to use the ZI-HGT in your work and code to reproduce all figures, 
real data results, and simulation results from the paper _A Zero-Inflated Hierarchical Generalized Transformation Model to Address Non-Normality in 
Spatially-Informed Cell-Type Deconvolution_ by Hunter J. Melton, Jonathan R. Bradley, and Chong Wu.

## Using the ZI-HGT

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

## Instructions for Replication and Organization of the files

There are four directories that contain reproduction code, plus one that contains the figures.  Each code directory
contains a README file that briefly describes the files within it, and each code file includes instructions for running
the code.

To fully replicate the results shown in the paper, follow the steps given below.

1) Download the data.  You can find the OSCC ST data from Arora et al. [here](https://figshare.com/articles/dataset/Spatial_transcriptomics_reveals_distinct_and_conserved_tumor_core_and_edge_architectures_that_predict_survival_and_targeted_therapy_response_/20304456/1) and the corresponding scRNAseq data from Puram et al. [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322).

2) Go to the _Real_Data_Analysis_ directory.  Process the scRNAseq data with _ProcessReference.R_ (and using the _Data/cluster_identification.csv_ file.  Then run the real data analysis by running _Real_Data_Analysis.R_.  On the logon node of your local HPC, run sbatch _Real_Data_Analysis.sh_ (you'll need to make a few edits to the .sh file with your desired job name, partition name, etc.  Once you've generated all of the results and they're saved to the _Results_ subdirectory, you're in good shape here.

3) You can now plot all of our real data results.  Go to the _Figures_Code_ directory.  You can run the code in _Figure2.R_, _Figure3.R_, and _Figure5_S7-S41.R_ to generate the corresponding figures.

4) Go to the _Simulations_ directory.  You can follow along running the scripts as they're described in the _Simulations/README.md_ file.

5) You can now plot all of the simulations results.  Go to the _Figures_Code_ directory.  You can run the code in _Figure4_S4-S6.R_, _FigureS1.R_, _FigureS2.R_, and _FigureS3.R_ to generate the corresponding figures.

## Figures

The figures from the manuscript are contained within this directory.

## Figures_Code

Code to generate the figures is contained in this directory.  Each individual file has instructions and code for prepping
the results for the plot and for making the plot.  Oftentimes the results are too big to store on GitHub, but the 
prepared results will be available.  Any real data results not found here may be generated following the code in the 
Real_Data_Analysis directory or Simulations directory, or they may be found on [OSF](https://osf.io/kygsx/).

## Real_Data_Analysis

Code to generate the results based on the real data, 12 OSCC TME samples collected by Arora et al.  We ran this code on
the university's high performance computer, which uses a Slurm workload manager.  For exact specifications of the job we
ran, see the .sh file.  The real data is unfortunately too large to store on GitHub, but it may be found on [OSF](https://osf.io/kygsx/).

## Simulations

Code to generate and analyze all simulationed data included in the project.  Most simulations are built to mimic the structure of the OSCC TME ST data, though we followed CARD's exact simulations as well for additional analysis.  We ran all code on the university's high performance computer, which uses a Slurm workload manager.  For exact specifications of the jobs we ran, see the .sh file.

## Utilities

All of the functions used to preprocess and analyze data, build and analyze simulations, and in general get everything done for the project, are in one of the two files in this directory.



