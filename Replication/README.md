## Instructions for Replication and Organization of the files

There are four directories that contain reproduction code, plus one that contains the figures.  Each code directory
contains a README file that briefly describes the files within it, and each code file includes instructions for running
the code.

To fully replicate the results shown in the paper, follow the steps given below.

1) Download the data.  You can find the OSCC ST data from Arora et al. [here](https://figshare.com/articles/dataset/Spatial_transcriptomics_reveals_distinct_and_conserved_tumor_core_and_edge_architectures_that_predict_survival_and_targeted_therapy_response_/20304456/1) and the corresponding scRNAseq data from Puram et al. [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322).

2) Go to the _Real_Data_Analysis_ directory.  Process the scRNAseq data with _ProcessReference.R_ (and using the _Data/cluster_identification.csv_ file.  Then run the real data analysis by running _Real_Data_Analysis.R_.  On the logon node of your local HPC, run sbatch _Real_Data_Analysis.sh_ (you'll need to make a few edits to the .sh file with your desired job name, partition name, etc.  Once you've generated all of the results and they're saved to the _Results_ subdirectory, you're in good shape here.

3) You can now plot all of our real data results.  Go to the _Figures_Code_ directory.  You can run the code in _Figure2.R_, _Figure3.R_, and _Figure5-6_S8-S41.R_, _Figure_S42.R_, _Figure_S43.R_, _Figure_S44.R_, and _Figure_S45.R_ to generate the corresponding figures.

4) Go to the _Simulations_ directory.  You can follow along running the scripts as they're described in the _Simulations/README.md_ file.

5) You can now plot all of the simulations results.  Go to the _Figures_Code_ directory.  You can run the code in _Figure4_S4-S7.R_, _FigureS1.R_, _FigureS2.R_, and _FigureS3.R_ to generate the corresponding figures.

## Figures

The figures from the manuscript are contained within this directory.

## Figures_Code

Code to generate the figures is contained in this directory.  Each individual file has instructions and code for summarizing
the results for the plot and for making the plot.  Oftentimes the results are too big to store on GitHub, but the 
summarized results will be available.  Any real data results not found here may be generated following the code in the 
Real_Data_Analysis directory or Simulations directory, or they may be found on [OSF](https://osf.io/kygsx/).

## Real_Data_Analysis

Code to generate the results based on the real data, 12 OSCC TME samples collected by Arora et al.  We ran this code on
the university's high performance computer, which uses a Slurm workload manager.  For exact specifications of the job we
ran, see the .sh file.  The real data is unfortunately too large to store on GitHub, but it may be found on [OSF](https://osf.io/kygsx/).

## Simulations

Code to generate and analyze all simulationed data included in the project.  Most simulations are built to mimic the structure of the OSCC TME ST data, though we followed CARD's exact simulations as well for additional analysis.  We ran all code on the university's high performance computer, which uses a Slurm workload manager.  For exact specifications of the jobs we ran, see the .sh file.

## Utilities

All of the functions used to preprocess and analyze data, build and analyze simulations, etc. are in one of the two files in this directory.
