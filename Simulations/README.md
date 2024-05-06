# Simulations

This directory contains all the code needed to generate all of the simulation results presented in the manuscript.
Unfortunately, the simulated data and results are too large to comfortably store on GitHub.  
You may find them on our [OSF](https://osf.io/kygsx/).

We ran this analysis on our university's high performance computer, which uses the Slurm workload manager.
For the exact specifications of each job we ran, see the .sh files.  You may need to adjust some of the
options depending on the computational resources you have access to, and any options you'll need to fill in yourself have 
underscores before and after the option (for example, \_JOB\_NAME\_HERE\_).

There are quite a few files and subdirectories in this directory, so we will describe each of them below.

# ST_SPARSim_OSCC_Sample_2_Different_Sparsities.R

This file (along with the accompanying submission script ST_SPARSim_OSCC_Sample_2_Different_Sparsities.sh) generates
the simulated ST data for the comparison of the ZI-HGT + CARD vs. CARD alone on different ST datasets with different
sparsity levels for OSCC Sample 2.  The results from this analysis are displayed in Figure 3.  Exact specifications
of the job we ran can be found in the two files.

This relies on the simulated scRNAseq data found in the _scRNAseq_Simulations_ subdirectory and the OSCC Sample 2
reference information found in the _Reference_ subdirectory.

