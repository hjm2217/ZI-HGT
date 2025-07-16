# Real_Data_Analysis

This directory contains the code used to generate the real data results, which come from 12 OSCC TME ST samples collected 
by Arora et al.  We ran this analysis on our university's high performance computer, which uses the Slurm workload manager.
For the exact specifications of the job we ran, see the Real_Data_Analysis.sh file.  You may need to adjust some of the
options depending on the computational resources you have access to, and any options you'll need to fill in yourself have 
underscores before and after the option (for example, \_JOB\_NAME\_HERE\_).

Unfortunately, the data files are too large to comfortably store and share.  You may find them on our [OSF](https://osf.io/kygsx/).  The results generated using the WAIC-chosen hyperparameters are directly saved on our OSF, and the rest may be found in the Results.tar.gz archive.

Before conducting the cell-type deconvolution, you'll need to process the single-cell data from Puram et al. with ProcessReference.R.
