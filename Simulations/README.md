# Simulations

This directory contains all the code needed to generate all of the simulation results presented in the manuscript.
Unfortunately, the simulated data and results are too large to comfortably store on GitHub.  
You may find them on our [OSF](https://osf.io/kygsx/).

We ran this analysis on our university's high performance computer, which uses the Slurm workload manager.
For the exact specifications of each job we ran, see the .sh files.  You may need to adjust some of the
options depending on the computational resources you have access to, and any options you'll need to fill in yourself have 
underscores before and after the option (for example, \_JOB\_NAME\_HERE\_).

We describe each file and subdirectory below, split into three categories: 1) simulating the scRNA-seq data, 2) simulating the ST data, 3) subdirectories and data storage.

# Simulating the scRNA-seq data

### SPARSim_OSCC_2.R and SPARSim_OSCC_2_Bigger_Phi_only.R

These files (along with the accompanying submission scripts _SPARSim_OSCC_2.sh_ and _SPARSim_OSCC_2_Bigger_Phi_only.sh_) generate the simulated scRNAseq data used in constructing the simulated
ST data for the comparison of the ZI-HGT + CARD vs. CARD alone on different ST datasets with different
sparsity levels for OSCC Sample 2 (Figure 3).  The simulated scRNAseq datasets may be found in the _scRNAseq_sim_data/_ subdirectory on [OSF](https://osf.io/kygsx/), though not on GitHub as they are too large.
You'll need to run both files to generate all of the necessary simulated scRNAseq data.

### Split_SPARSim_OSCC_scRNAseq.R and Split_SPARSim_OSCC_scRNAseq_Bigger_Phi_only.R

These files (along with the accompanying submission scripts _Split_SPARSim_OSCC_scRNAseq.sh_ and _Split_SPARSim_OSCC_scRNAseq_Bigger_Phi_only.sh_) split the aforementioned simulated scRNAseq data into 
two parts.  The first split of the data is used to generate the simulated ST data, the second split is used
as the cell-type deconvolution reference dataset.  The split simulated scRNAseq datasets may be found in the _scRNAseq_sim_data/_ subdirectory on [OSF](https://osf.io/kygsx/), though not on GitHub as they are too large.
You'll need to run both files to split all of simulated scRNAseq data.

# Simulating the ST data

### ST_SPARSim_OSCC_Sample_2_Different_Sparsities.R

This file (along with the accompanying submission script ST_SPARSim_OSCC_Sample_2_Different_Sparsities.sh) generates
the simulated ST data for the comparison of the ZI-HGT + CARD vs. CARD alone on different ST datasets with different
sparsity levels for OSCC Sample 2.  The results from this analysis are displayed in Figure 3.  Exact specifications
of the job we ran can be found in the two files.

This relies on the simulated scRNAseq data found in the _scRNAseq_Simulations_ subdirectory and the OSCC Sample 2
reference information found in the _Reference_ subdirectory.

