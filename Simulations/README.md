# Simulations

This directory contains all the code needed to generate all of the simulation results presented in the manuscript.
Unfortunately, the simulated data and results are too large to comfortably store on GitHub.  
You may find them on our [OSF](https://osf.io/kygsx/).

We ran this analysis on our university's high performance computer, which uses the Slurm workload manager.
For the exact specifications of each job we ran, see the .sh files.  You may need to adjust some of the
options depending on the computational resources you have access to, and any options you'll need to fill in yourself have 
underscores before and after the option (for example, \_JOB\_NAME\_HERE\_).

We describe each file and subdirectory below, split into three categories: 1) simulating the scRNA-seq data, 2) simulating the ST data, 3) subdirectories and data storage.

# Simulating the scRNA-seq data (All contained in the scRNAseq_Simulations subdirectory)

### SPARSim_OSCC_2.R and SPARSim_OSCC_2_Bigger_Phi_only.R

These files (along with the accompanying submission scripts _SPARSim_OSCC_2.sh_ and _SPARSim_OSCC_2_Bigger_Phi_only.sh_) generate the simulated scRNAseq data used in constructing the simulated
ST data for the comparison of the ZI-HGT + CARD vs. CARD alone on different ST datasets with different
sparsity levels for OSCC Sample 2 (Figure 3).  The simulated scRNAseq datasets may be found in the _scRNAseq_sim_data/_ subdirectory on [OSF](https://osf.io/kygsx/), though not on GitHub as they are too large.
You'll need to run both files to generate all of the necessary simulated scRNAseq data.  Exact specifications
of the jobs we ran can be found in the files.

### Split_SPARSim_OSCC_scRNAseq.R and Split_SPARSim_OSCC_scRNAseq_Bigger_Phi_only.R

These files (along with the accompanying submission scripts _Split_SPARSim_OSCC_scRNAseq.sh_ and _Split_SPARSim_OSCC_scRNAseq_Bigger_Phi_only.sh_) split the aforementioned simulated scRNAseq data into 
two parts.  The first split of the data is used to generate the simulated ST data, the second split is used
as the cell-type deconvolution reference dataset.  The split simulated scRNAseq datasets may be found in the _scRNAseq_sim_data/_ subdirectory on [OSF](https://osf.io/kygsx/), though not on GitHub as they are too large.
You'll need to run both files to split all of simulated scRNAseq data.  Exact specifications
of the jobs we ran can be found in the files.

# Simulating the ST data

### ST_SPARSim_OSCC_Sample_2_Different_Sparsities.R

This file (along with the accompanying submission script ST_SPARSim_OSCC_Sample_2_Different_Sparsities.sh) generates
the simulated ST data for the comparison of the ZI-HGT + CARD vs. CARD alone on ST datasets with different
sparsity levels for OSCC Sample 2.  The results from the accompanying analysis of this simulated data are primarily displayed in Figure 3.  Exact specifications
of the job we ran can be found in the two files.

This relies on the simulated scRNAseq data found in the _scRNAseq_Simulations_ subdirectory and the OSCC Sample 2
reference information found in the _Reference_ subdirectory.  The simulated ST data are saved in the _SPARSIM_ST_Data/_ directory, which is empty here on GitHub but contains all simulated data on [OSF](https://osf.io/kygsx/).

### ST_SPARSim_OSCC_Sample_2_Different_Numbers_Cell_Types.R

This file (along with the accompanying submission script ST_SPARSim_OSCC_Sample_2_Different_Numbers_Cell_Types.sh) generates
the simulated ST data for the comparison of the ZI-HGT + CARD vs. CARD alone with different
numbers of cell types on "realistically sparse" simulated ST data for OSCC Sample 2.  The results from the accompanying analysis of this simulated data are primarily displayed in Supplementary Figure 1.  Exact specifications of the job we ran can be found in the two files.

This relies on the simulated scRNAseq data found in the _scRNAseq_Simulations_ subdirectory and the OSCC Sample 2
reference information found in the _Reference_ subdirectory.  The simulated ST data are saved in the _SPARSIM_ST_Data/_ directory, which is empty here on GitHub but contains all simulated data on [OSF](https://osf.io/kygsx/).

# Analyzing the Simulated ST data

### runHGTCARD_OSCC_SPARSims_Sample_2_Different_Sparsities.R

This file (along with the accompanying submission script _runHGTCARD_OSCC_SPARSims_Sample_2_Different_Sparsities.sh_) runs the ZI-HGT + CARD and CARD alone on the simulated ST data for OSCC Sample 2 with different sparsity levels.  The results from this simulation are shown in Figure 3.  Exact specifications of the job can be found in the two files.  We note that, as currently constructed, you will run 3 simulations each across 1500 arrays on your HPC.  On our university's HPC, we could only run 100 arrays at any given time, so we would adjust the submission script's array argument to read --array=1-100.  After running this, we would run it again with array=101-200, and so on.

The results from these simulations are saved in the _SPARSim_Results/_ directory, which is empty here on GitHub but contains all results on [OSF](https://osf.io/kygsx/).

### runHGTCARD_OSCC_SPARSims_Sample_2_Different_Numbers_Cell_Types.R

This file (along with the accompanying submission script _runHGTCARD_OSCC_SPARSims_Sample_2_Different_Numbers_Cell_Types.sh_) runs the ZI-HGT + CARD and CARD alone with different numbers of cell types on "realistically sparse" simulated ST data for OSCC Sample 2.  The results from this simulation are shown in Supplementary Figure 1.  Exact specifications of the job can be found in the two files.  We note that, as currently constructed, you will run 3 simulations each across 1200 arrays on your HPC.  On our university's HPC, we could only run 100 arrays at any given time, so we would adjust the submission script's array argument to read --array=1-100.  After running this, we would run it again with array=101-200, and so on.

The results from these simulations are saved in the _SPARSim_Results/_ directory, which is empty here on GitHub but contains all results on [OSF](https://osf.io/kygsx/).



