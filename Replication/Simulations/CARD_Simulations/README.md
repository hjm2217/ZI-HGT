# CARD Simulations

This directory contains all of the code needed to generate the simulation results following the exact
simulation setup in [CARD](https://www.nature.com/articles/s41587-022-01273-7).  
Unfortunately, the simulated data and results are too large to comfortably share.  

## Simulating the ST Data

### Simulate_CARD_Data.R

This file (along with the accompanying submission script Simulate_CARD_Data.sh) generates
the simulated ST data for the comparison of the ZI-HGT + CARD vs. CARD alone on ST datasets with different
noisiness levels and on different simulation scenarios following the CARD simulations.  The results from the
accompanying analysis of this simulated data are displayed in Figure 4 and Supplementary Figures 4-6.  
Exact specifications of the job we ran can be found in the two files.

This relies on the scRNAseq data found in the _Reference_ subdirectory and the mouse olfactory bulb (MOB)
reference information also found in the _Reference_ subdirectory.  The simulated ST data are saved in the _Simulated_ST_Data/_ directory, which is empty here.

## Analyzing the Simulated ST Data

### runHGTCARD_CARD_Sc1-4.R

This file (along with the accompanying submission script runHGTCARD_CARD_Sc1-4.sh) runs the ZI-HGT + CARD and CARD alone on the simulated MOB ST data following the CARD simulations for all noisiness levels (Mixnoise0-3) and the first four simulation scenarios.  The results from this analysis are shown in Figure 4 and Supplementary Figures 4-6. Exact specifications of the job can be found in the two files. We note that, as currently constructed, you will run 2 simulations each across 1800 arrays on your HPC. On our university's HPC, we could only run 100 arrays at any given time, so we would adjust the submission script's array argument to read --array=1-100. After running this, we would run it again with array=101-200, and so on.

The results from analyzing these simulations are saved in the Results_CARD/ directory, which is empty here.


### runHGTCARD_CARD_Sc5.R

This file (along with the accompanying submission script runHGTCARD_CARD_Sc5.sh) runs the ZI-HGT + CARD and CARD alone on the simulated MOB ST data following the CARD simulations for all noisiness levels (Mixnoise0-3) and simulation scenario 5 (the mismatched single-cell reference).  The results from this analysis are shown in Figure 4 and Supplementary Figures 4-6. Exact specifications of the job can be found in the two files. We note that, as currently constructed, you will run 3 simulations each across 1200 arrays on your HPC. On our university's HPC, we could only run 100 arrays at any given time, so we would adjust the submission script's array argument to read --array=1-100. After running this, we would run it again with array=101-200, and so on.

The results from analyzing these simulations are saved in the Results_CARD/ directory, which is empty here.


### runCARD.R

This file (along with the accompanying submission script runCARD.sh) runs CARD  on the simulated MOB ST data following the CARD simulations for all noisiness levels (Mixnoise0-3), all simulation scenarios (1-5), and all different sequencing depths (which we control by controlling the number of single cells at each spatial location). The results from this analysis are shown in Figure 4 and Supplementary Figures 4-6. Exact specifications of the job can be found in the two files. 

The results from analyzing these simulations are saved in the Results_CARD/ directory, which is empty here.

