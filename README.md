ZI-HGT Reproducibility
================

This GitHub repository contains code to reproduce all figures, real data results, and simulation results from the paper 
_A Zero-Inflated Hierarchical Generalized Transformation Model to Address Non-Normality in Spatially-Informed Cell-Type Deconvolution_
by Hunter J. Melton, Jonathan R. Bradley, and Chong Wu.

## Organization

There are four directories that contain reproduction code, plus one that contains the figures.  Each code directory
contains a README file that briefly describes the files within it, and each code file includes instructions for running
the code.

## Figures

The figures from the manuscript are contained within this directory.

## Figures_Code

Code to generate the figures is contained in this directory.  Each individual file has instructions and code for prepping
the results for the plot and for making the plot.  Oftentimes the raw results are too big to store on GitHub, but the 
prepared results will be available.  Any raw results not found here may be generated following the code in the 
Real_Data_Analysis directory or Simulations directory, or they may be found on [OSF](https://osf.io/kygsx/).

## Real_Data_Analysis

Code to generate the results based on the real data, 12 OSCC TME samples collected by Arora et al.  We ran this code on
the university's high performance computer, which uses a Slurm workload manager.  For exact specifications of the job we
ran, see the .sh file.  The real data is unfortunately too  large to store on GitHub, but it may be found on [OSF](https://osf.io/kygsx/).

## Simulations

## Utilities




