#!/bin/bash                                                                    
#SBATCH -J _JOB_NAME_HERE_
#SBATCH -p _PARTITION_NAME_HERE_
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH -t 4:00:00
#SBATCH --array=1-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=_YOUR_EMAIL_HERE_

#SBATCH -C "YEAR2017|YEAR2018|YEAR2019"                                         

export PATH=$PATH:_PATH_TO_SOFTWARE_HERE_

module load gnu/9.1.1
module load gnu-openmpi

module load gnu
module load R/4.2.0

R CMD BATCH --no-save --no-restore SPARSim_OSCC_2_Bigger_Phi_only.R ./log/SPARSim_OSCC_2_Bigger_Phi_only_$SLURM_ARRAY_TASK_ID.txt

