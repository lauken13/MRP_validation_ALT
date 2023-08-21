#!/bin/sh

#SBATCH -p skylake
#SBATCH --job-name=mse_mrp_3_2
#SBATCH --time=04:00:00
#SBATCH --mail-user=lauren.a.kennedy@adelaide.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100


module use /apps/skl/modules/all/
module load R/4.2.3-foss-2021b

R CMD BATCH --no-save --no-restore section_crps_fullyobs.R script_$SLURM_ARRAY_TASK_ID