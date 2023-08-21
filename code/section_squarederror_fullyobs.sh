#!/bin/sh

#SBATCH -p skylake
#SBATCH --job-name=mse_mrp_3_1
#SBATCH --time=04:00:00
#SBATCH --mail-user=lauren.a.kennedy@adelaide.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=32000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100


module use /apps/skl/modules/all/
module load R/4.1.2

R CMD BATCH --no-save --no-restore section_squarederror_fullyobs.R script_$SLURM_ARRAY_TASK_ID