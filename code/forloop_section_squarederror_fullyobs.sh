#!/bin/sh

#SBATCH -p skylake
#SBATCH --job-name=mse_mrp
#SBATCH --time=35:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB

module use /apps/skl/modules/all/
module load R/4.1.2

R CMD BATCH --no-save --no-restore forloop_section_squarederror_fullyobs.R script_$SLURM_ARRAY_TASK_ID