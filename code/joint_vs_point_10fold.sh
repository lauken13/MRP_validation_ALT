#!/bin/sh

#SBATCH --job-name=jointvpoint
#SBATCH --time=00:30:00
#SBATCH --mail-user=lauren.kennedy1@monash.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4
#SBATCH --array=1-100

module load R

R CMD BATCH --no-save --no-restore joint_vs_point_10fold.R script_$SLURM_ARRAY_TASK_ID