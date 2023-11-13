#!/bin/sh

#SBATCH -p skylake
#SBATCH --job-name=sae_scores
#SBATCH --time=06:00:00
#SBATCH --mail-user=lauren.a.kennedy@adelaide.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100

module use /apps/skl/modules/all/
module load R/4.2.3-foss-2021b

R CMD BATCH --no-save --no-restore score_model_sae.R sae_score_$SLURM_ARRAY_TASK_ID