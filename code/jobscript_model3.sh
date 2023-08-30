#!/bin/sh

#SBATCH -p skylake
#SBATCH --job-name=score_mrp_model3
#SBATCH --time=04:00:00
#SBATCH --mail-user=lauren.a.kennedy@adelaide.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=32000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100


module use /apps/skl/modules/all/
module load R/4.2.3-foss-2021b

R CMD BATCH --no-save --no-restore score_model3.R script_model3_$SLURM_ARRAY_TASK_ID