#!/bin/sh

#SBATCH -p skylake
#SBATCH --job-name=mse_mrp_3_1
#SBATCH --time=02:00:00
#SBATCH --mail-user=lauren.a.kennedy@adelaide.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=4
#SBATCH --array=1-2


module use /apps/skl/modules/all/
  module load R/4.1.2

R CMD BATCH --no-save --no-restore test4.R script_test_job$SLURM_ARRAY_TASK_ID