#!/bin/sh

#SBATCH -p skylake
#SBATCH --job-name=mse_mrp_3_1
#SBATCH --time=36:00:00
#SBATCH --mail-user=lauren.a.kennedy@adelaide.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem-per-cpu=24000
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100


module use /apps/skl/modules/all/
module load R/4.2.3-foss-2021b

R CMD BATCH --no-save --no-restore section_squarederror_fullyobs2.R script_sqrderror_fullobs2_$SLURM_ARRAY_TASK_ID