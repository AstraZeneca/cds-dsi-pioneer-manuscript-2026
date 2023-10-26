#!/bin/bash -l
#
#SBATCH -J early-predict-sbc
#SBATCH --partition=core
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-05:30:00
#SBATCH --output=../temp/log/sbc.log

module load R-core

cd .. # Get back to root so R uses correct renv project
Rscript src/sbc.R -c ${SLURM_CPUS_PER_TASK} -n 520 --append

scontrol show job ${SLURM_JOB_ID}

echo Job ${SLURM_JOB_ID} completed.