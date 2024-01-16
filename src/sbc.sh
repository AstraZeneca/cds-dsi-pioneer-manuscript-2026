#!/bin/bash -l
#
#SBATCH -J early-predict-sbc
#SBATCH --partition=core
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-05:00:00
#SBATCH --output=/wscratch/%u/adc-early-predict/log/%x_%j.log

module load R-core

cd .. # Get back to root so R uses correct renv project
Rscript src/sbc.R ${SLURM_CPUS_PER_TASK} $* #1000 sbc_ic_ignore --censor-intervals=3,4,6,9

scontrol show job ${SLURM_JOB_ID}

echo Job ${SLURM_JOB_ID} completed.
