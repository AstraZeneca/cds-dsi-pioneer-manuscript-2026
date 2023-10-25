#!/bin/bash -l
#
#SBATCH -J early-predict-sbc
#SBATCH --partition=core
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=0.5G
#SBATCH --time=0-00:30:00

module load R-core

Rscript sbc.R -c 12 -n 12