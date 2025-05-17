#!/bin/bash
#SBATCH --job-name=test_eta
#SBATCH --output=eta_%A_%a.out
#SBATCH --error=eta_%A_%a.err
#SBATCH --array=1-21
#SBATCH --time=04:00:00

# Define eta based on the array task ID
eta_values=(30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70)

export eta=${eta_values[$SLURM_ARRAY_TASK_ID-1]}

# Run the R script with eta as an argument
module load r/4.3.1

Rscript grid_search_optimal.R

