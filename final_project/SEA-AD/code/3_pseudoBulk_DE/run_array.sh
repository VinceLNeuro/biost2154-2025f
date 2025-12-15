#!/usr/bin/env bash

######## Slurm resource allocation ########
#SBATCH --job-name=3_pseudoBulk_DE_array
#SBATCH --cluster=teach
#SBATCH --account=biost2154-2025f
#SBATCH --time=4:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8   # Reduced cores per job
#SBATCH --mem=64GB          # 64GB per cell type should be sufficient
#SBATCH --array=1-24        # 24 cell types
#SBATCH --mail-user=til177@pitt.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --output=./array_logs/%x_slurm%A_%a.out

######## Load software into environment ########
set -ev
module purge
module load r/4.5.0

# Record the start time
start_time=$(date +%s)

# Run R script with array task ID
Rscript --vanilla "3_pseudoBulk_DE_array.R" ${SLURM_ARRAY_TASK_ID}

# Record the end time
end_time=$(date +%s)

# Calculate and record the duration
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
