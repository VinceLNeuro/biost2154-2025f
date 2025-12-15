#!/usr/bin/env bash

######## Slurm resource allocation ########
#SBATCH --job-name=subset_VLMC
#SBATCH --cluster=teach
#SBATCH --account=biost2154-2025f
#SBATCH --time=4:00:00

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=23   # Reduced cores per job
#SBATCH --mem=184GB          # 184GB per cell type should be sufficient
#SBATCH --mail-user=til177@pitt.edu
#SBATCH --mail-type=END,FAIL

#SBATCH --output=./%x_slurm%A.out

######## Load software into environment ########
module purge
# module load r/4.5.0
source ~/conda_init.sh
conda activate singleCell

set -ev
# Confirm Python version and conda environment
echo "Using Python: $(which python)"
python -V
echo "Conda env: $CONDA_DEFAULT_ENV"


start_time=$(date +%s)
# Rscript --vanilla "3_pseudoBulk_DE_array.R" ${SLURM_ARRAY_TASK_ID}
python 0_subset_VLMC.py
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Duration: $duration seconds"
