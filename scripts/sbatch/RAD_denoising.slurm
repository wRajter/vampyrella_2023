#!/bin/bash
#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=12
#SBATCH --mem=46g
#SBATCH --time=300:00:00
#SBATCH --account=ag-hess
#SBATCH --output=out_%x.o%A_%a.log

# Load Julia module
module load julia

# Activate your environment (if needed)
# conda activate your_conda_environment

# Run the Julia script
julia RAD_denoising.jl
