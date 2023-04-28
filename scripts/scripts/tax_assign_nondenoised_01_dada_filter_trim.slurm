#!/bin/bash
#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=12
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=%x.o%A_%a

# Activate your Conda environment
source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

# Run your R script
Rscript tax_assign_nondenoised_01_dada_filter_trim_slurm.R
