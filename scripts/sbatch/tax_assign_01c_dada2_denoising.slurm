#!/bin/bash
#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --partition=mesca
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500gb
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=dada2_job_%j.out

# Activate your Conda environment
source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

# Run your R script
Rscript tax_assign_01c_dada2_denoising.r
