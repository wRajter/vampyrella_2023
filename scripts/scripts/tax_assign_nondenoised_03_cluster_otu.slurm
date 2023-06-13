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



####################
## OTU CLUSTERING ##
####################

# De novo clustering reads into OTUs using vsearch

# Variables
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
RAW_READS_DIR="PacBio/filtered"
CLUST_DIR="otu_clust"
THREADS="32"
THRESHOLD=0.99


###################
## USING VSEARCH ##
###################


# per sample

mkdir -p ${CLUST_DIR}

SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')



for SAMPLE in ${SAMPLES}
do

  # Cleaning
  rm -f ${CLUST_DIR}/otu_${SAMPLE}.fasta

  # Clustering
  vsearch --cluster_fast ${RAW_READS_DIR}/${SAMPLE}.fastq.gz \
          --id ${THRESHOLD} \
          --threads "${THREADS}" \
          --consout ${CLUST_DIR}/otu_${SAMPLE}.fasta

done
