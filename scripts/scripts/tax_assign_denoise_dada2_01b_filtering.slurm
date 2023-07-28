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

#############
### Setup ###
#############

library(dada2)
library(stringi)

#################
### Variables ###
#################

cell <- "cell"
marker <- "rDNA"
project <- "Jamy_2022"
raw_data <- "raw_data"
min_len <- 3000
max_len <- 6000
raw_reads_suffix <- ".fastq.gz"
# Primers
forward <- "GGCAAGTCTGGTGCCAG"
reverse <- "GACGAGGCATTTGGCTACCTT"
# long fragment:
# F3nd: GGCAAGTCTGGTGCCAG
# R21: GACGAGGCATTTGGCTACCTT
# 18S:
# f1: CTGGTTGATYCTGCCAGT
# EukBr: TGATCCTTCTGCAGGTTCACCTAC


#################
### Functions ###
#################

# Dada2 reverse complement function
rc <- dada2:::rc

# Specify path to the raw reads
path <- sprintf("%s/PacBio/%s_%s/%s", raw_data, project, marker, cell)
path_raw_reads <- sprintf("%s/raw", path)

# Getting names and the path to the raw PacBio reads
raw_seqs <- list.files(path = path_raw_reads)
raw_seqs_names <- gsub(raw_reads_suffix, "", raw_seqs)


####################
#### Filtering  ####
####################

# Input and output files
nop <- file.path(path, "noprimers", sprintf("%s_subreads_nonprimer.fastq.qz", raw_seqs_names))
dir.create(sprintf("%s/filtered", path))
filt <- file.path(path, "filtered", basename(fn))


track <- dada2:::filterAndTrim(nop,
                               filt,
                               minQ=NULL,
                               minLen=min_len,
                               maxLen=max_len,
                               maxN=0,
                               rm.phix=FALSE,
                               maxEE=2,
                               verbose=TRUE)

print("filtering done")

# Checking the number of sequences filtered out
stats <- cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2])
# Save
write.table(stats,file="reads_num_track.csv", row.names=TRUE)
file.copy("reads_num_track.csv", sprintf("%s/denoise/%s/%s/%s/dada2/reads_num_track.csv", raw_data, project, marker, cell))
file.remove("reads_num_track.csv")
