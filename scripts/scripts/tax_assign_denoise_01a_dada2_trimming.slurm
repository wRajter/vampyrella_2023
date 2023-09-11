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


########################
### Trimming primers ###
########################

# Specify path to the raw reads
path <- sprintf("%s/PacBio/%s_%s/%s", raw_data, project, marker, cell)
path_raw_reads <- sprintf("%s/raw", path)


# Getting names and the path to the raw PacBio reads
raw_seqs <- list.files(path = path_raw_reads)
raw_seqs_names <- gsub(raw_reads_suffix, "", raw_seqs)


fn <-list.files(path = path_raw_reads, full.names = TRUE)

# Output path for fastq with trimmed primers
dir.create(sprintf("%s/noprimers", path))

nop <- file.path(path, "noprimers", sprintf("%s_subreads_nonprimer.fastq.qz", raw_seqs_names))

print(fn)

# Remove primers and orient reads
prim <- dada2:::removePrimers(fn, nop, primer.fwd=forward, primer.rev=dada2:::rc(reverse), orient=TRUE, verbose=TRUE)


# Inspect length distribution of the non-primer sequences
lens.fn <- lapply(nop, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)
file.copy("Rplots.pdf", sprintf("%s/denoise/%s/%s/%s/dada2/reads_length_dist_plot.pdf", raw_data, project, marker, cell))
file.remove("Rplots.pdf")
