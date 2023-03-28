#!/usr/bin/env -S Rscript --vanilla

# Setup
library(dada2)

print(paste("your working directory is: ", getwd()))


# Variables
cell <- "cell"
marker <- "rDNA"
project <- "Jamy_2019_rDNA"
raw_data <- "../../raw_data"

# Specify paths and primers
path <- sprintf("%s/PacBio/%s/%s", raw_data, project, cell)

# TODO: gunzip the raw reads in this script
# TODO: create loop to get the path automatically
# Path to the raw reads
fn <- c(sprintf("%s/ERR2355431_subreads.fastq", path),
        sprintf("%s/ERR2355432_subreads.fastq", path),
        sprintf("%s/ERR2355433_subreads.fastq", path))



# TODO: create loop to get the path automatically
# Output path for fastq with trimmed primers
nop <- c(sprintf("%s/noprimers/ERR2355431_subreads_nonprimer.fastq", path),
         sprintf("%s/noprimers/ERR2355432_subreads_nonprimer.fastq", path),
         sprintf("%s/noprimers/ERR2355433_subreads_nonprimer.fastq", path))



# Primers
F3nd <- "GGCAAGTCTGGTGCCAG"
R21 <- "GACGAGGCATTTGGCTACCTT"

# Dada2 reverse complement function
rc <- dada2:::rc


#### Remove primers and orient reads
prim <- dada2:::removePrimers(fn, nop, primer.fwd=F3nd, primer.rev=dada2:::rc(R21), orient=TRUE, verbose=TRUE)


#### Filter
filt <- file.path(path, "filtered", basename(fn))

track <- dada2:::filterAndTrim(nop, filt, minQ=3, minLen=3000, maxLen=6000, maxN=0, rm.phix=FALSE, maxEE=4, verbose=TRUE)
