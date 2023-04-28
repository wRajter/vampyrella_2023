#!/usr/bin/env -S Rscript --vanilla

# Setup
library(dada2)


# Specify paths and primers
path <- "PacBio"

# Getting the path to the raw PacBio reads
raw_seqs <- list.files(path = path)
raw_seqs_names <- gsub(".fastq.gz","",raw_seqs)

fn <- c()

for (i in raw_seqs) {
  fn <- append(fn, sprintf("%s/%s", path, i))
}

# Output path for fastq with trimmed primers
dir.create(sprintf("%s/noprimers", path))

nop <-c()

for (i in raw_seqs_names) {
  nop <- append(nop, sprintf("%s/noprimers/%s_subreads_nonprimer.fastq", path, i))
}


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
