#!/usr/bin/env -S Rscript --vanilla

# Setup
library(dada2)
library(stringi)

print(paste("your working directory is: ", getwd()))

# Functions
# dada_to_fasta (https://github.com/vmikk/metagMisc/blob/master/R/dada_to_fasta.R)
dada_to_fasta <- function(seqtab, out = "DADA2.fasta", hash = "sha1", ...){

  # require(dada2)
  # require(openssl)

  # prepare sequence names in USEARCH and VSEARCH-style
  seq_uniq <- dada2::getUniques(seqtab)   # integer vector named by unique sequence and valued by abundance.

  if(hash == "sha1"){ hh <- openssl::sha1(names(seq_uniq)) }
  if(hash == "sha256"){ hh <- openssl::sha256(names(seq_uniq)) }
  if(hash == "md5"){ hh <- openssl::md5(names(seq_uniq)) }

  seq_names <- paste(as.character(hh),
                     ";size=",
                     seq_uniq,
                     ";",
                     sep="")

  # Export sequence as fasta
  dada2::uniquesToFasta(seq_uniq, fout = out, ids = seq_names, ...)

  invisible(seq_names)
}


# Variables
cell <- "cell"
marker <- "rDNA"
project <- "Suthaus_2022"
raw_data <- "../../raw_data"

# Specify paths
path_filt_seq <- sprintf("%s/PacBio/%s_rDNA/%s/test", raw_data, project, cell)
print("fastq files:")
list.files(path_filt_seq)

# All sample paths
samples <- list.files(path = path_filt_seq, full.names = TRUE)

# Looping through the samples
for (sample in samples)
{
  # Printing current sample
  sample_path <- sprintf("%s/PacBio/%s_rDNA/%s/test/", raw_data, project, cell)
  sample_suffix <- ".hifi_reads.fastq.gz"
  sample_name <- stri_replace_all_regex(sample,
                                        pattern=c(sample_path, sample_suffix),
                                        replacement=c("", ""),
                                        vectorize=FALSE)
  print("Working on sample:")
  print(sample_name)

  # Dereplicate
  drp2 <- derepFastq(sample, verbose=TRUE)

  # Learn errors (https://rdrr.io/bioc/dada2/man/learnErrors.html)
  err2 <- learnErrors(drp2, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

  saveRDS(err2, file = sprintf("%s/denoise/%s/%s/%s/%s_err2.rds", raw_data, project, marker, cell, sample_name))


  plotErrors(err2, nominalQ=TRUE)
  # How to read the plots:
  # The error rates for each possible transition (A→C, A→G, …) are shown in the plots.
  # Points are the observed error rates for each consensus quality score.
  # The black line shows the estimated error rates after convergence of the machine-learning algorithm.
  # The red line shows the error rates expected under the nominal definition of the Q-score.
  # So if the estimated error rates (black line) are a good fit to the observed rates (points), they overlap,
  # and the error rates should drop with increased quality.

  file.copy("Rplots.pdf", sprintf("%s/denoise/%s/%s/%s/%s_err_plot.pdf", raw_data, project, marker, cell, sample_name))



  # Denoise
  dd2 <- dada(drp2, err=err2, BAND_SIZE=32, multithread=TRUE)
  saveRDS(dd2, file = sprintf("%s/denoise/%s/%s/%s/%s_err2.rds", raw_data, project, marker, cell, sample_name))

  # Inspecting the returned dada-class object
  dd2[2]

  # Creating sequence table from the dada-class object
  seqtab <- makeSequenceTable(dd2)

  # Write DADA sequences to fasta
  dada_to_fasta(seqtab, out = sprintf("%s/denoise/%s/%s/%s/%s_asv.fasta", raw_data, project, marker, cell, sample_name), hash = "sha1")
}
