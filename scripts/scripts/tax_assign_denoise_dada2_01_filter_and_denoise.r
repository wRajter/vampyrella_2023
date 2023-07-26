#!/usr/bin/env -S Rscript --vanilla

#############
### Setup ###
#############

library(dada2)
library(stringi)
library(openssl)

#################
### Variables ###
#################

cell <- "cell"
marker <- "rDNA"
project <- "Suthaus_2022"
raw_data <- "../../raw_data"
min_len <- 3000
max_len <- 6000
raw_reads_suffix <- ".hifi_reads.fastq.gz"
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

# dada_to_fasta (https://github.com/vmikk/metagMisc/blob/master/R/dada_to_fasta.R)
dada_to_fasta <- function(seqtab, out = "DADA2.fasta", hash = "sha1", ...){
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

# Dada2 reverse complement function
rc <- dada2:::rc


########################
### Trimming primers ###
########################

# Specify path to the raw reads
path <- sprintf("%s/PacBio/%s_%s/%s/test", raw_data, project, marker, cell)
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


####################
#### Filtering  ####
####################

dir.create(sprintf("%s/filtered", path))
filt <- file.path(path, "filtered", basename(fn))


track <- dada2:::filterAndTrim(nop,
                               filt,
                               minQ=3,
                               minLen=min_len,
                               maxLen=max_len,
                               maxN=0,
                               rm.phix=FALSE,
                               maxEE=4,
                               verbose=TRUE)

print("filtering done")

# Checking the number of sequences filtered out
stats <- cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2])
# Save
write.table(stats,file="reads_num_track.csv", row.names=TRUE)
file.copy("reads_num_track.csv", sprintf("%s/denoise/%s/%s/%s/dada2/reads_num_track.csv", raw_data, project, marker, cell))
file.remove("reads_num_track.csv")



# #################
# ### Denoising ###
# #################

# Specify paths
path_filt_seq <- sprintf(sprintf("%s/filtered", path))
print("Filtered fastq files:")
list.files(path_filt_seq)

# All sample paths
samples <- list.files(path = path_filt_seq, full.names = TRUE)

# Looping through the samples
for (sample in samples)
{
  # Printing current sample
  sample_path <- sprintf("%s/filtered", path)
  sample_name <- stri_replace_all_regex(sample,
                                        pattern=c(sample_path, raw_reads_suffix),
                                        replacement=c("", ""),
                                        vectorize=FALSE)
  print("Working on sample:")
  print(sample_name)

  # Dereplicate
  drp2 <- derepFastq(sample, verbose=TRUE)

  # Learn errors (https://rdrr.io/bioc/dada2/man/learnErrors.html)
  err2 <- learnErrors(drp2, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

  # saveRDS(err2, file = sprintf("%s/denoise/%s/%s/%s/%s_err2.rds", raw_data, project, marker, cell, sample_name))

  plotErrors(err2, nominalQ=TRUE)
  # How to read the plots:
  # The error rates for each possible transition (A→C, A→G, …) are shown in the plots.
  # Points are the observed error rates for each consensus quality score.
  # The black line shows the estimated error rates after convergence of the machine-learning algorithm.
  # The red line shows the error rates expected under the nominal definition of the Q-score.
  # So if the estimated error rates (black line) are a good fit to the observed rates (points), they overlap,
  # and the error rates should drop with increased quality.

  file.copy("Rplots.pdf", sprintf("%s/denoise/%s/%s/%s/dada2/%s_err_plot.pdf", raw_data, project, marker, cell, sample_name))


  # Denoise
  dd2 <- dada(drp2, err=err2, BAND_SIZE=32, multithread=TRUE)
  saveRDS(dd2, file = sprintf("%s/denoise/%s/%s/%s/dada2/%s_err2.rds", raw_data, project, marker, cell, sample_name))

  # Creating sequence table from the dada-class object
  seqtab <- makeSequenceTable(dd2)

  # Write DADA sequences to fasta
  dada_to_fasta(seqtab, out = sprintf("%s/denoise/%s/%s/%s/dada2/%s_asv.fasta", raw_data, project, marker, cell, sample_name), hash = "sha1")
}
