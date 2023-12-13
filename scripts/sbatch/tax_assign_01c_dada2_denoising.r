#!/usr/bin/env -S Rscript --vanilla

# libraries
library(dada2)

#############
# Variables #
#############

marker <- 'rDNA'
project <- 'Jamy_2022'

####################
# Helper Functions #
####################

make_path <- function(...) {
  file.path(..., fsep = '/')
}

# Function to measure and print the duration
measure_duration <- function(start_time, end_time, step_name) {
  duration <- end_time - start_time
  cat(sprintf("Time taken for %s: %s\n", step_name, duration))
}

################
# Define paths #
################

raw_data_path <- make_path('..', '..', 'raw_data')
path_filt_seq <- make_path(raw_data_path, 'dada2', project, marker, 'filtered')

samples <- list.files(path = path_filt_seq, full.names = TRUE)
test_sample_path <- make_path(path_filt_seq, 'ERR6454467.fastq.gz')

output_dir_path <- make_path(raw_data_path, 'dada2', project, marker, 'denoised')

# Define the log file path
log_file <- make_path(output_dir_path, 'log_file.log')

# Start logging
sink(log_file)
sink(log_file, type = "message")

#############
# Denoising #
#############

print("Working on sample:")
print(basename(test_sample_path))


# Start timer for dereplication
start_time <- Sys.time()

# Dereplicate
drp2 <- derepFastq(test_sample_path, verbose=TRUE)

# End timer for dereplication and print duration
end_time <- Sys.time()
measure_duration(start_time, end_time, "Dereplication")


# Start timer for error learning
start_time <- Sys.time()

# Learn errors (https://rdrr.io/bioc/dada2/man/learnErrors.html)
err2 <- learnErrors(drp2, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

# End timer for error learning and print duration
end_time <- Sys.time()
measure_duration(start_time, end_time, "Error Learning")


plotErrors(err2, nominalQ=TRUE)
# How to read the plots:
# The error rates for each possible transition (A→C, A→G, …) are shown in the plots.
# Points are the observed error rates for each consensus quality score.
# The black line shows the estimated error rates after convergence of the machine-learning algorithm.
# The red line shows the error rates expected under the nominal definition of the Q-score.
# So if the estimated error rates (black line) are a good fit to the observed rates (points), they overlap,
# and the error rates should drop with increased quality.

file.copy("Rplots.pdf", make_path(output_dir_path, 'ERR6454467_err_plot.pdf'))

# Start timer for denoising
start_time <- Sys.time()

# Denoise
dd2 <- dada(drp2, err=err2, BAND_SIZE=32, multithread=TRUE)
saveRDS(dd2, file = make_path(output_dir_path, 'ERR6454467_err2.rds'))

# End timer for denoising and print duration
end_time <- Sys.time()
measure_duration(start_time, end_time, "Denoising")

# Creating sequence table from the dada-class object
seqtab <- makeSequenceTable(dd2)

# Write DADA sequences to fasta
dada_to_fasta(seqtab, out = make_path(output_dir_path, 'ERR6454467_asv.fasta'), hash = "sha1")

# Remember to close the logging
sink(type = "message")
sink()
