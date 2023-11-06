#!/usr/bin/env -S Rscript --vanilla


library(dada2)
library(stringi)

####################
# Helper Functions #
####################

make_path <- function(...) {
  file.path(..., fsep = '/')
}

create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

#################
# Main Function #
#################

trim_reads <- function(raw_data_path, project, marker, cell, fwd_primer, rev_primer) {
  # Define paths
  raw_reads_path <- make_path(raw_data_path, "PacBio", project, marker, cell)
  dada2_path <- make_path(raw_data_path, "dada2")
  dada2_noprimers_path <- make_path(dada2_path, project, marker, "noprimers")

  # Ensure output directory exists
  create_dir(dada2_noprimers_path)

  # List raw read files
  individual_raw_reads_paths <- list.files(raw_reads_path, pattern = "*.fastq.gz", full.names = TRUE)
  output_files <- basename(individual_raw_reads_paths)
  # Prepare output file names
  dada2_output_files_paths <- make_path(dada2_noprimers_path, output_files)

  # Remove primers and orient reads
  removePrimers(individual_raw_reads_paths, dada2_output_files_paths,
                primer.fwd = fwd_primer, primer.rev = dada2:::rc(rev_primer),
                orient = TRUE, verbose = TRUE)
}

#################
# Script Config #
#################

cell <- "cellCombined"
marker <- "Full18S"
project <- "Suthaus_2022"
raw_data <- "../../raw_data"

# Primers
forward_primer <- "CTGGTTGATYCTGCCAGT"
reverse_primer <- "TGATCCTTCTGCAGGTTCACCTAC"
# long fragment:
# F3nd: GGCAAGTCTGGTGCCAG
# R21: GACGAGGCATTTGGCTACCTT
# 18S:
# f1: CTGGTTGATYCTGCCAGT
# EukBr: TGATCCTTCTGCAGGTTCACCTAC

# Run the main function
trim_reads(raw_data, project, marker, cell, forward_primer, reverse_primer)
