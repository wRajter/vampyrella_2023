#!/usr/bin/env -S Rscript --vanilla

library(dada2)

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
# Script Config #
#################

min_len <- 700
max_len <- 2500
cell <- "cellCombined"
marker <- "Full18S"
project <- "Suthaus_2022"
raw_data <- make_path("..", "..", "raw_data")
noprimers_path = make_path(raw_data, "dada2", project, marker, "noprimers")
noprimers_files <- list.files(path = noprimers_path)
noprimers_paths <- make_path(noprimers_path, noprimers_files)
filt_path <- make_path(raw_data, "dada2", project, marker, "filtered")
filt_paths <- make_path(filt_path, noprimers_files)


# Ensure output directory exists
create_dir(filt_path)

track <- dada2:::filterAndTrim(noprimers_paths,
                               filt_paths,
                               minQ=2,
                               minLen=min_len,
                               maxLen=max_len,
                               maxN=0,
                               rm.phix=FALSE,
                               maxEE=2,
                               verbose=TRUE)

print("filtering done")
