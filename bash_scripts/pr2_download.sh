#!/bin/bash
#
# This shell script downloads and trims the PR2 database for the stampa pipeline
# adjusted from Frederic Mahe's script: https://github.com/frederic-mahe/stampa


# 1. download the UTAX version of the PR2 database
VERSION="4.14.0"
URL="https://github.com/pr2database/pr2database/releases/download"
SOURCE="pr2_version_${VERSION}_SSU_UTAX.fasta"
LOCATION="../raw_data/pr2"

wget -P ${LOCATION} "${URL}/v${VERSION}/${SOURCE}.gz"

gunzip -k ${LOCATION}/${SOURCE}.gz

# 2. extract the 18S from the downloaded PR2 database
PRIMER_F="CTGGTTGATYCTGCCAGT"
PRIMER_R="TGATCCTTCTGCAGGTTCACCTAC"
OUTPUT="${SOURCE/_UTAX*/}_${PRIMER_F}_${PRIMER_R}.fas"
LOG="${OUTPUT/.fas/.log}"
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))
CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH}"

dos2unix < "${LOCATION}/${SOURCE}" | \
    sed '/^>/ s/;tax=k:/ /
         /^>/ s/,[dpcofgs]:/|/g
         /^>/ ! s/U/T/g' | \
    ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOCATION}/${LOG}" | \
    ${CUTADAPT} -a "${PRIMER_R}" -O "${MIN_R}" - 2>> "${LOCATION}/${LOG}" > "${LOCATION}/${OUTPUT}"
