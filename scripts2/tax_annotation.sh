#!/bin/bash -
#
# This shell script splits a multiple fasta file in individual fasta
# files, each file is then pairwised aligned with a reference dataset
# using vsearch. The best identity (in terms of percentage of
# identity) is kept.

# Variables
DATA="../raw_data/stampa"
DATABASE="../raw_data/stampa/reference.fas"
INPUT_FILE="../raw_data/stampa/input.fas"
THRESHOLD="10000"
STAMPA_FOLDER="stampa_analysis"
OUTPUT_PREFIX="fasta."

## Verify the uniqueness of reference sequence names
duplicates=$(grep "^>" "${DATABASE}" | cut -d " " -f 1 | sort -d | uniq -d)
if [[ "${duplicates}" ]] ; then
	echo -e "WARNING!\nThe reference database contains duplicated accession numbers\n${duplicates}\n" 1>&2
fi

## Compute the number of jobs
AMPLICON_NUMBER=$(grep -c "^>" "${INPUT_FILE}")
echo "AMPLICON_NUMBER: ${AMPLICON_NUMBER}"
if (( AMPLICON_NUMBER % THRESHOLD == 0 )) ; then
	MAX=$((2 * AMPLICON_NUMBER / THRESHOLD))
else
	MAX=$((2 * AMPLICON_NUMBER / THRESHOLD + 1 ))
fi

## The upper limit to the number of amplicons is 10,000 * THRESHOLD / 2
if [[ ${MAX} -gt 10000 ]] ; then
	echo -e "Too many amplicons!\nChange the threshold value or further split your dataset." 1>&2
	exit 1
fi

## Remove old analysis and create a work folder
if  [[ -d "${DATA}/${STAMPA_FOLDER}" ]] ; then
	echo "Removing old stampa analysis." 1>&2
	rm -rf "${DATA}/${STAMPA_FOLDER}/"
fi
mkdir "${DATA}/${STAMPA_FOLDER}"
cd "${DATA}/${STAMPA_FOLDER}/"
pwd
ls

## Split the input fasta file into chuncks (convert vsearch abundance-style to swarm style)
split --numeric-suffixes \
    --lines="${THRESHOLD}" \
    --suffix-length=5 \
    <(sed -e '/^>/ s/;size=/_/' -e '/^>/ s/;$//' "../input.fas") "${OUTPUT_PREFIX}"
