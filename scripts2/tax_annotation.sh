#!/bin/bash -


# variables
VSEARCH="../packages/vsearch-2.22.1/bin/vsearch"
INPUT_FILE="../raw_data/reads/ASV/X17007_asv2.fasta"
DATABASE="../raw_data/stampa/ref_test.fas"
INDEX=$(printf "%05g\n" $(( $LSB_JOBINDEX - 1 )))
QUERIES="fasta.${INDEX}"
ASSIGNMENTS="${QUERIES/fasta./hits.}"
IDENTITY="0.5"
MAXREJECTS=32
THREADS=16
NULL="/dev/null"


# compare environmental sequences to known reference sequences
"${VSEARCH}" --usearch_global "${INPUT_FILE}" \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --userfields query+id1+target \
    --maxaccepts 0 \
    --maxrejects "${MAXREJECTS}" \
    --top_hits_only \
    --output_no_hits \
    --db "${DATABASE}" \
    --id "${IDENTITY}" \
    --iddef 1 \
    --userout "${ASSIGNMENTS}" > "${NULL}" 2> "${NULL}"

rm -f "${QUERIES}"