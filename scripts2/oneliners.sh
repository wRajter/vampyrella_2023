#!/bin/bash

# unaligng the already aligned sequences

sed 's/-//g' vamp_refseq_2022.fasta | grep "\S" > vamp_refseq_2022_unalign.fasta


# insert taxpath to the ref alignment
sed 's/>*_/&;tax=k:Eukaryota,d:Rhizaria,p:Cercozoa,c:Endomyxa,o:Vampyrellida,f:,g:,s:  /' review_vamp_ref_2022.fasta

# Extract fasta sequences from a file using a list in another file
seqtk subseq test.fa test.txt
