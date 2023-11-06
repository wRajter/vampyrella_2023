#!/bin/bash

# Formatting reference sequences and taxonomy data for taxonomic assignement.
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
RAW_DATA="../../raw_data"
PR2_DIR="${RAW_DATA}/reference_alignments/pr2"
REF_ALIGN="pr2_version_4.14.0_SSU_UTAX_plus_review_vamp_2023.fasta"



#########################################
## CONVERTING FASTA INTO QIIME2 FORMAT ##
#########################################

# Converting reference sequences in fasta format to Qiime2 artifact (gza) format

# Note:
# Each sequence sould have the accession ID as the FASTA ID, e.g.:
# >KY676659.1 Eurycea nerea isolate MO42 cytochrome oxidase subunit 1 (CO1)
# the GenBank ID is the FASTA ID and everything after a space is a description
# In our fasta file we need to substitute the ; character for whitespace:
# >KC174719;tax=k:Eukaryota ----> >KC174719 tax=k:Eukaryota

echo "Reformating reference sequence fasta file and converting it to Qiime2 qza file"

# Cleaning
rm -f ${PR2_DIR}/ref_align_2023_formated.fasta \
      ${PR2_DIR}/ref_sequences_2023.qza

# Format reference alignment
sed 's/;/ /' ${PR2_DIR}/${REF_ALIGN} > \
             ${PR2_DIR}/ref_align_2023_formated.fasta

# Create reference alignment Qiime2 qza file
qiime tools import \
  --input-path ${PR2_DIR}/ref_align_2023_formated.fasta \
  --output-path ${PR2_DIR}/ref_sequences_2023.qza \
  --type 'FeatureData[Sequence]'



############################################
## CONVERTING TAXONOMY INTO QIIME2 FORMAT ##
############################################

# Converting taxonomic path from the fasta file into Qiime2 artifact (gza) format

# Notes:
# we need to make a TSV mapping file with the genbank/accession ID to the taxonomic string you want to see
# two possible source-formats:
#   TSVTaxonomyFormat which has this header: Feature ID<tab>Taxon on the first line
#   HeaderlessTSVTaxonomyFormat which has no header.
# They both are just a TSV file with IDs, tab, then taxonomy string (usually delimited by ;)

# Correct format example:
# Feature ID	Taxon
# 229854  k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
# 367523	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__
# 203525	k__Bacteria; p__OP11; c__OP11-1; o__; f__; g__; s__

# Our format:
# >MF039355.1.1540_U tax=k:Eukaryota,d:Apusozoa,p:Apusomonadidae,c:Apusomonadidae_Group-4,o:Apusomonadidae_Group-4_X,f:Apusomonadidae_Group-4_XX,g:Collodictyon,s:Collodictyon_triciliatum
# >KC174719 tax=k:Eukaryota,d:Rhizaria,p:Cercozoa,c:Endomyxa,o:Vampyrellida,f:Leptophryidae,g:Leptophrys,s:Leptophrys_solarihaesitans_strain_LAK2013
# >HE609038_LV03 tax=k:Eukaryota,d:Rhizaria,p:Cercozoa,c:Endomyxa,o:Vampyrellida,f:Leptophryidae,g:*,s:Leptophrys_vorax

# We need to:
# delete the > character
# delete *
# delete tax=
# replace whitespace with tab
# replace the : character with __
# replace , to ;
# in the second column add whitespace after ;

echo "Reformating taxonomic path from fasta file and converting it to Qiime2 qza file"

# Cleaning
rm -f ${PR2_DIR}/ref_taxonomy_2023.qza

# Format taxonomic path from the fasta file into tsv file
grep '>' ${PR2_DIR}/ref_align_2023_formated.fasta | sed 's/[>*]//g' | \
sed 's/tax=//g' | \
sed 's/\s/\t/g' | \
sed 's/:/__/g' | \
sed 's/,/;/g' | \
sed 's/ //g' | \
sed 's/;/& /g' > ${PR2_DIR}/ref_taxonomy_2023_formatted.tsv

Add header
sed  -i '1i Feature ID\tTaxon' ${PR2_DIR}/ref_taxonomy_2023_formatted.tsv

# Converting taxonomic path tsv file into qiime qza file
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path ${PR2_DIR}/ref_taxonomy_2023_formatted.tsv \
  --output-path ${PR2_DIR}/ref_taxonomy_2023.qza

rm -f ${PR2_DIR}/ref_taxonomy_2023_formatted.tsv \
      ${PR2_DIR}/ref_align_2023_formated.fasta
