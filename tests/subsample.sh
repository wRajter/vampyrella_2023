# make a subsample to test the taxonomic assignment

# Variables
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
RAW_DATA="../raw_data"
PR2_DIR="${RAW_DATA}/reference_alignments/pr2"
REF_ALIGN="pr2_version_4.14.0_SSU_UTAX_plus_review_vamp_2022.fasta"
QUERY_SEQ="../raw_data/extracted_18S/Jamy_2019/rDNA/cell/sim99/extracted_18S_seqs_trimmed.fasta"
QUERY_TAB="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}/otu_table.qza"

# ########################
# ## SUBSAMPLE DATABASE ##
# ########################

# # Format reference alignment
# sed 's/;/ /' ${PR2_DIR}/${REF_ALIGN} > ref_align_formated.fasta

# # subsample query sequences to 100 sequences
# ../raw_data/packages/seqtk/seqtk sample -s100 ref_align_formated.fasta 5000 > sub_ref_alignment.fasta

# rm ref_align_formated.fasta

# Create subsampled reference alignment Qiime2 qza file
# qiime tools import \
#   --input-path sub_ref_alignment.fasta \
#   --output-path sub_ref_alignment.qza \
#   --type 'FeatureData[Sequence]'


# # Format taxonomic path from the fasta file into tsv file
# grep '>' sub_ref_alignment.fasta | sed 's/[>*]//g' | \
# sed 's/tax=//g' | \
# sed 's/\s/\t/g' | \
# sed 's/:/__/g' | \
# sed 's/,/;/g' | \
# sed 's/ //g' | \
# sed 's/;/& /g' > sub_ref_taxonomy.tsv

# # Add header
# sed  -i '1i Feature ID\tTaxon' sub_ref_taxonomy.tsv

# # Converting taxonomic path tsv file into qiime qza file
# qiime tools import \
#   --type 'FeatureData[Taxonomy]' \
#   --input-path sub_ref_taxonomy.tsv \
#   --output-path sub_ref_taxonomy.qza



###############################
## SUBSAMPLE QUERY SEQUENCES ##
###############################

rm -f sub_seqs_interm.fasta \
      sub_seqs.fasta \
      sub_seqs.qza

# subsample database sequences to 100 sequences
../raw_data/packages/seqtk/seqtk sample -s100 ${QUERY_SEQ} 200 > sub_seqs_interm.fasta
../raw_data/packages/seqtk/seqtk trimfq -L 1400 sub_seqs_interm.fasta > sub_seqs.fasta


rm sub_seqs_interm.fasta

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path sub_seqs.fasta \
  --output-path sub_seqs.qza
