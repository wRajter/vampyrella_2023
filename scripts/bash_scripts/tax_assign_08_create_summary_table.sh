#!/bin/bash

# create a tsv table with ASVs and their stats:
  # ASV_ID
  # Sample_ID
  # Abundance
  # Taxonomy
  # Percent similarity
  # Sequence


# Variables:
CELL="cell1"
MARKER="Full18S"
SIM="sim97"
RAW_DATA="../../raw_data"
FILT_OTU_DIR="${RAW_DATA}/OTU_filtered/${MARKER}/${CELL}/${SIM}"
OTU_SUMMARY_DIR="${RAW_DATA}/OTU_summary_tables"
ASSIGNMENT_DIR="${RAW_DATA}/tax_assign_results/${MARKER}/${CELL}/${SIM}"
TAX_LEVELS=";d__ ;p__ ;c__ ;o__ ;f__ ;g__ ;s__"
RAW_READS_DIR="${RAW_DATA}/PacBio/Suthaus${MARKER}/${CELL}"
SAMPLE_LABELS=$(ls ${RAW_READS_DIR}/*reads.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '_' '{ print $1 }' | \
          awk 'BEGIN { ORS = "\t" } { print $1 }')
TAX_LABELS="Kingdom Domain  Phyllum Class Order Family  Genus Species"


#####################
## PREPARE TABLE 1 ##
#####################

# Create feature (OTU) abundance table
echo 'Creating OTU abundance table...'

# Cleaning
mkdir -p ${OTU_SUMMARY_DIR}/
rm -f ${OTU_SUMMARY_DIR}/feature-table.biom \
      ${OTU_SUMMARY_DIR}/otu_table_incomplete.tsv \
      ${OTU_SUMMARY_DIR}/otu_table.tsv

# Creating the feature-table.biom file
qiime tools export \
  --input-path ${FILT_OTU_DIR}/otu_table_filtered.qza \
  --output-path ${OTU_SUMMARY_DIR}/

# Converting biome format to tsv
biom convert \
  --to-tsv -i ${OTU_SUMMARY_DIR}/feature-table.biom \
  -o ${OTU_SUMMARY_DIR}/otu_table_incomplete.tsv

rm ${OTU_SUMMARY_DIR}/feature-table.biom

# Sum all abundances across samples for each OTU
# Add the sum as the last column representing total abundance
awk '{ for(i=2; i<=NF; i++){ sum+= $i } \
    print $0"\t"sum; sum = 0 }' \
    ${OTU_SUMMARY_DIR}/otu_table_incomplete.tsv > ${OTU_SUMMARY_DIR}/otu_table.tsv
rm ${OTU_SUMMARY_DIR}/otu_table_incomplete.tsv



#####################
## PREPARE TABLE 2 ##
#####################

# Create a metadata table that contains ID, DNA sequence and a taxonomy assignment with percent similarity for each OTU:
echo 'Creating a metadata table with taxonomy assignment and DNA sequence for each OTU...'

# Cleaning
rm -f ${OTU_SUMMARY_DIR}/viz.qzv

# Merging assignment results and otu sequences into one table
qiime metadata tabulate \
  --m-input-file ${FILT_OTU_DIR}/otu_seqs_filtered.qza \
  --m-input-file ${ASSIGNMENT_DIR}/taxonomy.tsv \
  --o-visualization ${OTU_SUMMARY_DIR}/viz.qzv

# Cleaning
rm -fr ${OTU_SUMMARY_DIR}/metadata
rm -f ${OTU_SUMMARY_DIR}/metadata_incomplete*

qiime tools export \
  --output-path ${OTU_SUMMARY_DIR}/metadata \
  --input-path ${OTU_SUMMARY_DIR}/viz.qzv
mv ${OTU_SUMMARY_DIR}/metadata/metadata.tsv ${OTU_SUMMARY_DIR}/metadata_incomplete.tsv
rm -r ${OTU_SUMMARY_DIR}/metadata/
rm -f ${OTU_SUMMARY_DIR}/viz.qzv

# Change the column order: ID, seqs, consensus, taxonomic path
awk 'BEGIN { OFS="\t" } \
    { print $1,$2,$4,$3 }' \
    ${OTU_SUMMARY_DIR}/metadata_incomplete.tsv > \
    ${OTU_SUMMARY_DIR}/metadata_incomplete2.tsv
rm -f ${OTU_SUMMARY_DIR}/metadata_incomplete.tsv

# Split taxonomy into tab separated fields and delete a header

# Cleaning
rm -f ${POST_DENOISE_DIR}/metadata.tsv

echo 'Spliting taxonomic path into separete columns...'
for TAX_LEVEL in ${TAX_LEVELS}
do
  sed -i "s/${TAX_LEVEL}/ /" ${OTU_SUMMARY_DIR}/metadata_incomplete2.tsv
done
sed -i 's/k__//' ${OTU_SUMMARY_DIR}/metadata_incomplete2.tsv

# Finding the widest record/row:
max=$(awk 'max < NF { max = NF } \
      END { print max }' \
      ${OTU_SUMMARY_DIR}/metadata_incomplete2.tsv)
# Using the input from the previous command (widest line) when filling out the mising columns in taxonomic path:
awk -v max=$max '{ for(i=NF+1; i<=max; i++) \
                  $i = "NA"; print }' \
                  ${OTU_SUMMARY_DIR}/metadata_incomplete2.tsv > \
                  ${OTU_SUMMARY_DIR}/metadata.tsv
rm -f ${OTU_SUMMARY_DIR}/metadata_incomplete2.tsv



#########################
## MERGE TABLE 1 and 2 ##
#########################

# Merge abundance table and metadata table using sort and join commands
# Deleting and adding a new header
echo 'Merging OTU abundance table and metadata table...'

# Cleaning
rm -f ${OTU_SUMMARY_DIR}/table_*
rm -f ${OTU_SUMMARY_DIR}/header
rm -f ${OTU_SUMMARY_DIR}/otu_summary_table.tsv

tail -n +3 ${OTU_SUMMARY_DIR}/otu_table.tsv | sort > ${OTU_SUMMARY_DIR}/table_1
tail -n +3 ${OTU_SUMMARY_DIR}/metadata.tsv | sort > ${OTU_SUMMARY_DIR}/table_2
join -1 1 -2 1 ${OTU_SUMMARY_DIR}/table_1 ${OTU_SUMMARY_DIR}/table_2 > ${OTU_SUMMARY_DIR}/table_3

# Creating a header
echo "ID  ${SAMPLE_LABELS}  Total Sequence  Percent_similarity  ${TAX_LABELS}" > ${OTU_SUMMARY_DIR}/header

# Combine the header with the table
cat ${OTU_SUMMARY_DIR}/header ${OTU_SUMMARY_DIR}/table_3 > ${OTU_SUMMARY_DIR}/otu_summary_table.tsv

rm -f ${OTU_SUMMARY_DIR}/table_* \
      ${OTU_SUMMARY_DIR}/metadata.tsv \
      ${OTU_SUMMARY_DIR}/otu_table.tsv \
      ${OTU_SUMMARY_DIR}/header

# Checking if all the rows have the same number of the fields
echo 'Checking if all the rows have the same number of the columns...'
field_num=$(awk 'max < NF { max = NF } END { print max }' ${OTU_SUMMARY_DIR}/otu_summary_table.tsv)
field_num_output=$(awk -v field_num=$field_num 'NF == field_num {print NR}' ${OTU_SUMMARY_DIR}/otu_summary_table.tsv | wc -l)
num_asvs=$(cat ${OTU_SUMMARY_DIR}/otu_summary_table.tsv | wc -l)

if [ ${num_asvs} -eq ${field_num_output} ]
then
  echo 'Great, your otu summary table has the correct number of rows :)'
else
  echo 'Oops, it seems that the otu summary table has inconsistent number of columns'
fi

# Rearrange fields
echo 'Rearranging columns in the final ASV stats table...'

# Cleaning
rm -f ${OTU_SUMMARY_DIR}/otu_summary_table_${MARKER}_${CELL}_${SIM}.tsv

# Removing convert Windows line endings to Unix line endings and rearanging the the columns
sed -e "s/\r//g" ${OTU_SUMMARY_DIR}/otu_summary_table.tsv | \
awk 'BEGIN { OFS="\t" } \
           { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$15,$16,$17,$18,$19,$20,$21,$22,$14,$13}' \
              > ${OTU_SUMMARY_DIR}/otu_summary_table_${MARKER}_${CELL}_${SIM}.tsv

rm -f ${OTU_SUMMARY_DIR}/otu_summary_table.tsv
echo "Done. Your OTU summary table is saved as otu_summary_table_${MARKER}_${CELL}_${SIM}.tsv"
