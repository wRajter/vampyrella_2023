#!/bin/bash

# create a tsv table with ASVs and their stats:
  # ASV_ID
  # Sample_ID
  # Abundance
  # Taxonomy
  # Percent similarity
  # Sequence


# Variables:
CELL="cellCombined"
MARKER="Full18S"
PROJECT_DIR="/home/lubo/code/wRajter/vampyrella_2023"
#for local machine at uni: "/home/lubomir/projects/vampyrella_2023"
QIIME_DIR="${PROJECT_DIR}/raw_data/qiime_output"
DENOISE_DIR="${QIIME_DIR}/dada2_denoise"
POST_DENOISE_DIR="${QIIME_DIR}/dada2_post_denoise"
ASSIGNMENT_DIR="${QIIME_DIR}/assignment_results/${MARKER}_${CELL}"
TAX_LEVELS=";d__ ;p__ ;c__ ;o__ ;f__ ;g__ ;s__"



# prepare data:
# create feature (ASV) abundance table
echo 'Creating an ASV abundance table...'
rm -f ${POST_DENOISE_DIR}/feature-table.biom
rm -f ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}.*

qiime tools export \
  --input-path ${DENOISE_DIR}/table_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --output-path ${POST_DENOISE_DIR}/
mv ${POST_DENOISE_DIR}/feature-table.biom ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}.biom

biom convert \
  --to-tsv -i ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}.biom \
  -o ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}_incomplete.tsv
rm ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}.biom
awk '{for(i=2; i<=NF; i++){ sum+= $i } \
    print $1"\t"sum"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; sum = 0}' \
    ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}_incomplete.tsv > ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}.tsv
rm ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}_incomplete.tsv


# create a metadata table that contains ID, DNA sequence and a taxonomy assignment with percent similarity for each AVS:
echo 'Creating a metadata table with taxonomy assignment and DNA sequence for each ASV...'
rm -f ${POST_DENOISE_DIR}/viz.qzv
qiime metadata tabulate \
  --m-input-file ${DENOISE_DIR}/asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --m-input-file ${ASSIGNMENT_DIR}/taxonomy_ALLSAMPLES.tsv \
  --o-visualization ${POST_DENOISE_DIR}/viz.qzv

rm -fr ${POST_DENOISE_DIR}/metadata
rm -f ${POST_DENOISE_DIR}/metadata_incomplete*
qiime tools export \
  --output-path ${POST_DENOISE_DIR}/metadata \
  --input-path ${POST_DENOISE_DIR}/viz.qzv
mv ${POST_DENOISE_DIR}/metadata/metadata.tsv ${POST_DENOISE_DIR}/metadata.tsv
mv ${POST_DENOISE_DIR}/metadata.tsv ${POST_DENOISE_DIR}/metadata_incomplete.tsv
rm -r ${POST_DENOISE_DIR}/metadata/
rm -f ${POST_DENOISE_DIR}/viz.qzv

awk '{print $1"\t"$2"\t"$4"\t"$3}' ${POST_DENOISE_DIR}/metadata_incomplete.tsv > ${POST_DENOISE_DIR}/metadata_incomplete2.tsv
rm -f ${POST_DENOISE_DIR}/metadata_incomplete.tsv

# split taxonomy into tab separated fields and delete a header
rm -f ${POST_DENOISE_DIR}/metadata.tsv

echo 'Spliting taxonomy path into separete columns...'
for TAX_LEVEL in ${TAX_LEVELS}
do
  sed -i "s/${TAX_LEVEL}/ /" ${POST_DENOISE_DIR}/metadata_incomplete2.tsv
done
sed -i 's/k__//' ${POST_DENOISE_DIR}/metadata_incomplete2.tsv

# finding the widest line:
max=$(awk 'max < NF { max = NF } END { print max }' ${POST_DENOISE_DIR}/metadata_incomplete2.tsv)
# using the input from the previous command (widest line) when filling out the mising columns:
awk -v max=$max '{ for(i=NF+1; i<=max; i++) $i = "NA"; print }' ${POST_DENOISE_DIR}/metadata_incomplete2.tsv > ${POST_DENOISE_DIR}/metadata.tsv
rm -f ${POST_DENOISE_DIR}/metadata_incomplete*


# merge abundance table and metadata table using sort and join commands
# deleting and adding a new header
echo 'Merging ASV abundance table and metadata table...'
rm -f ${POST_DENOISE_DIR}/table_1
rm -f ${POST_DENOISE_DIR}/table_2
rm -f ${POST_DENOISE_DIR}/asv_stats_table.tsv

tail -n +3 ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}.tsv | sort > ${POST_DENOISE_DIR}/table_1
tail -n +3 ${POST_DENOISE_DIR}/metadata.tsv | sort > ${POST_DENOISE_DIR}/table_2
join -1 1 -2 1 ${POST_DENOISE_DIR}/table_1 ${POST_DENOISE_DIR}/table_2 | \
sed '1i ID\tTotal\tA3\tMock\tNH1\tNH4\tSim17\tSim22\tTh16\tTh38\tTh40\tX17007\tSequence\tPercent_similarity\tKingdom\tDomain\tPhyllum\tClass\tOrder\tFamily\tGenus\tSpecies' > ${POST_DENOISE_DIR}/asv_stats_table.tsv
rm -f ${POST_DENOISE_DIR}/table_1 ${POST_DENOISE_DIR}/table_2
rm -f ${POST_DENOISE_DIR}/metadata.tsv
rm -f ${POST_DENOISE_DIR}/asv_table_${MARKER}_${CELL}.tsv

# checking if all the rows have the same number of the fields
echo 'Checking if all the rows have the same number of the columns...'
field_num=$(awk 'max < NF { max = NF } END { print max }' ${POST_DENOISE_DIR}/asv_stats_table.tsv)
field_num_output=$(awk -v field_num=$field_num 'NF == field_num {print NR}' ${POST_DENOISE_DIR}/asv_stats_table.tsv | wc -l)
num_asvs=$(cat ${POST_DENOISE_DIR}/asv_stats_table.tsv | wc -l)

if [ ${num_asvs} -eq ${field_num_output} ]
then
  echo 'Great, your asv stats table has the correct number of rows :)'
else
  echo 'Oops, it seems that the asv stats table has inconsistent number of columns'
fi

# rearrange fields
rm -f ${POST_DENOISE_DIR}/asv_stats_${MARKER}_${CELL}.tsv
echo 'Rearranging columns in the final ASV stats table...'
sed -e "s/\r//g" ${POST_DENOISE_DIR}/asv_stats_table.tsv | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$14"\t"$13}' > ${POST_DENOISE_DIR}/asv_summary_table_${MARKER}_${CELL}.tsv

rm -f ${POST_DENOISE_DIR}/asv_stats_table.tsv
echo "Done. Your ASV stats table is saved as asv_stats_${MARKER}_${CELL}.tsv"
