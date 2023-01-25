#!/bin/bash

# create a tsv table with ASVs and their stats:
  # ASV_ID
  # Sample_ID
  # Biome
  # Abundance
  # Num_Reads
  # Taxonomy
  # Percent similarity
  # Sequence length
  # Sequence


# Variables:
CELL="cell2"
MARKER="Full18S"
PROJECT_DIR="/home/lubo/code/wRajter/vampyrella_2023"
#for local machine at uni: "/home/lubomir/projects/vampyrella_2023"
ASV_DIR="${PROJECT_DIR}/raw_data/qiime_output/ASV"
DADA2_DIR="${PROJECT_DIR}/raw_data/qiime_output/dada2_output"
ASSIGNMENT_DIR="${PROJECT_DIR}/raw_data/qiime_output/assignment_results/${MARKER}_${CELL}"
SAMPLES="A3 Mock"
TAX_LEVELS=";d__ ;p__ ;c__ ;o__ ;f__ ;g__ ;s__"



# prepare data:
# create feature (ASV) abundance table
echo 'Creating an ASV abundance table...'
rm -f ${DADA2_DIR}/feature-table.biom ${DADA2_DIR}/asv_table_${MARKER}_${CELL}.biom ${DADA2_DIR}/asv_table_${MARKER}_${CELL}.tsv
qiime tools export --input-path ${DADA2_DIR}/table_${MARKER}_${CELL}.qza --output-path ${DADA2_DIR}/
mv ${DADA2_DIR}/feature-table.biom ${DADA2_DIR}/asv_table_${MARKER}_${CELL}.biom
biom convert \
  --to-tsv -i ${DADA2_DIR}/asv_table_${MARKER}_${CELL}.biom \
  -o ${DADA2_DIR}/asv_table_${MARKER}_${CELL}_incomplete.tsv
rm ${DADA2_DIR}/asv_table_${MARKER}_${CELL}.biom
awk '{for(i=2; i<=NF; i++){ sum+= $i } \
    print $1"\t"sum"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; sum = 0}' \
    ${DADA2_DIR}/asv_table_${MARKER}_${CELL}_incomplete.tsv > ${DADA2_DIR}/asv_table_${MARKER}_${CELL}.tsv
rm ${DADA2_DIR}/asv_table_${MARKER}_${CELL}_incomplete.tsv


# create a metadata table that contains ID, DNA sequence and a taxonomy assignment with percent similarity for each AVS:
echo 'Creating a metadata table...'
rm -f ${DADA2_DIR}/viz.qzv
qiime metadata tabulate \
  --m-input-file ${DADA2_DIR}/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --m-input-file ${ASSIGNMENT_DIR}/taxonomy_ALLSAMPLES.tsv \
  --o-visualization ${DADA2_DIR}/viz.qzv

rm -fr ${DADA2_DIR}/metadata
qiime tools export --output-path ${DADA2_DIR}/metadata --input-path ${DADA2_DIR}/viz.qzv
mv ${DADA2_DIR}/metadata/metadata.tsv ${DADA2_DIR}/metadata.tsv
mv ${DADA2_DIR}/metadata.tsv ${DADA2_DIR}/metadata_incomplete.tsv
rm -r ${DADA2_DIR}/metadata/
rm -f viz.qzv

awk '{print $1"\t"$2"\t"$4"\t"$3}' ${DADA2_DIR}/metadata_incomplete.tsv > ${DADA2_DIR}/metadata_incomplete2.tsv

# split taxonomy into tab separated fields and delete a header
echo 'Spliting taxonomy path into separete columns...'
for TAX_LEVEL in ${TAX_LEVELS}
do
  sed -i "s/${TAX_LEVEL}/ /" ${DADA2_DIR}/metadata_incomplete2.tsv
done
sed -i 's/k__//' ${DADA2_DIR}/metadata_incomplete2.tsv

# finding the widest line:
max=$(awk 'max < NF { max = NF } END { print max }' ${DADA2_DIR}/metadata_incomplete2.tsv)
# using the input from the previous command (widest line) when filling out the mising columns:
awk -v max=$max '{ for(i=NF+1; i<=max; i++) $i = "NA"; print }' ${DADA2_DIR}/metadata_incomplete2.tsv > ${DADA2_DIR}/metadata.tsv
rm -f ${DADA2_DIR}/metadata_incomplete*


# merge abundance table and metadata table using sort and join commands
# deleting and adding a new header
echo 'Merging ASV abundance table and metadata table...'
rm -f ${DADA2_DIR}/table_1 ${DADA2_DIR}/table_2 ${DADA2_DIR}/asv_stats_table.tsv
tail -n +3 ${DADA2_DIR}/asv_table_Full18S_cell2.tsv | sort > ${DADA2_DIR}/table_1
tail -n +3 ${DADA2_DIR}/metadata.tsv | sort > ${DADA2_DIR}/table_2
join -1 1 -2 1 ${DADA2_DIR}/table_1 ${DADA2_DIR}/table_2 | \
sed '1i ID\tTotal\tA3\tMock\tNH1\tNH4\tSim17\tSim22\tTh16\tTh38\tTh40\tX17007\tSequence\tPercent_similarity\tKingdom\tDomain\tPhyllum\tClass\tOrder\tFamily\tGenus\tSpecies' > ${DADA2_DIR}/asv_stats_table.tsv
rm -f ${DADA2_DIR}/table_1 ${DADA2_DIR}/table_2 ${DADA2_DIR}/metadata.tsv

# checking if all the rows have the same number of the fields
echo 'Checking if all the rows have the same number of the columns...'
field_num=$(awk 'max < NF { max = NF } END { print max }' ${DADA2_DIR}/asv_stats_table.tsv)
field_num_output=$(awk -v field_num=$field_num 'NF == field_num {print NR}' ${DADA2_DIR}/asv_stats_table.tsv | wc -l)
num_asvs=$(cat ${DADA2_DIR}/asv_stats_table.tsv | wc -l)

if [ ${num_asvs} -eq ${field_num_output} ]
then
  echo 'Great, your asv stats table has the correct number of rows :)'
else
  echo 'Oops, it seems that the asv stats table has inconsistent number of columns'
fi

# rearrange fields
rm -f ${DADA2_DIR}/asv_stats_${MARKER}_${CELL}.tsv
echo 'Rearranging columns in the final ASV stats table...'
sed -e "s/\r//g" ${DADA2_DIR}/asv_stats_table.tsv | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$14"\t"$13}' > ${DADA2_DIR}/asv_summary_table_${MARKER}_${CELL}.tsv

rm -f ${DADA2_DIR}/asv_stats_table.tsv
echo "Done. Your ASV stats table is saved as asv_stats_${MARKER}_${CELL}.tsv"
