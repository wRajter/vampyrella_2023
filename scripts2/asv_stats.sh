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
ASV_DIR="${PROJECT_DIR}/raw_data/qiime_output/ASV"
DADA2_DIR="${PROJECT_DIR}/raw_data/qiime_output/dada2_output"
ASSIGNMENT_DIR="${PROJECT_DIR}/raw_data/qiime_output/assignment_results/${MARKER}_${CELL}"
SAMPLES="A3 Mock"
TAX_LEVELS=";d__ ;p__ ;c__ ;o__ ;f__ ;g__ ;s__"


# table with tax assignment (id, taxon, perc._sim.): ${ASSIGNMENT_DIR}/taxonomy_${SAMPLE}.tsv
# table with abundance (id, sampleX, sampleY, ...): ${DADA2_DIR}/asv_table_${MARKER}_${CELL}.tsv


# prepare data:
# create feature (ASV) abundance table
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
rm -f ${DADA2_DIR}/viz.qzv
qiime metadata tabulate \
  --m-input-file ${DADA2_DIR}/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --m-input-file ${ASSIGNMENT_DIR}/taxonomy_ALLSAMPLES.tsv \
  --o-visualization ${DADA2_DIR}/viz.qzv

rm -fr ${DADA2_DIR}/metadata
qiime tools export --output-path ${DADA2_DIR}/metadata --input-path ${DADA2_DIR}/viz.qzv
mv ${DADA2_DIR}/metadata/metadata.tsv ${DADA2_DIR}/metadata.tsv
rm -r ${DADA2_DIR}/metadata/
rm -f viz.qzv

# split taxonomy into tab separated fields and delete a header
# for TAX_LEVEL in ${TAX_LEVELS}
# do
#   sed -i "s/${TAX_LEVEL}/  ${TAX_LEVEL}/" ${DADA2_DIR}/metadata.tsv
# done
# sed -i 's/k__//' ${DADA2_DIR}/metadata.tsv


# merge abundance table and metadata table
# sort and join
# delete headers
# rearrange fields
# new header

rm ${DADA2_DIR}/table_1 ${DADA2_DIR}/table_2 ${DADA2_DIR}/asv_stats_table.tsv
tail -n +3 ${DADA2_DIR}/asv_table_Full18S_cell2.tsv | sort > ${DADA2_DIR}/table_1
tail -n +3 ${DADA2_DIR}/metadata.tsv | sort > ${DADA2_DIR}/table_2
join -1 1 -2 1 ${DADA2_DIR}/table_1 ${DADA2_DIR}/table_2 | \
sed '1i ID\tTotal\tA3\tMock\tNH1\tNH4\tSim17\tSim22\tTh16\tTh38\tTh40\tX17007\tSequence\tTaxonomy\tPercent_similarity' > ${DADA2_DIR}/asv_stats_table.tsv
awk 'BEGIN { ORS = "\t" } {print $1 $2 $14}' asv_stats_table.tsv


# create new header: ID\Total\tA3\tMock\tNH1\tNH4\tSim17\tSim22\tTh16\tTh38\tTh40\tX17007\tKingdom\tDomain\tPhyllum\tClass\tOrder\tFamily\tGenus\tSpecies\tPercent_similarity\tSequence
# sed '1i ID\tSample\tTaxonomy\tPercent_similarity\tSequence' > ${DADA2_DIR}/asv_stats.tsv



# # lets add sample ID
# # second field: Sample ID
# # removing the first two lines that contained only metadata
# # adding header
# rm -f ${DADA2_DIR}/asv_stats.tsv
# awk '{print $1"\t""total""\t"$3"\t"$4"\t"$2}' ${DADA2_DIR}/metadata/metadata.tsv |  \
# tail -n +3 | \
# sed '1i ID\tSample\tTaxonomy\tPercent_similarity\tSequence' > ${DADA2_DIR}/asv_stats.tsv

# the table still lack frequency per sample, so let's added it from feature (ASV) abundance table:
# join --header -1 1 -2 1 -t "\t"


# for SAMPLE in ${SAMPLES}
# do
#   grep '>' ${ASV_DIR}/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.fasta | \
#   sed 's/>//' | \ # first field: ASV ID
#   awk -v sample="${SAMPLE}" '{print $0"\t"sample}' # second field: Sample ID
# done



#TODO: copy asvs IDs as field 1
#TODO: copy asvs SampleID as field 2
#TODO: copy biome type as field 3
#TODO: copy abundance as field 4
#TODO: copy num_reads as field 5
#TODO: copy taxonomy as field 6
#TODO: copy Percent similarity as field 6
#TODO: copy sequences as field 7
