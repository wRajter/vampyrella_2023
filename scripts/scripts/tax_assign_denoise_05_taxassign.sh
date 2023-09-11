#!/bin/bash

# Taxonomic assignment using vsearch

# Variables
THREADS=12
PROJECT="Suthaus_2022"
MARKER="Full18S"
CELL="cellCombined"
SIM="sim95"
DENOISE_METHOD="RAD"
RAW_DATA="../../raw_data"
DATABASE="${RAW_DATA}/reference_alignments/pr2/pr2_version_4.14.0_SSU_UTAX_plus_review_vamp_2023.fasta"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"
SEQTK="${RAW_DATA}/packages/seqtk"
IDENTITY=0.6
# Input directory
QUERY_SEQ="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}/${DENOISE_METHOD}"
# Output directory
ASSIGNMENT_DIR="${RAW_DATA}/tax_assign_results/${PROJECT}/${MARKER}/${CELL}/${SIM}/${DENOISE_METHOD}"
VAMP_SEQ_DIR="${RAW_DATA}/vamp_specific_seqs/${PROJECT}/${MARKER}/${CELL}/${SIM}/${DENOISE_METHOD}"


################
## ASSIGNMENT ##
################


# samples
SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')



echo "Samples used:"
echo "$SAMPLES"

mkdir -p ${ASSIGNMENT_DIR}/

# Taxonomic assignment analysis:

for SAMPLE in ${SAMPLES}
do
  echo -e "\n\n\n###### Working on sample: $SAMPLE ######"
  rm -f ${ASSIGNMENT_DIR}/blast6_${SAMPLE}.tab
  vsearch --usearch_global ${QUERY_SEQ}/${SAMPLE}_otu.fasta \
         --dbmask none \
         --qmask none \
         --db ${DATABASE} \
         --id ${IDENTITY} \
         --iddef 3 \
         --blast6out ${ASSIGNMENT_DIR}/blast6_${SAMPLE}.tab
  echo -e "\nVampyrellids found:"
  grep 'Vamp' ${ASSIGNMENT_DIR}/blast6_${SAMPLE}.tab | wc -l
done

# if using the exrtacted 18S from the long fragment use this path
# --usearch_global ${QUERY_SEQ}/extracted_18S_${SAMPLE}.fasta
#  # --usearch_global ${OTU_CHIM_FILT}/otu_${SAMPLE}.fasta \

###############################
## PULLING OUT THE VAMP SEQS ##
###############################


mkdir -p ${VAMP_SEQ_DIR}

for SAMPLE in ${SAMPLES}
do
  rm -f ${ASSIGNMENT_DIR}/vamp_ids_${SAMPLE}.txt \
        ${VAMP_SEQ_DIR}/${SAMPLE}_otu.fasta

  grep 'Vampyrellida' ${ASSIGNMENT_DIR}/blast6_${SAMPLE}.tab | \
  awk '{print $1}' > \
  ${ASSIGNMENT_DIR}/vamp_ids_${SAMPLE}.txt

  ${SEQTK}/seqtk subseq \
    ${QUERY_SEQ}/${SAMPLE}_otu.fasta \
    ${ASSIGNMENT_DIR}/vamp_ids_${SAMPLE}.txt > \
    ${VAMP_SEQ_DIR}/${SAMPLE}_otu.fasta

  rm -f ${ASSIGNMENT_DIR}/vamp_ids_${SAMPLE}.txt
done



# #################################################################
# ## CREATING A SINGEL FASTA FILE WITH ALL SEQS FROM ALL SAMPLES ##
# #################################################################


# # Add the sample id to the each sequence in the fasta files
# for SAMPLE in ${SAMPLES}
# do
#   rm -f ${VAMP_SEQ_DIR}/otu_vamp_labeled_${SAMPLE}.fasta
#   sed "/;seqs=/a ;${SAMPLE}" ${VAMP_SEQ_DIR}/otu_vamp_${SAMPLE}.fasta | \
#   sed "/;seqs=/ { N; s/\n// }" > \
#   ${VAMP_SEQ_DIR}/otu_vamp_labeled_${SAMPLE}.fasta
# done

# # Combine all fasta files into a single fasta file
# cat ${VAMP_SEQ_DIR}/otu_vamp_labeled_* > ${VAMP_SEQ_DIR}/${PROJECT}_${MARKER}_${SIM}_all_seqs.fasta
# rm ${VAMP_SEQ_DIR}/otu_vamp_labeled_*
