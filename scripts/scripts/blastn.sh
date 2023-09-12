#!/bin/bash

# Path to the BLASTN tool
BLASTN_PATH="/path/to/blastn"

# Path to your database
DB_PATH="/path/to/your/blast/database/nt"

# Directory containing your input sequences
INPUT_DIR="/path/to/your/input/sequences"

# Directory to save BLAST results
OUTPUT_DIR="/path/to/save/output"

for file in $INPUT_DIR/*.fasta; do
    base=$(basename $file .fasta)
    $BLASTN_PATH -query $file -db $DB_PATH -out $OUTPUT_DIR/${base}_blast_results.txt -outfmt 6 -max_target_seqs 5
done
