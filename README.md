# Vampyrellid diversity exploration 2023

## Project description

This project aims to explore Vampyrellid's diversity using long-read metabarcoding of the eukaryotic rDNA operon.
In total, we collected nine environmental samples and one mock sample (Table 1). All environmental DNA has been extraced using the Qiagen's PowerSoil Powerlyser Kit, and sent to Integrated Microbiome Resource (IMR) for sequencing. IMR sequenced the full 18S and ITS rRNA amplicon fragments by PacBio Sequel 2.
The details on librabry preparation, and sequencing are described on the IMR website: https://imr.bio/protocols.html.

## Initial data

### Raw reads

We received demultiplexed raw reads in fastq format with phred+33 encoding from IMR.
Files: /raw_data/PacBio/{SuthausFull18S, SuthausFullITS}/{cell1, cell2}
The reads were from full 18S and ITS rRNA amplicon fragments. For each sample, we received two version of the reads - cell1 and cell2. Cell1 samples were sequenced with 1x depth and cell2 with 2x depth.

### Reference alignments and trees

For **taxonomic assignment** of the environmental sequences, we used full PR2 reference alignment (version 4.14.0). We manually updated the PR2 alignment with TODO:number vampyrellid sequences (TODO: create a list of seqs - names and GB ids) to have more accurate and refine taxonomic assignment of vampyrellids.
For **phylogenetic assignment**, we... TODO

## Project folder hierarchy

**analysis:** the main output files from the analyses are stored here. (Not present at GitHub for storage size and possession purposes). Please, contact me if you want to access to any of these data.
**notebooks:** all the jupyter notebooks are stored here.
**raw_data:** the initial data and intermediate results are stored here. (Not present at GitHub for storage size and possession purposes). Please, contact me if you want to access to any of these data.

- fastqc_out: output from the quality assessment of the raw reads
- PacBio: raw reads in fastq format received from IMR
- packages: all the programms installed from binaries
- qiime_input: input data (mostly manifest files) needed for qiime2 analyses
- qiime_output: all the intermediate output files from qiime2 analyses (artifact files, tax and assignment intermediate results)
- reference alignments: PR2 databases and reference alignments for phylogenetic placement
  **scripts:** the code written during this project is stored here:
- bash_scripts: script that are executed by the bash program
- python_functions: all the python files
- r_scripts: all the R files
  **test:** directory for testing the code.

## Chronology of the analyses

### Inspecting read quality

Using inspect_reads_quality.sh script to generate quality reports employing fastqc and multiqc.
_Input:_ fastq files from raw reads.
_Path:_ /rawdata/Pacbio/
_Output:_ quality reports in html format for each fastq file and one aggregated report(multiqc_report).
_Path:_ /raw_data/fastqc_data/

### Importing raw reads (FASTQ files) into qiime2 artifact files

For importing the raw reads into qiime2 artifacts, I had to prepare manifest files. Those represent tsv files with a specific header, sample names, and the full paths for the raw sequences (review qiime2 manual https://docs.qiime2.org/2022.11/tutorials/importing/ for more details).
raw*data/qiime_input/PacBioCCSmanifest*\[molecular*marker\]*\[cell\].tsv
Then, I used qiime*reads2asv.sh script to import the reads to qiime2 artifact format and generate reads summary file:
*Input:* manifest files.
*Path:* /raw_data/qiime_input/PacBioCCSmanifest*\[molecular*marker\]*\[cell\].tsv
_Output:_ qiime formated reads and the summary file.
_Path:_ /raw*data/qiime_output/reads_qza/raw_reads*\[molecular*marker\]*\[cell\].qza
_Path:_ /raw*data/qiime_output/reads_qza/raw_reads_summary*\[molecular*marker\]*\[cell\].qzv

### Denoising reads and creating ASVs using DADA2 via qiime2

We used qiime denoiser to remove primers, filter sequences based on the length, correct reads and get amplicon sequence variants (ASVs). The denoising step is part of the qiime*reads2asv.sh script.
*Input:* raw reads in qiime2 format.
*Path:* /raw_data/qiime_output/reads_qza/raw_reads*\[molecular*marker\]*\[cell\].qza
_Output._ Feature table (similar to abundance ASV table), ASVs in qiime2 format, and stats file
_Path:_ /raw*data/qiime_output/dada2_output/table*\[molecular*marker\]*\[cell\].qza
_Path:_ /raw*data/qiime_output/dada2_output/representative_sequences*\[molecular*marker\]*\[cell\]_ALLSAMPLES.qza
*Path:* /raw_data/qiime_output/dada2_output/stats_\[molecular*marker\]*\[cell\].qza

### Post denoising summary statistics and visualizations

Afer the denoising the raw reads we summarized DADA2 output, visualized the representative sequences, and converted ASV sequences from qza to fasta format. All these analyses are part of the qiime_reads2asv.sh script.

1. Summarizing DADA2 output:
   _Input:_ qiime2 feature table.
   _Path:_ /raw*data/qiime_output/dada2_output/table*\[molecular*marker\]*\[cell\].qza
   _Output:_ qiime2 graphic visualization file (qzv).
   _Path:_ /raw*data/qiime_output/dada2_output/dada2_table_summary*\[molecular*marker\]*\[cell\].qzv
2. Visualization of the representative sequences:
   _Input:_ ASVs in qiime2 format.
   _Path:_ /raw*data/qiime_output/dada2_output/representative_sequences*\[molecular*marker\]*\[cell\]_ALLSAMPLES.qza
   *Output:* qiime2 graphic visualization file (qzv).
   *Path:* /raw_data/qiime_output/dada2_output/representative_sequences_${MARKER}_${CELL}\_ALLSAMPLES.qzv
3. Converting ASV sequences from qza to fasta format.
   _Input:_ ASVs in qiime2 format.
   _Path:_ /raw*data/qiime_output/dada2_output/representative_sequences*\[molecular*marker\]*\[cell\]_ALLSAMPLES.qza
   *Output:* fasta files with the ASVs
   *Path:* /raw_data/qiime_output/ASV/representative_sequences_\[molecular*marker\]*\[cell\]\_ALLSAMPLES.fasta

### Filtering ASVs sequences into individual samples

In this step, I split the ASVs based on the sample they belong to, and then created a fasta file for each sample.
To best of my knowledge, this step is not trivial and there is no a single commnand in the qiime2 environmnet to separate the ASV based on the samples and then generate fasta files. Nevertheless, there are several ways how to tackle it.
I transposed the feature table to flip the sample and feature axes, so the x-axis contained samples and y-axis encomapassed ASVs. In this way, it was possible to use the transposed table as feature metadata, and filter ASVs based on the samples they belong to using the SQL where clause syntax. The ASVs were at first filtered into the qiime2 format and then into the fasta files. All the steps are part of the qiime_reads2asv.sh script.

1. Transposing the qiime2 ASV table in the qza and qzv formats.
   _Input:_ qiime2 ASV table.
   _Path:_ raw*data/qiime_output/dada2_output/table*\[molecular*marker\]*\[cell\].qza
   _Output:_ transposed qiime2 ASV table.
   _Path:_ raw*data/qiime_output/dada2_output/transposed_table*\[molecular*marker\]*\[cell\].qza
   _Path:_ raw*data/qiime_output/dada2_output/transposed_table*\[molecular*marker\]*\[cell\].qzv
2. Filtering ASVs into the individual samples using the transposed ASV table.
   _Input:_ ASVs in qiime2 format and metadata file in form of the transposed table.
   _Path:_ /raw*data/qiime_output/dada2_output/representative_sequences*\[molecular*marker\]*\[cell\]_ALLSAMPLES.qza
   *Path:* raw_data/qiime_output/dada2_output/transposed_table_\[molecular*marker\]*\[cell\].qza
   _Output:_ ASVs in qiime2 format for each sample in the qza and qzv formats.
   _Path:_ raw*data/qiime_output/dada2_output/representative_sequences*\[molecular*marker\]*\[cell\]_\[sample\].qza
   *Path:* raw_data/qiime_output/dada2_output/representative_sequences_\[molecular*marker\]*\[cell\]\_\[sample\].qzv
3. Converting ASV sequences from qza to fasta format.
   _Input:_ ASVs in qiime2 format for each sample.
   _Path:_ raw*data/qiime_output/dada2_output/representative_sequences*\[molecular*marker\]*\[cell\]_\[sample\].qza
   *Output:* ASVs in fasta files for each sample.
   *Path:* raw_data/qiime_output/ASV/representative_sequences\_\_\[molecular_marker\]_\[cell\]\_\[sample\].fasta

### Converting reference sequences and taxonomy into the qiime2 artifact files

TODO

### Taxonomic assignment using v-search through qiime2

As taxonomic assignment is computationally intensive, it was done at CHEOPS cluster (https://rrzk.uni-koeln.de/hpc-projekte/hpc). I used full PR2 reference alignment (version 4.14.0) enriched by vampirellid sequences (TODO: add list of the seqs names and GB ids). I run two scripts at the cluster: one for all samples combined (_qiime_taxassign.slurm_), and the other one that created taxonomic assignemnt for each sample separatelly (_qiime_taxassign_persample.slurm_).
_Input:_ ASVs in qiime2 format, reference sequences and taxonomy in qiime2 format.
_Path:_ raw*data/qiime_output/dada2_output/representative_sequences*\[molecular*marker\]*\[cell\]_\[sample\].qza
*Path:* raw_data/reference_alignments/ref_sequences.qza
*Path:* raw_data/reference_alignments/ref_taxonomy.qza
*Output:* vsearch taxonomy file in qiime2 format.
*Path:* raw_data/qiime_output/assignment_results/\[molecular_marker\]_\[cell\]/vsearch*taxonomy*\[molecular*marker\]*\[cell\]\_\[sample\].qza

### Post-taxonomic assignment analyses

After the taxonomic assignment was done, I used the _qiime_post_taxassign.sh_ script to convert the output file into qiime visualization format and the taxonomic assignemnt tsv table.
_Input:_ vsearch taxonomy file in qiime2 format.
_Path:_ raw*data/qiime_output/assignment_results/\[molecular_marker\]*\[cell\]/vsearch*taxonomy*\[molecular*marker\]*\[cell\]_\[sample\].qza
*Output:* taxonomic assignment qiime2 visualization file and taxonomic assignment tsv table
*Path:* raw_data/qiime_output/assignment_results/\[molecular_marker\]_\[cell\]/tax*viz*\[sample\].qzv
_Path:_ raw*data/qiime_output/assignment_results/\[molecular_marker\]*\[cell\]/taxonomy\_\[sample\].tsv

### Filtering Vampyrella-specific ASVs

I filtered ASVs that were assigned into Vampyrellida and created separate fasta files using the _artefact2fasta.sh_ script.
_Input:_ ASVs in fasta file and taxonomic assignment tsv table.
_Path:_ raw*data/qiime_output/ASV/representative_sequences*\[molecular*marker\]*\[cell\]_\[sample\].fasta
*Path:* raw_data/qiime_output/assignment_results/\[molecular_marker\]_\[cell\]/taxonomy*\[sample\].tsv
*Output:* fasta files with vampyrella-specific ASVs for each sample.
*Path:* raw_data/qiime_output/dada2_output/ASVs/Vampyrellida*\[molecular*marker\]*\[cell\]\_\[sample\].fasta

### Creating an ASV summary table

_Input:_
_Path:_
_Output:_
_Path:_
