# Vampyrellid diversity exploration 2023

## Project description

This project aims to explore Vampyrellid's diversity using long-read metabarcoding of the eukaryotic rDNA operon.
In total, we collected nine environmental samples and one mock sample (Table 1). All environmental DNA has been extraced using the Qiagen's PowerSoil Powerlyser Kit, and sent to Integrated Microbiome Resource (IMR) for sequencing. IMR sequenced the full 18S and ITS rRNA amplicon fragments by PacBio Sequel 2.
The details on librabry preparation, and sequencing are described on the IMR website: https://imr.bio/protocols.html.

## Initial data

### Raw reads

We received demultiplexed raw reads in fastq format with phred+33 encoding from IMR. \
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

In this step, we generate quality reports from the raw PacBio reads. Input files represent `fastq` files from raw reads. The quality reports are genereted using FastQC (TODO: add version) and MultiQC (TODO: add version) packages. The output files (quality reports) are in HTML format - one output file per each `fastq` input file and one aggregated report for all samples (multiqc report).
Script: `reads_quality.sh`
Output directory: raw_data/fastqc_data/

### Converting raw reads (in FASTQ format) into Qiime2 artifact files

For importing the raw reads into Qiime2 artifacts, we need the raw reads in `fastq` format and so-called manifest file. The manifest file has to contain a specific header and two columns: i) sample name and ii) full path for raw reades (review qiime2 manual https://docs.qiime2.org/2022.11/tutorials/importing/ for more details). After the manifest file is prepared, Qiime2 (version 2022.11) is used to convert the reads to Qiime2 artifact format (`qza` extension) and generate reads summary file (`qzv` extension).
Script: `reads2asv.sh`
Manifest files: raw_data/qiime_input/
Output directory: raw_data/qiime_output/reads_qza/

### Denoising reads and creating ASVs using DADA2 (via Qiime2)

The Qiime2's DADA2 denoiser is used to correct reads and get amplicon sequence variants (ASVs). This step also removes primers and filters sequences based on the length. The input files represent raw reads in Qiime2 format that were created in the previous step. The output represent three files: 1) abundance ASV table (also called a feature table), 2) ASVs with corresponding sequences in Qiime2 `qza` format (similar to fasta format), and 3) summary statistic file.
Script: `reads2asv.sh`
Output directory: raw_data/qiime_output/dada2_denoise/

### Post denoising summary statistics and visualizations

Afer the denoising the raw reads into ASVs, we make the following steps:

1. Summarizing DADA2 output using the Qiime2 feature table from the denoising step.
2. Visualizing the representative (ASV) sequences, which basicly means to convert the representative sequence `qza` files to `qzv` files.
3. Converting ASV sequences from the `qza` to fasta format.

Script: `qiime_reads2asv.sh`
Output direcoty: raw_data/qiime_output/dada2_post_denoise/

### Filtering ASVs sequences into individual samples

In this step, the ASVs are split based on the sample they belong to. Then, a fasta file and an abundance table is created for each sample.
This step is not trivial. To best of my knowledge, there is no a single commnand in the Qiime2 environmnet to separate the ASV based on the samples and then generate fasta files. Nevertheless, there are several ways how to tackle it.
First, the feature table is transposed flipping the sample and feature axes. After transposing the table, the x-axis contained samples and y-axis contained ASVs. In this way, it is possible to use the transposed table as feature metadata, and filter ASVs based on the samples they belong to using the SQL `where` clause syntax. The ASVs are at first filtered into the Qiime2 format and then into the fasta files. The procedure is summarized by the following steps:

1. Transposing the qiime2 ASV table in the qza and qzv formats.
2. Filtering ASVs into the individual samples using the transposed ASV table.
3. Converting ASV sequences from qza to fasta format.

Script: `qiime_reads2asv.sh`
Output direcoty: raw_data/qiime_output/dada2_post_denoise/

### Converting reference sequences and taxonomy into the qiime2 artifact files

TODO: add this part

### Taxonomic assignment using v-search through qiime2

As taxonomic assignment is computationally intensive, it is done at CHEOPS cluster (https://rrzk.uni-koeln.de/hpc-projekte/hpc). We use the full PR2 reference alignment (version 4.14.0) enriched by vampyrellid sequences (TODO: add list of the seqs names and GB ids). We run two scripts at the cluster: one for all samples combined (`qiime_taxassign.slurm`), and the other that creates taxonomic assignemnt for each sample separatelly (`qiime_taxassign_persample.slurm`). The input files for the taxonomic assignement are ASV sequences, reference sequences and reference taxonomy, all in Qiime2 (`gza`) format. The output is vsearch taxonomy file in Qiime2 format.
Scripts: `qiime_taxassign.slurm` and `qiime_taxassign_persample.slurm`
Output directory: raw_data/qiime_output/assignment_results/

### Post-taxonomic assignment

After the taxonomic assignment is done, we convert the output file into Qiime2 visualization format (`gzv`) and the taxonomic assignemnt tsv table.
Output directory: raw_data/qiime_output/assignment_results/
Script: `qiime_post_taxassign.sh`

### Filtering Vampyrella-specific ASVs

In this step, we filter sequences that were assigned into Vampyrellida. The input files are ASVs in fasta files and taxonomic assignment tsv table. The vampyrellid sequences are outputed in form of fasta files for all samples and per each sample separetelly.
Script: `filter_vamp_seqs.sh`
Output directory: raw_data/qiime_output/dada2_post_denoise/fasta

### Creating an ASV summary table

As the ASV taxonomic assignemnt part is finished, we summarizes all the output data we have yielded so far into a final table. This table contains the following information: ASV ID, abundance per sample, taxonomic assignemnt at various taxonomic level, percentage similarity of the assignment, and the sequence for each ASV. As those pieces of information are scattered in several output files, we created a custome script to process the files and assembled the information in a summary tsv table (asv_summary_table). This table is used for data visualization in the following steps.
Script: `create_summary_table.sh`
Output directory: raw_data/qiime_output/dada2_post_denoise/
