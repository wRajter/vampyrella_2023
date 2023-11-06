# tax_assign_denoise.jl

using Pkg

# Pkg.add(url="https://github.com/MurrellGroup/NextGenSeqUtils.jl")
# Pkg.add(url="https://github.com/MurrellGroup/DPMeansClustering.jl")
# Pkg.add(url="https://github.com/MurrellGroup/RobustAmpliconDenoising.jl")

# Load required packages
using NextGenSeqUtils, RobustAmpliconDenoising, CodecZlib

# Define project information
project = "Suthaus_2022"
marker = "Full18S"
denoised_method = "RAD"
raw_data = "../../raw_data"
suffix = ".fastq"

# Construct the raw_reads_dir path
reads_dir = joinpath(raw_data, "dada2", project, marker, "filtered")

# List files in the raw_reads_dir directory
files = [file for file in readdir(reads_dir) if endswith(file, suffix)]


# Loop through each sample and perform denoising
for file in files

  # Construct the file path for the current sample
  filt_seq_path = joinpath(reads_dir, file)

  sample_name = match(r"^(.*?)\.(.*)", file).captures[1]

  # Read and denoise the sequences
  seqs, QVs, seq_names = read_fastq(filt_seq_path)
  templates, template_sizes, template_indices = denoise(seqs)

  # Construct the save_fasta path for the denoised sequences
  save_fasta = joinpath(raw_data, "denoised", project, marker, denoised_method, "$sample_name" * "_asv.fasta")

  # Write denoised sequences to a fasta file
  write_fasta(save_fasta, templates, names = ["seq$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])

  # Print a message indicating completion for the current sample
  println("Denoising for $sample_name done. Fasta file save to $save_fasta")

end
