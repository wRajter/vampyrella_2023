# tax_assign_denoise.jl

# using Pkg

# Pkg.add(url="https://github.com/MurrellGroup/NextGenSeqUtils.jl")
# Pkg.add(url="https://github.com/MurrellGroup/DPMeansClustering.jl")
# Pkg.add(url="https://github.com/MurrellGroup/RobustAmpliconDenoising.jl")

# Load required packages
using NextGenSeqUtils, RobustAmpliconDenoising

# Define project information
project = "Jamy_2022"
marker = "rDNA"
cell = "cell"
raw_data = "../../raw_data"
suffix = ".fastq"

# Construct the raw_reads_dir path
raw_reads_dir = joinpath(raw_data, "PacBio", "$project" * "_" * "$marker", cell, "filtered")

# List files in the raw_reads_dir directory
files = readdir(raw_reads_dir)


# Extract sample names from filenames using the specified suffix
samples = [match(r"^(.*?)\.(.*)", file).captures[1] for file in files if endswith(file, suffix)]


# Loop through each sample and perform denoising
for sample in samples

    # Construct the file path for the current sample
    filt_seqs = "$raw_reads_dir/$sample$suffix"

    # Read and denoise the sequences
    seqs, QVs, seq_names = read_fastq(filt_seqs)
    templates, template_sizes, template_indices = denoise(seqs)

     # Construct the save_fasta path for the denoised sequences
    save_fasta = joinpath(raw_data, "denoise", project, marker, cell, "FAD", "$sample" * "_asv.fasta")

    # Write denoised sequences to a fasta file
    write_fasta(save_fasta, templates, names = ["seq$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])

    # Print a message indicating completion for the current sample
    println("Denoising for $sample done. Fasta file save to $save_fasta")

end
