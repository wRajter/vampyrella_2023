# tax_assign_denoise.jl

using Pkg

Pkg.add(url="https://github.com/MurrellGroup/NextGenSeqUtils.jl")
Pkg.add(url="https://github.com/MurrellGroup/DPMeansClustering.jl")
Pkg.add(url="https://github.com/MurrellGroup/RobustAmpliconDenoising.jl")

using NextGenSeqUtils, RobustAmpliconDenoising

# Accept user input for the file path

filt_seqs = "../../raw_data/PacBio/Suthaus_2022_Full18S/cellCombined/filtered/Mock_18S.hifi_reads.fastq"

seqs, QVs, seq_names = read_fastq(filt_seqs)

templates, template_sizes, template_indices = denoise(seqs)

save_fasta = "Mock_21R_asv.fasta"

write_fasta(save_fasta, templates, names = ["seq$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])
