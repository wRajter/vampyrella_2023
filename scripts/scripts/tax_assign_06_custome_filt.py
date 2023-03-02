# custome based filtering of rare OTUs per sample

# load modules
import qiime2
import pandas as pd
import math

# note: qiime2 is install in the qiime2-2022.11 conda environment

# filtering rare OTUs
threshold = 0.001
datasets_labels = [
    "cell1_sim97", "cell1_sim99", "cell2_sim97", "cell2_sim99",
    "cellCombined_sim97", "cellCombined_sim99"
]

for dataset in datasets_labels:
    # load data into dataframe
    cell, sim = dataset.split("_")
    raw_data = "../../raw_data"
    otu_table_fp = f"{raw_data}/OTU_clust/Full18S/{cell}/{sim}/otu_table.qza"
    otu_table_qiime = qiime2.Artifact.load(otu_table_fp)
    otu_table = otu_table_qiime.view(pd.DataFrame)
    # otu abundance (total observation) per sample
    abundance_per_sample = otu_table.sum(axis=1)
    # minimal OTU frequency threshold per sample
    sample_thresholds = {}
    for i in range(len(abundance_per_sample)):
        sample_thresholds[abundance_per_sample.index[i]] = math.floor(
            abundance_per_sample[i] * threshold)
    # changing abundance of otus below threshold to zero
    for sample in sample_thresholds:
        current_threshold = sample_thresholds[sample]
        otu_table[otu_table.loc[[sample], ] < current_threshold] = 0
    # transpose the table
    otu_table = otu_table.transpose()
    # save to tsv
    path = f"{raw_data}/OTU_filtered/Full18S/{cell}/{sim}"
    otu_table.to_csv(f"{path}/table_rarefilt.tsv", sep="\t")
