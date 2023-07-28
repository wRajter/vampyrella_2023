# Creating stats for tracking the numbor of reads after the individual filtering steps in the pipeline

# Imports
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Variables
project = 'Suthaus_2022'
cell = 'cell'
marker = 'rDNA'
sim = 'sim97'
raw_data = '../../raw_data'
suffix = '.hifi_reads.fastq.gz'
denoise_path = f'{raw_data}/denoise/{project}/{marker}/{cell}/dada2'
results_path = f'{raw_data}/figures_and_tables_generated/{project}/{marker}/{cell}'

# Script

# Stats for css, primers, and filtered sequences
stats1 = pd.read_table(f'{denoise_path}/reads_num_track.csv',
                       sep=" ",
                       index_col=0)
stats1.index = stats1.index.str.replace(suffix, '', regex=False)

# Stats for denoise sequences
all_files = os.listdir(denoise_path)
all_asv_files = [file for file in all_files if '.fasta' in file]
all_asv_paths = [denoise_path + '/' + file_name for file_name in all_asv_files]

dict_counter = {}

for i in all_asv_files:
    dict_counter[i] = None

dict_counter

for asv_path in all_asv_files:
    with open(denoise_path + '/' + asv_path, 'rt') as f:
        counter_ = 0
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                counter_ += 1
        dict_counter[asv_path] = counter_

stats2 = pd.DataFrame.from_dict(dict_counter,
                                orient='index',
                                columns=['denoised'])
stats2.index = stats2.index.str.replace('_asv.fasta', '', regex=False)

# Combine stats1 and stats2 tables
track = pd.concat([stats1, stats2], axis=1)

# Save final stats table as csv
if not os.path.exists(results_path):
    os.makedirs(results_path)
    print(f'\nThe new directory {results_path} is created!\n')

track.to_csv(f'{results_path}/reads_num_track.csv')
print(f'\nThe table "reads_num_track.csv" is saved to {results_path}\n')

# Plotting the reads count through filtering steps (per sample)

plotting_barplot1 = input(
    '\nDo you want to save the "Reads count through filtering steps (per sample)" plot? (y/n): '
)

if plotting_barplot1 == 'y':

    # Adjust data for plotting
    track_plot = track.reset_index()
    track_plot = track_plot.rename(columns={'index': 'sample'})
    track_plot = track_plot.melt(
        id_vars=['sample'],
        value_vars=['ccs', 'primers', 'filtered', 'denoised'])

    # Plotting stats per sample
    sns.set(rc={'figure.figsize': (11.7, 8.27)})
    sns.set_style("white")
    plot = sns.barplot(data=track_plot, y='variable', x="value", hue='sample')
    plot.set(ylabel='',
             xlabel='number of reads',
             title='Reads count through filtering steps (per sample)')
    sns.despine()

    # adjusting legend
    h, l = plot.get_legend_handles_labels()
    plot.legend(handles=h, title='Sample', loc='lower right', frameon=False)

    # save
    plt.savefig(f'{results_path}/reads_num_track_per_sample.png',
                dpi=300,
                transparent=False)

    print(
        f'\nThe plot reads_num_track_per_sample.png is saved to {results_path}\n'
    )

    # Clear the current figure to avoid any carryover of settings to the next plot
    plt.clf()

# Plotting the reads count through filtering steps (all sample)

plotting_barplot2 = input(
    '\nDo you want to save the "Reads count through filtering steps (all samples pooled)" plot? (y/n): '
)

if plotting_barplot2 == 'y':

    # Adjusting data for plotting
    track_sum = track.sum(axis=0)
    track_sum = track_sum.to_frame()
    track_sum = track_sum.reset_index()
    track_sum = track_sum.rename(columns={0: 'reads_num', 'index': 'step'})

    # Plotting stats, all samples
    sns.set(rc={'figure.figsize': (11.7, 8.27)})
    sns.set_style("white")
    plot = sns.barplot(data=track_sum, y='reads_num', x="step")
    plot.set(xlabel='',
             ylabel='number of reads',
             title='Reads count through filtering steps (all samples pooled)')
    sns.despine()

    # save
    plt.savefig(f'{results_path}/reads_num_track_all_sample.png',
                dpi=300,
                transparent=False)

    print(
        f'\nThe plot reads_num_track_all_sample.png is saved to {results_path}\n'
    )
