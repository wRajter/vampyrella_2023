# Functions used in tax_assignment_downstream Jupyter notebook

# Imports:
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def standardize_taxopath(taxopath):
    parts = taxopath.split(';')
    # Ensure there are 4 parts (order, family, genus, species)
    while len(parts) < 4:
        parts.append('unassigned')
    return ';'.join(parts)


def create_count_plot_all_samples(dataframe, tax_level, save_dir, xlabels_rotation=0, xlabels_ha='center'):
    '''
    Create and save a count plot for the specified taxonomic level, excluding 'Mock' samples.

    This function takes a pandas DataFrame, the taxonomic level to plot (e.g., 'family'), and a
    directory path where the plot will be saved. It filters out any samples labeled as 'Mock',
    creates a count plot for the remaining samples at the specified taxonomic level, and saves
    the plot as a PNG file.

    Parameters:
    - dataframe (pd.DataFrame): The DataFrame containing the taxonomic data.
    - tax_level (str): The taxonomic level for which to create the plot (e.g., 'family', 'genus').
    - save_dir (str): The directory path where the plot image will be saved.

    The function does not return any value. It saves the plot to the specified directory with
    a file name that reflects the taxonomic level and notes that 'Mock' samples are excluded.

    The plot is styled with a 'white' background, colorblind-friendly palette, and horizontal
    gridlines for readability. Bar labels are adjusted to replace specific text, and bar width
    is set for visual clarity.
    '''

    # Create path to the save directory if not exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Filter out the mock sample
    tax_assign_filtered = dataframe[dataframe['sample'] != 'Mock']

    # Set the style of the visualization
    sns.set_style('white')

    # Create a count plot
    plt.figure(figsize=(10, 6))
    bar_width = 0.7
    ax = sns.countplot(x=tax_level,
                    data=tax_assign_filtered,
                    order = tax_assign_filtered[tax_level].value_counts().index,
                    hue = tax_level,
                    legend = False,
                    palette='colorblind',
                    edgecolor='.2',
                    linewidth=0.7,
                    width=bar_width)

    # Retrieve current x-axis label texts
    labels = [label.get_text() for label in ax.get_xticklabels()]

    # Set the new labels to the x-axis
    ax.set_xticklabels(labels, rotation=xlabels_rotation, ha=xlabels_ha, fontsize=10)

    # Add horizontal gridlines for better readability
    plt.grid(axis='y', linestyle='-', alpha=0.7)

    sns.despine()

    # Optional: Set plot title and labels
    plt.title(f'Count Plot of Sequences at {tax_level.capitalize()} Level')
    plt.ylabel('Number of OTUs')
    plt.xlabel('')

    # Customizing the x-axis to have ticks from 0 to max value with an interval of 1
    max_count = tax_assign_filtered[tax_level].value_counts().max()
    plt.yticks(range(0, max_count + 1, 1))

    plt.savefig(os.path.join(save_dir, f'count_plot_otus_at_{tax_level}_level_all_samples_without_mock.png'), bbox_inches='tight', dpi=600)


def rename_samples(ax, rename_dict, xlabels_rotation=30, xlabels_ha='right', xlabels_fontsize=9):
    '''
    Rename sample names in the x-axis of a plot.

    Parameters:
    - ax: The axis object of the plot.
    - rename_dict: A dictionary mapping old sample names to new sample names.

    Returns:
    - ax: The updated axis object with renamed x-axis labels.
    '''
    # Retrieve current x-axis label texts
    labels = [label.get_text() for label in ax.get_xticklabels()]
    # Replace labels based on provided dictionary
    new_labels = [rename_dict.get(label, label) for label in labels]
    # Set the new labels to the x-axis
    ax.set_xticklabels(new_labels, rotation=xlabels_rotation, ha=xlabels_ha, fontsize=xlabels_fontsize)
    return ax


def tax_composition_per_sample(pivot_table,
                               taxonomic_level,
                               rename_dict,
                               save_dir,
                               xlabels_rotation=30,
                               xlabels_ha='right',
                               xlabels_fontsize=9):
    '''
    Create and save a stacked bar plot of taxonomic composition per sample.
    '''

    # Define a color palette with enough colors
    colors = sns.color_palette('colorblind', n_colors=pivot_table.shape[1])

    # Plotting
    plt.figure(figsize=(10, 6))
    ax = pivot_table.plot(kind='bar', stacked=True, color=colors, edgecolor='black', linewidth=0.7)

    # Styling
    plt.title(f'Taxonomic Composition Per Sample at {taxonomic_level.capitalize()} Level')
    plt.xlabel('')
    plt.ylabel('Number of OTUs')

    # Customizing the x-axis to have ticks from 0 to max value with an interval of 1
    # max_count = dataframe[taxonomic_level].value_counts().max()
    max_count = pivot_table.sum(axis=1).max()
    plt.yticks(range(0, int(max_count) + 1, 1))

    plt.legend(loc='upper right', frameon=False, title=taxonomic_level.capitalize(), fontsize='small')
    sns.despine()
    plt.tight_layout()

    # Rename x-axis labels as per the dictionary provided
    rename_samples(ax, rename_dict, xlabels_rotation=xlabels_rotation, xlabels_ha=xlabels_ha, xlabels_fontsize=xlabels_fontsize)

    # Check if save_dir exists and create if not
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Save the plot
    plt.savefig(os.path.join(save_dir, f'tax_comp_at_{taxonomic_level}_level_per_sample.png'), bbox_inches='tight', dpi=600)


def tax_composition_per_sample_percentages(pivot_table,
                                           taxonomic_level,
                                           rename_dict,
                                           save_dir,
                                           xlabels_rotation=30,
                                           xlabels_ha='right',
                                           xlabels_fontsize=9):
    '''
    Create and save a stacked bar plot of taxonomic composition per sample as percentages.
    '''
    # Define a color palette with enough colors
    colors = sns.color_palette('colorblind', n_colors=pivot_table.shape[1])

    # Normalize the pivot_table by row to get percentages
    pivot_table_percent = pivot_table.div(pivot_table.sum(axis=1), axis=0) * 100

    # Plotting
    plt.figure(figsize=(10, 6))
    ax = pivot_table_percent.plot(kind='bar', stacked=True, color=colors, edgecolor='black', linewidth=0.7)

    # Styling
    plt.title(f'Taxonomic Composition Per Sample at {taxonomic_level.capitalize()} Level (Percentages)')
    plt.xlabel('')
    plt.ylabel('Percentage of OTUs (%)')

    # Customizing the y-axis to have ticks at fixed intervals
    plt.yticks(range(0, 101, 10))

    plt.legend(loc='upper right', frameon=False, title=taxonomic_level.capitalize(), fontsize='small')
    sns.despine()
    plt.tight_layout()

    # Rename x-axis labels as per the dictionary provided
    rename_samples(ax, rename_dict, xlabels_rotation=xlabels_rotation, xlabels_ha=xlabels_ha, xlabels_fontsize=xlabels_fontsize)

    # Adjust the legend position
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False, title=taxonomic_level.capitalize(), fontsize='small')

    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    # Check if save_dir exists and create if not
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Save the plot
    plt.savefig(os.path.join(save_dir, f'tax_comp_at_{taxonomic_level}_level_per_sample_percent.png'), bbox_inches='tight', dpi=600)
