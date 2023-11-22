# Functions used in tax_assignment_downstream Jupyter notebook

# Imports:
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


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



def create_venn(save_dir, plot_data, project):

    '''
    Create and save a series of Venn diagrams comparing OTUs between two approaches.

    This function generates a 2x2 grid of Venn diagrams, each representing the overlap of OTUs
    between 'GAPPA' and 'VSEARCH' at different taxonomic levels (order, family, genus, species).
    The generated figure is saved in the specified directory with a filename based on the given project name.

    Parameters:
    save_dir (str): The directory path where the plot image will be saved. If the directory
                    does not exist, it will be created.
    plot_data (dict): A dictionary with keys in the format 'otus_{tax_level}_{approach}' where
                      'tax_level' is one of 'order', 'family', 'genus', or 'species', and
                      'approach' is either 'gappa' or 'vsearch'. Each key maps to a set of OTUs.
    project (str): A project identifier that will be included in the filename of the saved plot.
    '''

    # Create a figure with subplots
    fig, axs = plt.subplots(2, 2, figsize=(8, 6))  # Adjust the size as needed

    # Add a main title to the figure
    fig.suptitle('Venn Diagram of OTUs between Gappa and VSEARCH', fontsize=14)

    # Define colors for the sets
    set_colors = ['#d95f02', '#7570b3']
    labels = ['GAPPA', 'VSEARCH']

    # Define a function to increase text size in the Venn diagrams
    def increase_text_size(venn):
        for text in venn.set_labels:
            if text: text.set_fontsize(0)  # Increase set label size
        for text in venn.subset_labels:
            if text: text.set_fontsize(14)  # Increase subset label size

    # Plot A: OTUs
    v = venn2([plot_data['otus_order_gappa'], plot_data['otus_order_vsearch']], set_labels=(None, None), ax=axs[0, 0], set_colors=set_colors)
    increase_text_size(v)
    axs[0, 0].set_title('A) Order-level', loc='left')

    # Plot B: Family-level
    v = venn2([plot_data['otus_family_gappa'], plot_data['otus_family_vsearch']], set_labels=(None, None), ax=axs[0, 1], set_colors=set_colors)
    increase_text_size(v)
    axs[0, 1].set_title('B) Family-level', loc='left')

    # Plot C: Genus-level
    v = venn2([plot_data['otus_genus_gappa'], plot_data['otus_genus_vsearch']], set_labels=(None, None), ax=axs[1, 0], set_colors=set_colors)
    increase_text_size(v)
    axs[1, 0].set_title('C) Genus-level', loc='left')

    # Plot D: Species-level
    v = venn2([plot_data['otus_species_gappa'], plot_data['otus_species_vsearch']], set_labels=(None, None), ax=axs[1, 1], set_colors=set_colors)
    increase_text_size(v)
    axs[1, 1].set_title('D) Species-level', loc='left')

    from matplotlib.patches import Patch
    legend_handles = [Patch(facecolor=set_colors[i], edgecolor=None, label=labels[i]) for i in range(len(labels))]
    fig.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=len(labels), frameon=False, fontsize=12)

    # Adjust layout for better spacing with the main title
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust the rect parameter as needed for your figure

    # Show the figure
    # plt.show()

    # Check if save_dir exists and create if not
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Save the plot
    plt.savefig(os.path.join(save_dir, f'venn_comp_gappa_and_vsearch_{project}.png'), bbox_inches='tight', dpi=600)



def count_tax_assigned_otus(melted_df, save_dir):


    '''
    Create and save a count plot for the number of taxonomically assigned OTUs.

    This function takes a DataFrame with taxonomic assignments at different levels and approaches,
    and creates a count plot comparing the number of OTUs assigned in each category. The plot
    is then saved to the specified directory.

    Parameters:
    melted_df (pd.DataFrame): A pandas DataFrame containing the melted data with columns
                              'taxonomic_level' and 'approach'.
    save_dir (str): The directory path where the plot image will be saved. If the directory
                    does not exist, it will be created.

    The DataFrame should have the following columns:
    - 'taxonomic_level': The level of taxonomy (order, family, genus, species).
    - 'approach': The approach used for taxonomic assignment ('GAPPA' or 'VSEARCH').
    - 'status': The assignment status ('Assigned' or 'Unassigned').

    The function will create a count plot with 'taxonomic_level' on the x-axis, the count of OTUs on
    the y-axis, and different hues representing each approach. The plot is saved with a filename
    'count_assigned_otus_gappa_vd_vsearch.png' at a resolution of 600 dpi.

    '''

    # Set the style of the visualization
    sns.set_style('white')

    # Now create the count plot
    plt.figure(figsize=(10, 6))
    ax = sns.countplot(data=melted_df,
                    x='taxonomic_level',
                    hue='approach',
                    palette=['#d95f02', '#7570b3'],
                    edgecolor = '.2',
                    dodge=True,
                    width= 0.7,
                    linewidth=0.7)

    # Add horizontal gridlines for better readability
    plt.grid(axis='y', linestyle='-', alpha=0.7)

    sns.despine()

    # Enhance the plot
    plt.title('Number of assigned OTUs (GAPPA vs. VSEARCH)')
    plt.ylabel('Count of taxonomically assigned OTUs')
    plt.xlabel('')

    # Adjust the subplot params to give some space for the external legend
    plt.subplots_adjust(right=0.75)

    # Adding legend
    plt.legend(title='Approach', loc='upper left', bbox_to_anchor=(1, 1), frameon=False, fontsize='medium')


    plt.tight_layout(rect=[0, 0, 0.85, 1])

    # plt.show()

    # Create path to the save directory if not exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    plt.savefig(os.path.join(save_dir, f'count_assigned_otus_gappa_vd_vsearch.png'), bbox_inches='tight', dpi=600)
