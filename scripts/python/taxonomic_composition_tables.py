import seaborn as sns
import matplotlib.pyplot as plt


def adjust_suthaus_2020_df(dataframe):
    '''
    Adjusts the dataframe based on the format of Suthaus (2020) to make it more interpretable and organized.

    This function performs multiple transformations:
    1. Removes the prefix 'centroid=' from the first column.
    2. Splits the first column into two separate columns ('OTU' and 'OTU_Num').
    3. Extracts a number representing the number of sequences from the 'OTU_Num' column.
    4. Drops the original first column.
    5. Extracts taxonomic hierarchy information (from 'Kingdom' to 'Species') from the second column.
    6. Drops the original second column.
    7. Extracts alignment information ('Pident' and 'Length') from subsequent columns.
    8. Rename the deep sea samples so they will appear next to each other on the plot.
    9. Rearranges the columns for better readability.

    Parameters:
    - dataframe : pd.DataFrame
        A dataframe based on the Suthaus (2020) format containing taxonomic and alignment data.

    Returns:
    - pd.DataFrame
        A transformed dataframe with extracted and organized columns.
    '''
    # Replace 'centroid=' with an empty string in the first column
    dataframe[dataframe.columns[0]] = dataframe[
        dataframe.columns[0]].str.replace('centroid=', '', regex=False)
    # Split the first column based on the ';' delimiter
    dataframe[['OTU', 'OTU_Num'
               ]] = dataframe[dataframe.columns[0]].str.split(';', expand=True)
    # Extract only the number from the 'OTU_Num' column (after 'seqs=')
    dataframe['OTU_Num'] = dataframe['OTU_Num'].str.split('=').str[1].astype(
        int)
    # Drop the original column
    dataframe = dataframe.drop(columns=dataframe.columns[0])
    # Create a reference ID column
    dataframe['Reference_ID'] = dataframe[1].str.split(';').str[0]
    # Creating a Kingdom columm
    dataframe['Kingdom'] = dataframe[1].str.split(':').str[1].str.split(
        ',').str[0]
    # Creating a Domain columm
    dataframe['Domain'] = dataframe[1].str.split(':').str[2].str.split(
        ',').str[0]
    # Creating a Phyllum columm
    dataframe['Phyllum'] = dataframe[1].str.split(':').str[3].str.split(
        ',').str[0]
    # Creating a Class columm
    dataframe['Class'] = dataframe[1].str.split(':').str[4].str.split(
        ',').str[0]
    # Creating an Order columm
    dataframe['Order'] = dataframe[1].str.split(':').str[5].str.split(
        ',').str[0]
    # Creating a Family columm
    dataframe['Family'] = dataframe[1].str.split(':').str[6].str.split(
        ',').str[0]
    # Creating a Genus columm
    dataframe['Genus'] = dataframe[1].str.split(':').str[7].str.split(
        ',').str[0]
    # Creating a Species columm
    dataframe['Species'] = dataframe[1].str.split(':').str[8].str.split(
        ',').str[0]
    # Dropping the original column
    dataframe = dataframe.drop(columns=dataframe.columns[0])
    # Creating a pident (percentage of identical positions) column
    dataframe['Pident'] = dataframe[2]
    # Creating a length (alignment length (sequence overlap)) column
    dataframe['Length'] = dataframe[3]
    # Drop unnecessary columns
    dataframe = dataframe.drop(
        dataframe.columns[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]], axis=1)
    # Rename the deep sea samples so they will appear next to each other on the plot
    dataframe['Sample'] = dataframe['Sample'].replace({
        'A3':
        'deep_sea_A3',
        'X17007':
        'deep_sea_X17007'
    })
    # Rearrange the columns
    cols = [col for col in dataframe if col != 'Sample'] + ['Sample']
    dataframe = dataframe[cols]
    return dataframe


def prepare_data_for_plotting(dataframe,
                              taxonomic_level,
                              threshold_percentage=5,
                              grouping='unique_otus'):
    '''
    Prepare data for taxonomic composition plotting.

    Parameters:
    - filtered_samples (pd.DataFrame): The filtered data containing taxonomic information.
    - taxonomic_level (str): The taxonomic level for the analysis (e.g., "Domain", "Phylum").
    - threshold_percentage (int, optional): The threshold below which taxa are grouped into 'Others'. Default is 5%.
    - grouping (str, optional): How to group the data. Options: 'unique_otus' or 'abundances'. Default is 'unique_otus'.

    Returns:
    - pivot_table (pd.DataFrame): The modified pivot table for plotting.
    - excluded_columns (list): A list of columns to be excluded from the legend.
    '''

    # Filter out the 'Mock' sample and 'Bacteria_X' domain
    filtered_samples = dataframe[(dataframe['Sample'] != 'Mock') &
                                 (dataframe[taxonomic_level] != 'Bacteria_X')]

    # Pivot table to get counts of sequences for each domain per sample
    if grouping == 'abundances':
        pivot_table = filtered_samples.groupby(
            ['Sample', taxonomic_level])['OTU_Num'].sum().unstack(fill_value=0)
    elif grouping == 'unique_otus':
        pivot_table = filtered_samples.groupby(
            ['Sample', taxonomic_level])['OTU'].nunique().unstack(fill_value=0)
    else:
        raise ValueError(
            "Invalid value for 'grouping'. Choose 'unique_otus' or 'abundances'."
        )

    # Calculate proportions within each sample
    proportions = pivot_table.divide(pivot_table.sum(axis=1), axis=0) * 100

    # Determine which entries fall below the threshold
    small_entries = proportions < threshold_percentage

    # Create 'Others' column
    pivot_table['Others'] = 0

    for sample in pivot_table.index:
        # Mask for the current sample's small entries
        mask = small_entries.loc[sample]

        # Sum up small groups into 'Others' for the current sample
        pivot_table.at[sample,
                       'Others'] = pivot_table.loc[sample,
                                                   mask.index[mask]].sum()

        # Zero out the original small group values for the current sample
        pivot_table.loc[sample, mask.index[mask]] = 0

    # Rename the 'Others' column to include the threshold percentage
    others_label = f'Others (Taxa below {threshold_percentage}%)'
    pivot_table.rename(columns={'Others': others_label}, inplace=True)

    # Identify columns to exclude from the legend
    excluded_columns = small_entries.all(
        axis=0)  # Find taxa that are consistently "small" across all samples
    excluded_columns = excluded_columns[excluded_columns].index.tolist()

    return pivot_table, excluded_columns


def plot_taxonomic_data(pivot_table,
                        excluded_columns,
                        taxonomic_level,
                        grouping,
                        sample_name_mapping=False,
                        save_path=False,
                        subgroup=False):
    '''
    Plot the taxonomic composition per sample based on the provided pivot table.

    Parameters:
    - pivot_table (pd.DataFrame): The pivot table containing the data to be plotted.
    - excluded_columns (list): A list of columns to be excluded from the legend.
    - taxonomic_level (str): The taxonomic level for the analysis (e.g., "Domain", "Phylum").
    - grouping (str): How to group the data. Options: 'unique_otus' or 'abundances'.

    Returns:
    - None
    '''
    # Define the color palette
    colors = sns.color_palette(
        'husl',
        len(pivot_table.columns) - len(excluded_columns))

    # Plotting
    ax = pivot_table.drop(columns=excluded_columns).plot(kind='bar',
                                                         stacked=True,
                                                         figsize=(10, 7),
                                                         color=colors)
    plt.title(f'Taxonomic Composition Per Sample at {taxonomic_level} Level')

    # Rename x labels
    if sample_name_mapping:
        new_labels = [
            sample_name_mapping.get(label.get_text(), label.get_text())
            for label in ax.get_xticklabels()
        ]
        ax.set_xticklabels(
            new_labels, rotation=45,
            ha='right')  # Adjusting rotation and alignment for readability
    plt.xlabel('Sample')

    # Disable the top and right borders around the figure
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    # Adjust the legend layout
    if subgroup:
        legend_title = subgroup
    else:
        legend_title = None

    plt.legend(loc='center left',
               bbox_to_anchor=(1.0, 0.5),
               frameon=False,
               title=legend_title)

    if grouping == 'abundances':
        plt.ylabel('OTU abundance')
    else:
        plt.ylabel('Total Number of unique OTUs')

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=600)

    plt.show()
