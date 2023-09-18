import pandas as pd
import re
from Bio import Entrez
import time


def extract_data_from_filename(blast_table, pattern):
    '''
    Extract relevant details from the filename using the provided regex pattern.
    '''
    match = pattern.match(blast_table)
    if not match:
        print(f'Filename {blast_table} did not match the pattern!')
        return None

    return match.groups()


def read_and_augment_csv(filename, sequence_id, taxon_name, score_tax_assign,
                         sample_name):
    '''
    Read the CSV into a dataframe and augment with additional columns.
    '''
    df = pd.read_csv(filepath_or_buffer=filename, sep='\t')
    df['Sample'] = sample_name
    df['Sequence_id'] = sequence_id
    df['Tax_assignment'] = taxon_name
    df['Score_tax_assign'] = score_tax_assign
    df = df.rename(columns={'Organism': 'BLASTn'})
    # Rearrange columns
    cols = [
        'Sequence_id', 'Sample', 'Tax_assignment', 'Score_tax_assign',
        'BLASTn', 'Max_Ident', 'Score', 'E_value', 'Hit_accession'
    ]
    # Extract only the first two words (representing species name) from the BLASTn column
    df['BLASTn'] = df['BLASTn'].apply(lambda x: ' '.join(x.split()[:2]))

    df = df[cols]
    return df


def get_taxonomy_from_genbank(genbank_id):
    '''
    Fetch the taxonomy lineage from GenBank using a given accession ID.

    Parameters:
    - genbank_id (str): The GenBank accession ID.

    Returns:
    - str: A string representation of the taxonomy lineage if successful.
    - None: If there's an error fetching the taxonomy.

    Notes:
    - This function makes use of the Entrez API and requires an email
      to be set using Entrez.email = 'your_email@example.com' before using.
    '''
    try:
        handle = Entrez.efetch(db='nucleotide',
                               id=genbank_id,
                               rettype='gb',
                               retmode='xml')
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.5)
        return record[0]['GBSeq_taxonomy']
    except Exception as e:
        print(f'Failed to fetch taxonomy for {genbank_id}. Error: {e}')
        return None  # or some default value or behavior you'd like


def process_hit_accession(x):
    '''
    Fetch the taxonomy lineage from GenBank for a given accession ID.

    Parameters:
    - x (str): The GenBank accession ID.

    Returns:
    - str: The taxonomy lineage corresponding to the accession ID.
           Returns None if fetching fails.

    Prints:
    - A message indicating the GenBank ID currently being processed.
    '''
    print(f'Processing GenBank ID: {x}')
    return get_taxonomy_from_genbank(x)


def add_taxopath_column(df):
    '''
    Add a 'Taxopath' column to a dataframe based on the 'Hit_accession' column.

    Parameters:
    - df (pandas.DataFrame): A dataframe with a 'Hit_accession' column
                             containing GenBank accession IDs.

    Returns:
    - pandas.DataFrame: The input dataframe augmented with a 'Taxopath' column
                        containing taxonomy lineage for each accession ID.
    '''
    df['Taxopath'] = df['Hit_accession'].apply(process_hit_accession)
    return df
