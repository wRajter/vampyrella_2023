import requests
import time
import xml.etree.ElementTree as ET
import os

MAX_RETRIES = 50
TIME_BETWEEN_RETRIES = 150  # in seconds
RAW_DATA = os.path.join('..', '..', 'raw_data')
OUTPUT_PATH = os.path.join(RAW_DATA, 'blast_results')


def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        sequence = ''
        header = None
        for line in file:
            if line.startswith('>'):
                if header:
                    sequences[header] = sequence
                header = line.strip().replace('>', '')
                sequence = ''
            else:
                sequence += line.strip()
        if header:
            sequences[header] = sequence
    return sequences


def run_blast(sequence):
    base_url = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'

    params = {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        'DATABASE': 'nt',
        'QUERY': sequence
    }
    response = requests.post(base_url, params=params)
    response.raise_for_status()

    rid = None
    for line in response.text.split('\n'):
        if 'RID = ' in line:
            rid = line.split('RID = ')[1].split(' ')[0].strip()
            break

    if not rid:
        raise ValueError('Failed to retrieve RID from BLAST response.')

    print(f'RID: {rid}')

    for attempt in range(MAX_RETRIES):
        print(f'Waiting for results... (Attempt {attempt+1}/{MAX_RETRIES})')
        time.sleep(TIME_BETWEEN_RETRIES)
        params = {'CMD': 'Get', 'FORMAT_TYPE': 'XML', 'RID': rid}
        result = requests.get(base_url, params=params)
        result.raise_for_status()

        if '<BlastOutput_iterations>' in result.text:
            return result.text

    raise ValueError(
        f'Failed to retrieve BLAST results after {MAX_RETRIES} attempts.')


def extract_blast_details_from_xml_to_tsv(xml_content):
    root = ET.fromstring(xml_content)

    # Prepare the header for the TSV file
    rows = ['Score\tE_value\tMax_Ident\tHit_accession\tOrganism']

    # Iterate over the top 20 hits
    for hit in root.findall('.//Hit')[:20]:
        hit_accession = hit.find('Hit_accession').text
        hit_def = hit.find('Hit_def').text
        score = hit.find('.//Hsp_bit-score').text
        e_value = hit.find('.//Hsp_evalue').text
        identities = float(hit.find('.//Hsp_identity').text)
        align_len = float(hit.find('.//Hsp_align-len').text)
        max_ident_percentage = (identities / align_len) * 100

        # Construct the row for this hit
        row = f'{score}\t{e_value}\t{max_ident_percentage:.2f}%\t{hit_accession}\t{hit_def}'
        rows.append(row)

    # Convert the list of rows to a single string with new line separation
    tsv_content = '\n'.join(rows)

    return tsv_content


if __name__ == '__main__':

    # Read FASTA files
    fasta_file = 'low_score_seqs.fasta'
    sequences = read_fasta(fasta_file)

    # BLASTing the FASTA files and saving them as the XML files
    for header, sequence in sequences.items():
        time.sleep(10)
        print(f'Running BLAST for {header}...')
        blast_output = run_blast(sequence)
        output_file_path = os.path.join(OUTPUT_PATH, f'{header}.xml')
        with open(output_file_path, 'w') as file:
            file.write(blast_output)
            print(f'{header} saved to {header}.xml...')

    # Extracting BLAST results from XML to TSV files:
    filenames = os.listdir(OUTPUT_PATH)
    for filename in filenames:
        if filename.endswith(".xml"):
            with open(os.path.join(OUTPUT_PATH, filename), 'r') as xml_file:
                xml_content = xml_file.read()

            tsv_content = extract_blast_details_from_xml_to_tsv(xml_content)
            file_name_without_ext = os.path.splitext(filename)[0]

            with open(
                    os.path.join(OUTPUT_PATH, f'{file_name_without_ext}.tsv'),
                    'w') as tsv_file:
                tsv_file.write(tsv_content)
