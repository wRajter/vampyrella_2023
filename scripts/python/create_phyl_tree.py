# Python functions for the create_phyl_tree Jypyter notebook:
# ../../notebooks/create_phyl_tree.ipynb

# Import modules
from Bio import AlignIO
import subprocess
import os
import re
import shutil

# Constants
GBLOCKS_PATH = os.path.join('..', 'raw_data', 'packages', 'Gblocks_0.91b', 'Gblocks')


# Function definitions

def run_mafft(input_file, output_file, log_file, method='localpair', maxiterate_num=1000):
    '''
    Run MAFFT alignment on the input file using the L-INS-i or G-INS-i strategy (recommended for <200 sequences).
    L-INS-i (method=localpair): probably most accurate; iterative refinement method incorporating local pairwise alignment information
    G-INS-i (method=globalpair): suitable for sequences of similar lengths; iterative refinement method incorporating global pairwise alignment information
    maxiterate_num: number cycles of iterative refinement are performed. Default: 1000
    More details in MAFFT manual: https://mafft.cbrc.jp/alignment/software/manual/manual.html
    '''

    cmd = ['mafft', f'--{method}', '--maxiterate', f'{maxiterate_num}', input_file]

    with open(output_file, 'w') as out_f, open(log_file, 'w') as log_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=log_f, text=True, check=True)

    if result.returncode == 0:
        print(f'Alignment was saved to: {output_file}')
        print(f'Log file was saved to: {log_file}')

    return result.returncode


def alignment_stats(alignment_file):
    '''
    Calculate and return various statistics for a given sequence alignment.

    This function computes the number of sequences, alignment length, average number of gaps per sequence,
    percentage of gaps, and percentage identity across all sequences in the alignment.

    Returns:
    num_sequences (int): The total number of sequences in the alignment.
    alignment_length (int): The length of the alignment (number of columns).
    average_gaps (float): The average number of gap characters ('-') per sequence.
    percentage_gaps (float): The percentage of the alignment that consists of gaps.
    percentage_identity (float): The average percentage identity across all pairwise sequence comparisons.

    Raises:
    ValueError: If the alignment length is zero, indicating an empty or invalid alignment file.

    Note:
    The percentage identity calculation here is a simplified approach suitable for small datasets. It may not be
    accurate for large or complex alignments, as it only considers direct matches in each column without considering
    phylogenetic relationships or alignment quality.
    '''

    # Load the alignment
    alignment = AlignIO.read(alignment_file, "fasta")

    # Number of sequences
    num_sequences = len(alignment)

    # Alignment length
    alignment_length = alignment.get_alignment_length()

    if alignment_length == 0:
        raise ValueError(f"The alignment in {alignment_file} has a length of zero. Please check the file.")

    # Calculate the total number of gaps across all sequences
    total_gaps = sum(str(record.seq).count('-') for record in alignment)

    # Average number of gaps per sequence
    average_gaps = total_gaps / num_sequences

    # Percentage of gaps across the sequences
    percentage_gaps = (average_gaps / alignment_length) * 100

    # Percentage identity (for simplicity, we'll just calculate this for pairwise comparisons)
    total_columns = alignment_length * num_sequences
    identical_columns = sum(col.count(col[0]) for col in zip(*alignment))
    percentage_identity = (identical_columns / total_columns) * 100

    return num_sequences, alignment_length, average_gaps, percentage_gaps, percentage_identity



def run_raxmlng_check(alignment, output_dir, model='GTR+G', prefix='T1'):
    '''
    Run raxml-ng to check the alignment and save the output in the specified directory.
    '''
    # Specify command for raxm-ng
    cmd = [
        'raxml-ng', '--check',
        '--msa', alignment,
        '--model', model,
        '--prefix', prefix
    ]
    # Run the command via the subrocess module
    results = subprocess.run(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                check=True)

    # Print output and error messages (if any)
    print(results.stdout)
    if results.stderr:
        print(results.stderr)


    # Move files to the output directory
    for file in os.listdir():
        if file.startswith(prefix):
            os.rename(file, os.path.join(output_dir, file))
    if results.returncode == 0:
        print(f'Output files were moved to: {output_dir}')

    return results.returncode


def run_raxmlng_tree(alignment, output_dir, num_threads=3, model='GTR+G', prefix='T2'):
    '''
    Run raxml-ng to compute the maximum likelihood tree and save the output in the specified directory.
    '''

    # Specify command for raxm-ng
    cmd = [
        'raxml-ng',
        '--msa', alignment,
        '--model', model,
        '--prefix', prefix,
        '--threads', str(num_threads),
        '--seed', '2',
        '--tree', 'pars{25},rand{25}'
    ]

    # Run the command via the subrocess module
    results = subprocess.run(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             text=True,
                             check=True)

    # Print output and error messages (if any)
    print(results.stdout)
    if results.stderr:
        print(results.stderr)

    # Move files to the output directory
    for file in os.listdir():
        if file.startswith(prefix):
            os.rename(file, os.path.join(output_dir, file))
    if results.returncode == 0:
        print(f'Output files were moved to: {output_dir}')

    return results.returncode


def run_raxmlng_bootstrap(alignment, output_dir, num_threads=3, model='GTR+G', prefix='T3'):
    '''
    """
    Runs RAxML-NG to perform a bootstrap analysis on a given sequence alignment.

    RAxML-NG is called with bootstrap mode enabled, using the specified number of threads,
    evolutionary model, and other parameters. The output files generated by RAxML-NG are
    then moved to the specified output directory.

    Parameters:
    alignment (str): Path to the alignment file in FASTA format.
    output_dir (str): Directory where RAxML-NG output files will be stored.
    num_threads (int, optional): Number of threads RAxML-NG will use for the computation. Defaults to 3.
    model (str, optional): The substitution model used for RAxML-NG. Defaults to 'GTR+G'.
    prefix (str, optional): The prefix for output files from RAxML-NG. Defaults to 'T3'.

    Returns:
    int: The return code from the RAxML-NG process. Zero indicates success, while non-zero indicates an error.

    Raises:
    CalledProcessError: An error occurred in RAxML-NG.
    '''

    # Specify command for raxm-ng
    cmd = [
        'raxml-ng',
        '--bootstrap',
        '--msa', alignment,
        '--model', model,
        '--prefix', prefix,
        '--threads', str(num_threads),
        '--seed', '2'
    ]

    # Run the command via the subrocess module
    results = subprocess.run(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             text=True,
                             check=True)

    # Print output and error messages (if any)
    print(results.stdout)
    if results.stderr:
        print(results.stderr)

    # Move files to the output directory
    for file in os.listdir():
        if file.startswith(prefix):
            os.rename(file, os.path.join(output_dir, file))
    if results.returncode == 0:
        print(f'Output files were moved to: {output_dir}')

    return results.returncode



def run_raxmlng_support(tree, output_dir, bootstraps, num_threads=3, prefix='T4'):
    '''
    Runs RAxML-NG to map bootstrap support values onto the best-scoring ML tree.

    This function takes a given maximum likelihood (ML) tree and a set of bootstrap trees,
    then uses RAxML-NG to calculate support values for each branch of the ML tree. The output,
    including the ML tree with mapped support values, is then saved to the specified output directory.

    Parameters:
    tree (str): Path to the best-scoring ML tree file.
    output_dir (str): Directory where RAxML-NG output files will be stored.
    bootstraps (str): Path to the bootstrap trees file.
    num_threads (int, optional): Number of threads RAxML-NG will use for the computation. Defaults to 3.
    prefix (str, optional): The prefix for output files from RAxML-NG. Defaults to 'T4'.

    Returns:
    int: The return code from the RAxML-NG process. Zero indicates success, while non-zero indicates an error.

    Raises:
    CalledProcessError: An error occurred in RAxML-NG.
    '''

    # Specify command for raxm-ng
    cmd = [
        'raxml-ng',
        '--support',
        '--tree', tree,
        '--bs-trees', bootstraps,
        '--prefix', prefix,
        '--threads', str(num_threads),
    ]

    # Run the command via the subrocess module
    results = subprocess.run(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             text=True,
                             check=True)

    # Print output and error messages (if any)
    print(results.stdout)
    if results.stderr:
        print(results.stderr)

    # Move files to the output directory
    for file in os.listdir():
        if file.startswith(prefix):
            os.rename(file, os.path.join(output_dir, file))
    if results.returncode == 0:
        print(f'Output files were moved to: {output_dir}')

    return results.returncode


def filter_fasta_file(file_path, exclude_pattern='XXXXXXXX'):
    '''
    Reads a FASTA file and filters out sequences whose headers contain a specified pattern.

    This function processes a FASTA file line by line. If a header line (starting with '>')
    contains the specified pattern, the function skips this line and its associated sequence.
    Lines not matching the exclusion pattern are retained.

    Args:
    file_path (str): The path to the FASTA file to be processed.
    exclude_pattern (str, optional): A string pattern to search for in the header lines.
                                     Defaults to 'XXXXXXXX'. If a header contains this pattern,
                                     the function excludes the entire sequence (header and sequence data).

    Returns:
    list: A list of strings, where each string is a line from the original FASTA file that
          does not include the exclusion pattern in its header.

    '''
    with open(file_path, 'r') as infile:
        lines = infile.readlines()

    filtered_lines = []
    skip = False

    for line in lines:
        if line.startswith('>'):
            skip = exclude_pattern in line
        if not skip:
            filtered_lines.append(line)
    return filtered_lines


def shorten_sequence_names(alignment):
    shorten_seq_names = []
    for line in alignment:
        if line.startswith('>'):
            tax_name = line.split('|')[-1].strip()
            gi = line.split('|')[0].strip('>').strip()
            new_name = '>' + tax_name + '_' + gi + '\n'
            shorten_seq_names.append(new_name)
        else:
            shorten_seq_names.append(line)

    return shorten_seq_names


def num_seqs(alignment_path):
    '''
    Returns number sequences from FASTA file that is provide as an input (alignment_path)
    '''
    with open(alignment_path, 'r') as infile:
        count_ = 0
        for line in infile:
            if line.startswith('>'):
                count_ += 1
    return count_


def run_gblocks(input_alignment, b1, b2, b3, b4, b5, name_specifier='', gblocks_path = GBLOCKS_PATH):

    # Variables
    directory_path = os.path.dirname(input_alignment)
    input_file_name, _ = os.path.splitext(os.path.basename(input_alignment))
    # Specify the original and desired output names
    original_fasta_output = input_alignment + ".gb"
    fasta_output = os.path.join(directory_path, f'{input_file_name}_gblocks{name_specifier}.fasta')
    original_html_output = input_alignment + ".gb.htm"
    html_output = os.path.join(directory_path, f'{input_file_name}_gblocks{name_specifier}.fasta.html')
    # Path for the output .log file
    log_file_path = os.path.join(directory_path, f'{input_file_name}_gblocks{name_specifier}.log')

    # Gblocks parameters
    cmd = [gblocks_path,
           input_alignment,
           '-t=d',
           f'-b1={b1}',
           f'-b2={b2}',
           f'-b3={b3}',
           f'-b4={b4}',
           f'-b5={b5}',
           '-e=.gb']

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Rename the files
    os.rename(original_fasta_output, fasta_output)
    os.rename(original_html_output, html_output)

    # Grab info from the gblocks HTML file and save it to the log file
    # Open the HTML file and read its content
    with open(html_output, 'r') as file:
        content = file.read()

        # Using regex to find content between <pre> and </pre>
        match = re.search(r'<pre>.*?<b>Parameters used</b>(.*?)</pre>', content, re.DOTALL)

        # Start with the Gblocks stdout content
        log_content = f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}\n\n\nParameters used:\n"

        if match:
            extracted_info = match.group(1).strip()  # Extract the content
            log_content += extracted_info  # Append the extracted info to the log content
        else:
            log_content += "Pattern not found in the HTML file!"

        # Write the combined content to the .log file
        with open(log_file_path, 'w') as log_file:
            log_file.write(log_content)


def run_gblocks_grid_search(mafft_alignment,
                            b1_values,
                            b2_values,
                            b3_values,
                            b4_values,
                            b5_values,
                            gblocks_path):
    '''
    Execute Gblocks using a grid of parameter settings to fine-tune the alignment process.

    Parameters:
    - mafft_alignment (str): Path to the alignment file generated using MAFFT.
    - b1_values (list): List of values for b1 parameter.
    - b2_values (list): List of values for b2 parameter.
    - b3_values (list): List of values for b3 parameter.
    - b4_values (list): List of values for b4 parameter.
    - b5_values (list): List of values for b5 parameter.
    - gblocks_path (str): Path to the Gblocks executable.

    Returns:
    None. However, the function will produce alignment files with various parameters settings based on the grid search.
    '''
    for b5_value in b5_values:
        for i in range(len(b1_values)):
            b1 = b1_values[i]
            b2 = b2_values[i]
            b3 = b3_values[i]
            b4 = b4_values[i]
            b5 = b5_value
            run_gblocks(input_alignment = mafft_alignment,
                        b1 = b1,
                        b2 = b2,
                        b3 = b3,
                        b4 = b4,
                        b5 = b5,
                        name_specifier=f'_grid_{b5_value}_{i}',
                        gblocks_path = gblocks_path)

def move_search_grid_files(directory_path, grid_search_path):
    '''
    Transfer Gblocks output files to a designated 'grid_search' directory and delete all associated .html files.

    Parameters:
    - directory_path (str): Path to the directory containing the output files from Gblocks.
    - grid_search_path (str): Destination path for the 'grid_search' directory where output files will be moved.

    Returns:
    None. However, the function will move the specified files and delete any .html files in the target directory.
    '''
    # Move all the output files to the separate 'grid_search' directory
    files_to_move = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if 'grid' in f.lower() and os.path.isfile(os.path.join(directory_path, f))]

    grid_search_path = os.path.join(directory_path, 'grid_search')

    # Check if directory already exists, if not, create it
    if not os.path.exists(grid_search_path):
        os.mkdir(grid_search_path)

    for file in files_to_move:
        file_name = os.path.basename(file)  # Extract just the filename from the full path
        destination_path = os.path.join(grid_search_path, file_name)  # Construct the destination path
        shutil.move(file, destination_path)

    # delete all the html files
    for file in os.listdir(grid_search_path):
        if file.endswith('.html'):
            file_path = os.path.join(grid_search_path, file)
            os.remove(file_path)
