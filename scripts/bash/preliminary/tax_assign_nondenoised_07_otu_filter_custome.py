# custome based filtering of rare OTUs per sample

## load modules
import os

## variables
project = 'Jamy_2022'
cell = 'cell'
marker = 'rDNA'
sim = 'sim97'
raw_data = '../raw_data'
otu_preclust_dir = f'{raw_data}/OTU_preclust/{project}/{marker}/{cell}/{sim}'
tax_assign_dir = f'{raw_data}/tax_assign_results/{project}/{marker}/{cell}/{sim}'
filt_otu_dir = f'{raw_data}/OTU_filtered/{project}/{marker}/{cell}/{sim}'
vamp_spec_dir = f'{raw_data}/vamp_specific_seqs/{project}/{marker}/{cell}/{sim}'
otu_chim_filt = f'{raw_data}/OTU_nonchimeric/{project}/{marker}/{cell}/{sim}'


## functions
def strip(file_path, output_path):
    '''
    Remove line breaks in FASTA files, so it will be easier to filter singeltons in the next step
    '''
    seqs = []
    # read a file and stripped it
    with open(file_path, 'rt') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                seqs.append('\n' + line)
            else:
                line_cleaned = line.rstrip()
                seqs.append(line_cleaned)
        seqs[0] = seqs[0].lstrip()
    # write the file into fasta
    # get sample name for our data:
    sample_name = file_path.split('/')[-1].lstrip('otu_').rstrip('.fasta')
    with open(f'{output_path}/otu_{sample_name}.fasta', 'w') as fp:
        for seq in seqs:
            fp.write(seq)


def filter_rare(path, threshold=1):
    seqs_nonsingl = []
    with open(path, 'rt') as f:
        lines = f.readlines()
        for line in lines[::-1]:
            if not line.startswith('>'):
                seq_ = line
            if ('>' in line) and (not line.endswith(f'seqs={threshold}\n')):
                seqs_nonsingl.append(line)
                seqs_nonsingl.append(seq_)
    # write the file into fasta
    # get sample name for our data:
    sample_name = path.split('/')[-1].lstrip('otu_').rstrip('.fasta')
    with open(f'{filt_otu_dir}/nonrare_otu_{sample_name}.fasta', 'w') as fp:
        for seq in seqs_nonsingl:
            fp.write(seq)


## code

# Get paths
all_files = os.listdir(otu_chim_filt)
all_paths = [otu_chim_filt + '/' + path for path in all_files]
all_paths

# filter singeltons out
for sample_path in all_paths:
    strip(file_path=sample_path, output_path=otu_chim_filt)
    sample_name = sample_path.split('/')[-1].lstrip('otu_').rstrip('.fasta')
    filter_rare(path=sample_path, threshold=1)

## Extract vamp sequences
# paths to fasta files
all_files = os.listdir(filt_otu_dir)
all_files = [sample for sample in all_files if "nonrare_otu_" in sample]
all_paths = [filt_otu_dir + '/' + path for path in all_files]

# getting vampyrellid-specific ids
for path in all_paths:
    sample_name = path.split('/')[-1].lstrip('nonrare_otu_').rstrip('.fasta')
    tax_assign_path = f'{tax_assign_dir}/blast6_{sample_name}.tab'
    with open(tax_assign_path, 'rt') as f:
        # get vamp ids from the tax assignment
        vamp_ids = []
        for line in f:
            if 'Vampyrellida' in line:
                vamp_ids.append(line.split('\t')[0])
        # use the vamp ids to extract the vamp-specific sequences from the fasta file
        filtered_seqs = []
        with open(path, 'rt') as f:
            lines = f.readlines()
            for line in lines[::-1]:
                if not line.startswith('>'):
                    sequence = line
                else:
                    if any([x in line for x in vamp_ids]):
                        filtered_seqs.append(line)
                        filtered_seqs.append(sequence)
        # write the filtered vampyrellid sequences into a fasta file
        with open(f'{vamp_spec_dir}/nonrare_otu_{sample_name}.fasta',
                  'w') as fp:
            for seq in filtered_seqs:
                fp.write(seq)
