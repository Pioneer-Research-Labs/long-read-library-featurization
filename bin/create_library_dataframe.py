### Create a featurized dataframe object from outputs of long read qc pipeline

### PSEUDOCODE
# - load barcodes and mapped inserts
# - calculate raw read lengths from fastplong output
# - qc barcodes and inserts based on mapping quality and barcode length (is_unmapped and mapping_quality)
# - merge barcodes and inserts on the read ID, use left join on barcodes
# - append raw read length
# - aggregate to barcode ID
# - QC: throw out inserts that have degenerate mapping, then retain median per barcode
# - save as .pq file

### Data scheme
# columns - barcode ID, mapped insert start/end (as pd.interval), insert sequence, insert length, empty boolean, median raw read length
# Throw away is_read1, is_read2, is_proper_pair

import pysam
import pandas as pd
from os.path import join, exists
from tqdm.auto import tqdm
import subprocess
import os
import warnings

MIN_BARCODE_LENGTH = 42
MAX_BARCODE_LENGTH = 46
MAX_EMPTY_INSERT_LENGTH = 100
MAPQ_CUTOFF = 60

def fastq_to_df(fastq_path):
    """
    Reads a FASTQ file using pysam and returns a DataFrame
    with read IDs and sequence lengths. Displays a progress bar.

    Parameters:
    - fastq_path (str): Path to the FASTQ file (can be gzipped).

    Returns:
    - pandas.DataFrame: A DataFrame with columns 'read_id' and 'length'.
    """
    read_ids = []
    lengths = []

    # Check if fastplong output exists
    if not exists(fastq_path):
        print('Fastplong output not found. Skipping raw read calculation.')
        return None
    else:
        print('Fastplong output found. Calculating raw read lengths.')

        with pysam.FastxFile(fastq_path) as fq:
            # Using tqdm with unknown total
            with tqdm(desc="Processing reads", unit="read") as pbar:
                for entry in fq:
                    read_ids.append(entry.name)
                    lengths.append(len(entry.sequence))
                    pbar.update(1)

        df = pd.DataFrame({
            'read_id': read_ids,
            'raw_read_length': lengths
        })

        return df.set_index('read_id')

def calculate_raw_read_lengths(library_path):
    """
    Calculate the raw read lengths from a FASTQ file.

    Parameters:
    - fastq_path (str): Path to the FASTQ file (can be gzipped).

    Returns:
    - pandas.DataFrame: A DataFrame with columns 'read_id' and 'length'.
    """
    return fastq_to_df(join(library_path, 'fastplong.fq'))

def load_bam_to_dataframe(bam_path, max_reads=None):
    """
    Load reads from a BAM file into a pandas DataFrame.

    Parameters:
        bam_path (str): Path to the BAM file.
        max_reads (int, optional): Max number of reads to load. Useful for large files.

    Returns:
        pd.DataFrame: DataFrame with read properties.
    """
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    records = []

    for i, read in enumerate(bamfile.fetch(until_eof=True)):
        if max_reads is not None and i >= max_reads:
            break

        records.append({
            'read_id': read.query_name,
            'reference_name': bamfile.get_reference_name(read.reference_id) if read.reference_id != -1 else None,
            'reference_start': read.reference_start,
            'reference_end': read.reference_end,
            'mapping_quality': read.mapping_quality,
            'is_unmapped': read.is_unmapped,
            'is_reverse': read.is_reverse,
            'insert_sequence': read.query_sequence
        })

    bamfile.close()
    return pd.DataFrame(records)

def load_data(library_path):
    '''
    Load data from long read library qc output (barcodes, mapped inserts, and empty inserts)
    '''
    # Load barcodes
    barcodes = pd.read_csv(join(library_path,
                            'barcodes.tsv'), sep='\t', index_col=0, names=['read_id','bc_sequence','blank','bc_length'])
    barcodes = barcodes.drop(columns=['blank'])

    # Mapped inserts
    # If the bam file is on s3, download to a temporary local file and then load
    if library_path.startswith('s3://'):
        # Extract the last folder in library_path to define as folder_name
        folder_name = library_path.rstrip('/').split('/')[-1]
        print('Input files are on s3. Downloading mapped inserts.bam to temporary local file.')
        temp_name = join('/tmp/', folder_name, 'mapped_inserts.bam')
        subprocess.run(['aws', 's3', 'cp', join(library_path, 'mapped_inserts.bam'), temp_name])
        mapped_inserts = load_bam_to_dataframe(temp_name).set_index('read_id')
        os.remove(temp_name)
        print('Mapped inserts.bam downloaded to temporary local file and loaded.')
    else:
        mapped_inserts = load_bam_to_dataframe(join(library_path, 'mapped_inserts.bam')).set_index('read_id')
    
    # Filter out barcodes that have weird lengths
    barcodes = barcodes[barcodes.bc_length.between(MIN_BARCODE_LENGTH, MAX_BARCODE_LENGTH)]
    
    # Drop unmapped inserts or inserts with poor mapping quality
    mapped_inserts = mapped_inserts[mapped_inserts['is_unmapped'] == False].drop(columns=['is_unmapped'])
    mapped_inserts = mapped_inserts[mapped_inserts['mapping_quality'] >= MAPQ_CUTOFF]
    mapped_inserts['insert_type'] = 'mapped'

    # Unmapped inserts -- load, drop blank column, label as "empty" or "other"
    unmapped_inserts = pd.read_csv(join(library_path,
        'sites.tsv'), sep='\t', index_col=0, names=['read_id','unmapped_insert_sequence','blank','unmapped_insert_length'])
    unmapped_inserts = unmapped_inserts.drop(columns=['blank'])
    unmapped_inserts.loc[unmapped_inserts.unmapped_insert_length > MAX_EMPTY_INSERT_LENGTH, 'insert_type'] = 'unmapped'
    unmapped_inserts.loc[unmapped_inserts.unmapped_insert_length <= MAX_EMPTY_INSERT_LENGTH, 'insert_type'] = 'empty'
    
    return barcodes, mapped_inserts, unmapped_inserts

def merge_and_qc_data(barcodes, mapped_inserts, unmapped_inserts, raw_read_lengths):
    '''
    Merge barcodes and inserts dataframes and perform QC
    '''

    # Merge reads to get data that contains both barcodes and inserts
    merge = pd.merge(
        barcodes, mapped_inserts, left_index=True, right_index=True, how='inner')
    
    # Add raw read lengths
    if raw_read_lengths is not None:
        merge = pd.merge(
            merge, raw_read_lengths, left_index=True, right_index=True, how='left')
    else:
        merge['raw_read_length'] = None

    # Calculate mapped insert length
    merge['insert_length'] = [len(x) if type(x) is str else None for x in merge['insert_sequence']]

    # Merge barcodes with unmapped inserts
    unmapped_inserts = pd.merge(
        barcodes, unmapped_inserts, left_index=True, right_index=True, how='inner')
    
    # Remove any unmapped inserts that have a shared barcode with a mapped insert
    unmapped_inserts = unmapped_inserts[~unmapped_inserts.index.isin(merge.index)]

    # Add unmapped inserts to merged dataframe
    merge = pd.concat([merge, unmapped_inserts])

    return merge

def aggregate_data(merge):
    '''
    Aggregate data to barcode level. For barcodes with multiple reads, use the modal insert. If degenerate,
    use the longest insert.
    '''
    tqdm.pandas()
    # Calculate reads per each barcode
    reads_per_barcode = merge.bc_sequence.value_counts()

    # Construct initial aggregate dataframe with barcodes that have single reads (or empty barcodes)
    single_read_barcodes = reads_per_barcode[reads_per_barcode == 1].index
    df_agg = merge[merge.bc_sequence.isin(single_read_barcodes)].copy()
    df_agg = df_agg.reset_index(drop=True) # Drop read ID column
    df_agg['n_reads'] = 1

    # For barcodes with multiple reads, check for a modal insert based on length
    multi_read_barcodes = reads_per_barcode[reads_per_barcode > 1].index
    merge_multi_reads = merge[merge.bc_sequence.isin(multi_read_barcodes)]
    
    def resolve_group(group):
        lengths = group['insert_length']
        mode_lengths = lengths.mode()

        if len(mode_lengths) == 1:
            mode_len = mode_lengths.iloc[0]
            modal_subset = group[group['insert_length'] == mode_len]
            # Aggregate modal insert reads
            # Ignore RunTime warnings about means of empty slices due to nan values
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=RuntimeWarning)
                return modal_subset.groupby('bc_sequence').agg({
                    'bc_length': 'first',
                    'reference_name': 'first',
                    'reference_start': 'median',
                    'reference_end': 'median',
                    'mapping_quality': 'median',
                    'is_reverse': 'first',
                    'insert_sequence': 'first',
                    'raw_read_length': 'median',
                    'insert_type': 'first',
                    'insert_length': 'median',
                    'unmapped_insert_sequence': 'first',
                    'unmapped_insert_length': 'median'
                }).assign(n_reads = len(modal_subset))
        else:
            # Degenerate: use the longest insert
            max_len = lengths.max()
            return group[group['insert_length'] == max_len].head(1).set_index('bc_sequence').assign(n_reads = 1)
        
    # Group by barcode and apply resolution logic
    resolved_multi = merge_multi_reads.groupby('bc_sequence', group_keys=False).progress_apply(resolve_group)

    # Reset index and combine
    resolved_multi = resolved_multi.reset_index()
    df_agg = pd.concat([df_agg, resolved_multi])

    return df_agg

def save_to_fasta(seqs, out_path, filename):
    """
    Save a list of DNA/RNA sequences to a FASTA file.

    Parameters:
        seqs (list of str): List of sequence strings.
        filename (str): Output FASTA file path.
    """
    # If the write path is local, write directly. Otherwise if the write path is on s3, upload to s3
    # by first writing to a temporary local file and then using aws s3 cp to upload to s3.
    if out_path.startswith('s3://'):
        out_name = join('/tmp/', filename)
    else:
        out_name = join(out_path, filename)
    # Write to file
    with open(out_name, 'w') as f:
        for i, seq in enumerate(seqs, 1):
            f.write(f">seq{i}\n{seq}\n")
    # If on s3, upload to s3 and delete temporary local file
    if out_path.startswith('s3://'):
        s3_name = join(out_path, filename)
        subprocess.run(['aws', 's3', 'cp', out_name, s3_name])
        os.remove(out_name)
