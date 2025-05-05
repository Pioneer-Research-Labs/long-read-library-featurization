### Create a featurized AnnData object from outputs of long read qc pipeline

import pysam
import pandas as pd
import numpy as np
from os.path import join
import anndata as ad
from tqdm import tqdm
from scipy.sparse import csr_matrix

MAPQ_CUTOFF = 60 # cutoff for mapping quality score
GENOME_PATH = 's3://pioneer-data/genomes/' # path to genomes

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

    with pysam.FastxFile(fastq_path) as fq:
        # Using tqdm with unknown total
        with tqdm(desc="Processing reads", unit="read") as pbar:
            for entry in fq:
                read_ids.append(entry.name)
                lengths.append(len(entry.sequence))
                pbar.update(1)

    df = pd.DataFrame({
        'read_id': read_ids,
        'length': lengths
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

def load_bedtools_intersection(file):
    '''
    Helper function to load bedtools intersect.out file into Pandas dataframe
    '''
    return pd.read_csv(file, 
                        sep='\t', 
                        header=None, 
                        names=['insert_chr',
                               'insert_start',
                               'insert_end',
                               'ID',
                               'mapping_quality',
                               'insert_sense',
                               'feature_chr',
                               'feature_start',
                               'feature_end',
                               'feature_name',
                               'intersect',
                               'feature_sense',
                               'overlap_bp'])

def load_gff_file(file):
    '''
    Helper function to load genome annotated .gff file into Pandas dataframe
    '''
    return pd.read_csv(file, 
                  sep='\t', 
                  header=None,
                  names=['genome','scheme','feature_type','feature_start','feature_end','score','strand','frame','info'])

def sparse_pivot_table(frame, row_label, col_label, val_label):
    '''
    Pivot table in pandas using sparse dtype for fill value and returns the sparse matrix
    '''
    from scipy.sparse import csr_matrix
    from pandas.api.types import CategoricalDtype
    
    row_cats = CategoricalDtype(sorted(frame[row_label].unique()), ordered=True)
    col_cats = CategoricalDtype(sorted(frame[col_label].unique()), ordered=True)
    
    row = frame[row_label].astype(row_cats).cat.codes
    col = frame[col_label].astype(col_cats).cat.codes
    
    sparse_matrix = csr_matrix((frame[val_label], (row, col)), \
                               shape=(row_cats.categories.size, col_cats.categories.size))



    return pd.DataFrame.sparse.from_spmatrix(sparse_matrix, index=row_cats.categories, columns=col_cats.categories)

def load_data(library_path, genome):
    '''
    Load data from long read library qc output and annotated genome
    '''
    # Load barcodes
    barcodes = pd.read_csv(join(library_path,
                            'barcodes.tsv'), sep='\t', index_col=0, names=['ID','sequence','','length'])
    # Mapped inserts
    mapped_inserts = load_bedtools_intersection(join(library_path,
                                                    'insert_intersect.out'))
    
    # Load gff file
    gff = load_gff_file(join(GENOME_PATH, genome, genome + '_genes.gff'))
    gff['gene_id'] = [x.split('ID=')[1].split('gene-')[1].split(';')[0] for x in gff['info']]
    gff['feature_name'] = [x.split('Name=')[1].split(';')[0] for x in gff['info']]

    # Add gene id column to mapped inserts dataframe
    mapped_inserts = pd.merge(mapped_inserts, 
                        gff[['gene_id','feature_start','feature_end','feature_name']], 
                        on=['feature_start','feature_end','feature_name'])
    
    # Get just the unique inserts in a separate dataframe without the features
    inserts_only = mapped_inserts[['ID','insert_chr','insert_start','insert_end','insert_sense','mapping_quality']].drop_duplicates().set_index('ID')
    
    return barcodes, mapped_inserts, inserts_only

def merge_and_qc_data(barcodes, inserts_only, raw_read_lengths):
    '''
    Merge barcodes and inserts dataframes and perform QC
    '''

    # Merge to get data that contains both barcodes and inserts
    merge = pd.merge(
        barcodes, inserts_only, left_index=True, right_index=True, how='inner').dropna()
    
    # Filter out reads with poor mapping quality
    merge = merge[merge.mapping_quality >= MAPQ_CUTOFF]

    # Calculate mapped insert length
    merge['insert_length'] = merge['insert_end'] - merge['insert_start']

    # Convert barcode length to int
    merge['length'] = merge['length'].astype(int)

    # Add raw read length
    merge['raw_read_length'] = [raw_read_lengths.loc[x].length for x in merge.index]

    # Also get the reads for empty barcodes with no inserts
    empty_barcodes = list(set(barcodes.sequence.unique()) - set(merge.sequence.unique()))
    empty_reads = barcodes[barcodes['sequence'].isin(empty_barcodes)]

    # Add raw read length to empty reads
    empty_reads['raw_read_length'] = [raw_read_lengths.loc[x].length for x in empty_reads.index]

    return merge, empty_reads

def create_features_and_metadata(merge, mapped_inserts, empty_reads):
    '''
    Create features and metadata dataframes for barcodes and features
    '''
    # Create the barcodes metadata dataframe and average by barcode
    bc_meta = merge[
        ['sequence','length','insert_chr','insert_start','insert_end','insert_length','insert_sense','raw_read_length']]
    bc_meta = bc_meta.rename(columns={'sequence':'bc_sequence',
                                    'length':'bc_length'})
    bc_meta_avg = bc_meta.groupby('bc_sequence').agg(lambda x: '|'.join(map(str, x)))

    # Create features metadata dataframe
    meta_cols = ['feature_start',
                'feature_end',
                'feature_name',
                'gene_id']
    features_meta = mapped_inserts[meta_cols].set_index('gene_id').drop_duplicates()

    # Subset to intersections present in barcode map
    mapped_inserts = mapped_inserts[mapped_inserts['ID'].isin(bc_meta.index)]

    # Add a column for fraction of gene contained
    mapped_inserts['overlap_frac'] = mapped_inserts['overlap_bp'] / (
        mapped_inserts['feature_end'] - mapped_inserts['feature_start'])

    # Add a column for one-hot encoding of overlap
    mapped_inserts['overlap_ohc'] = [1 if x == 1 else 0 for x in mapped_inserts.overlap_frac]
    
    # Merge barcode information and average overlap metrics across unique barcodes using median
    mapped_inserts['bc_sequence'] = mapped_inserts['ID'].map(bc_meta['bc_sequence'].to_dict())

    mapped_inserts_avg = mapped_inserts.groupby(['bc_sequence','gene_id'])[['overlap_frac','overlap_ohc']].median().reset_index()

    # Create a BC x feature matrix using one-hot encoded full/empty genes
    features_full = sparse_pivot_table(mapped_inserts_avg, 'bc_sequence', 'gene_id', 'overlap_ohc')

    # Using partial genes
    features_partial = sparse_pivot_table(mapped_inserts_avg, 'bc_sequence', 'gene_id', 'overlap_frac')

    # Save empty barcodes
    empty_reads = empty_reads.rename(columns={'sequence':'bc_sequence',
                                'length':'bc_length'})
    empty_barcodes = empty_reads.groupby('bc_sequence').agg(lambda x: '|'.join(map(str, x)))

    return features_full, features_partial, bc_meta_avg, features_meta, empty_barcodes

def create_anndata_object(features_full, features_partial, bc_meta_avg, features_meta, empty_barcodes):
    '''
    Create anndata object from features and metadata
    '''

    # Create AnnData object
    obj = ad.AnnData(X=csr_matrix(features_full.sparse.to_coo()),
                        obs=bc_meta_avg,
                        var=features_meta.reindex(features_full.columns),
                    uns={'empty_barcodes':empty_barcodes})
    obj.layers['partial'] = csr_matrix(features_partial.sparse.to_coo())

    return obj

def export_anndata_object(obj, out_path, filename):
    '''
    Export AnnData object to h5ad file
    '''
    breakpoint()
    # Export
    obj.write_h5ad(join(out_path, filename))