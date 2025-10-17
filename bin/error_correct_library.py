### Script to error correct library barcodes to base barcode Illumina sequencing clustered barcodes

import pandas as pd
from os.path import join
import os
from Levenshtein import distance
import subprocess
from utils import read_fasta_to_dataframe
from create_library_dataframe import save_to_fasta
import argparse
import shutil

GROUND_TRUTH_PATH = '/seqlrg/base_barcode_sequencing'
VSEARCH_ID_CUTOFF = '0.95'
MAX_DISTANCE = 2
N_THREADS = '64'

def load_library(library_path):
    '''
    Load library dataframe parquet file
    '''
    library = pd.read_parquet(join(library_path, 'library_dataframe.parquet'))

    # Check that there are no degenerate barcodes
    assert len(library) == len(library.drop_duplicates(subset='bc_sequence'))

    return pd.read_parquet(join(library_path, 'library_dataframe.parquet'))

def load_ground_truth_barcodes(construct_name):
    '''
    Load ground truth barcodes for a given construct
    '''
    ground_truth_barcodes = read_fasta_to_dataframe(join(GROUND_TRUTH_PATH, construct_name + '_barcodes.fasta'))
    return ground_truth_barcodes

def quantify_error_correction_metrics(library, base_barcodes):
    '''
    Quantify error correction metrics
    '''
    base_set = set(base_barcodes.sequence) # Set of base barcode sequences
    library_set = set(library.bc_sequence) # Set of library barcode sequences

    # Number of overlapping barcodes
    exact_match = library_set.intersection(base_set)
    print(f'Number of exact matches: {len(exact_match)}')

    # Number of barcodes that are in the library but not in the base set
    library_only = library_set - base_set
    print(f'Number of barcodes in library only: {len(library_only)}')

    # What fraction of barcodes are exact matches?
    fraction_exact_match = len(exact_match) / len(library_set)
    print(f'Fraction of barcodes that are exact matches: {fraction_exact_match}')

    return library_set, base_set, exact_match

def export_barcodes_for_error_correction(library_set, base_set, library_path, construct_name):
    '''
    Export barcodes for error correction
    '''
    # Get barcodes to correct
    barcodes_to_correct = list(library_set - base_set)

    # Save to library path
    save_to_fasta(barcodes_to_correct, library_path, f'{construct_name}_barcodes_to_correct.fasta')

    # Also export to a temporary local directory
    folder_name = library_path.rstrip('/').split('/')[-1]
    temp_dir = f'/tmp/vsearch_input_{folder_name}'
    os.makedirs(temp_dir, exist_ok=True)
    save_to_fasta(barcodes_to_correct, temp_dir, f'{construct_name}_barcodes_to_correct.fasta')

    return temp_dir

def run_vsearch_on_barcodes_to_correct(temp_dir, construct_name):
    '''
    Run vsearch usearch to correct barcodes with Python subprocess
    '''
    # Run vsearch usearch to correct barcodes
    temp_barcodes_to_correct_path = join(temp_dir, f'{construct_name}_barcodes_to_correct.fasta')
    vsearch_output_path = join(temp_dir, f'{construct_name}_vsearch_correction_output.tsv')

    cmds = [
        'vsearch',
        '--usearch_global',
        temp_barcodes_to_correct_path,
        '--db',
        f'{GROUND_TRUTH_PATH}/{construct_name}_barcodes.fasta',
        '--id',
        VSEARCH_ID_CUTOFF,
        '--uc',
        vsearch_output_path,
        '--threads',
        N_THREADS,
    ]

    # Run the vsearch command and print stdout in real time
    process = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"vsearch failed with exit code {process.returncode}")

    # Copy vsearch results to library path
    # If the library path is on s3, copy to s3. Otherwise, copy to local.
    output_path = join(library_path, f'{construct_name}_vsearch_correction_output.tsv')
    if library_path.startswith('s3://'):
        subprocess.run(['aws', 's3', 'cp', vsearch_output_path, output_path])
    else:
        shutil.copy(vsearch_output_path, output_path)

    return output_path

def load_vsearch_results(vsearch_output_path):
    '''
    Load vsearch results
    '''
    names = [
    'record_type',
    'target_seq_ordinal_number',
    'seq_length',
    'similarity',
    'match_orientation',
    'na1',
    'na2',
    'alignment',
    'query_label',
    'target_label'
    ]
    vsearch_results = pd.read_csv(vsearch_output_path, sep='\t', names=names)

    # How many hits do we have i.e. correctable barcodes?
    print('Number of correctable barcodes:', len(vsearch_results[vsearch_results['record_type'] == 'H']))

    # How many barcodes are not hit?
    print('Number of barcodes not hit:', len(vsearch_results[vsearch_results['record_type'] == 'N']))

    # How many barcodes are ambiguous?
    print('Number of ambiguous barcodes:', len(vsearch_results[vsearch_results['record_type'] == 'U']))

    return vsearch_results

def correct_barcodes(temp_dir, construct_name, base_barcodes, vsearch_results, library, exact_match):
    '''
    Correct barcodes
    '''
    # Load library barcodes
    barcodes_to_correct_df = read_fasta_to_dataframe(join(temp_dir, f'{construct_name}_barcodes_to_correct.fasta'))
    barcodes_to_correct_df = barcodes_to_correct_df.rename(columns={'id':'query_label',
                                                                    'sequence':'query_seq'})
    base_barcodes = base_barcodes.rename(columns={'id':'target_label',
                                                'sequence':'target_seq'})
    correction_map = vsearch_results[vsearch_results.record_type == 'H'][['query_label','target_label']]
    correction_map = pd.merge(correction_map, barcodes_to_correct_df, on='query_label')
    correction_map = pd.merge(correction_map, base_barcodes, on='target_label')

    # Calculate Levenshtein distances to verify correction
    correction_map['levenshtein_distance'] = [distance(x,y) for x,y in zip(correction_map.query_seq,
                                                                        correction_map.target_seq)]
    # Filter for Levenshtein distance <= 2
    correction_map = correction_map[correction_map.levenshtein_distance <= MAX_DISTANCE]

    # Make correction dictionary
    correction_dict = correction_map.set_index('query_seq')['target_seq'].to_dict() # Corrected barcodes

    # Create corrected library
    library['library_correction_status'] = 'uncorrected'
    library.loc[library.bc_sequence.isin(list(exact_match)),'library_correction_status'] = 'exact_match'
    library.loc[library.bc_sequence.isin(correction_map.query_seq),'library_correction_status'] = 'corrected'

    library['corrected_bc_sequence'] = None
    library.loc[library.library_correction_status == 'exact_match','corrected_bc_sequence'] = \
        library.loc[library.library_correction_status == 'exact_match','bc_sequence']
    library.loc[library.library_correction_status == 'corrected','corrected_bc_sequence'] = \
        library.loc[library.library_correction_status == 'corrected','bc_sequence'].map(correction_dict)

    # Drop unmatched barcodes
    library = library[library.library_correction_status != 'uncorrected']

    # Drop raw barcode column and rename corrected barcode column to bc_sequence
    library = library.drop(columns=['bc_sequence']).rename(columns={'corrected_bc_sequence':'bc_sequence'})

    # Save construct used for correction
    library['parent_construct_name'] = construct_name

    return library


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_path', type=str, required=True)
    parser.add_argument('--construct_name', type=str, required=True)
    args = parser.parse_args()

    library_path = args.library_path
    construct_name = args.construct_name

    # Load library
    print(f'Loading library from {library_path}/library_dataframe.parquet')
    library = load_library(library_path)
    print(f'Library loaded.')

    # Load ground truth barcodes
    print(f'Loading ground truth barcodes from {GROUND_TRUTH_PATH}/{construct_name}_barcodes.fasta')
    base_barcodes = load_ground_truth_barcodes(construct_name)
    print(f'Ground truth barcodes loaded.')

    # Quantify error correction metrics
    library_set, base_set, exact_match = quantify_error_correction_metrics(library, base_barcodes)

    # Export barcodes for error correction
    print(f'Exporting barcodes for error correction to {library_path}/{construct_name}_barcodes_to_correct.fasta')
    temp_dir = export_barcodes_for_error_correction(library_set, base_set, library_path, construct_name)
    print(f'Barcodes exported.')

    # Run vsearch on barcodes to correct
    print(f'Running vsearch on barcodes to correct')
    vsearch_output_path = run_vsearch_on_barcodes_to_correct(temp_dir, construct_name)
    print(f'Vsearch completed.')

    # Load vsearch results
    vsearch_results = load_vsearch_results(vsearch_output_path)
    print(f'Vsearch results loaded.')

    # Correct library barcodes
    print(f'Correcting library barcodes')
    library = correct_barcodes(temp_dir, construct_name, base_barcodes, vsearch_results, library, exact_match)
    print(f'Library barcodes corrected.')

    # Drop duplicate barcodes (except for empty inserts)
    # First sort by match type and the number of reads and keep the first occurrence
    library = library.sort_values(by=['library_correction_status', 'n_reads'], ascending=[False, False])
    idx_to_drop = (library.duplicated(subset='bc_sequence', keep='first')) & (library.insert_type != 'empty')
    library = library[~idx_to_drop]
    n_dropped_barcodes = idx_to_drop.sum()
    n_mapped_inserts = len(library[library.insert_type == 'mapped'])
    n_empty_inserts = len(library[library.insert_type == 'empty'])
    n_unmapped_inserts = len(library[library.insert_type == 'unmapped'])
    print(f'Dropped {n_dropped_barcodes} duplicate barcodes.')
    print(f'Library barcode number after dropping duplicates: {len(library)}')
    print(f'Number of mapped inserts: {n_mapped_inserts}')
    print(f'Number of empty inserts: {n_empty_inserts}')
    print(f'Number of unmapped inserts: {n_unmapped_inserts}')

    # Save corrected library
    print(f'Saving corrected library to {library_path}/corrected_library_dataframe.parquet')
    library.to_parquet(join(library_path, 'corrected_library_dataframe.parquet'))
    print(f'Corrected library saved.')

    # Delete temporary directory
    shutil.rmtree(temp_dir)