import argparse
from os.path import join
from os import listdir
from tqdm import tqdm
import os
import pandas as pd
from create_library_dataframe import (load_data, 
    calculate_raw_read_lengths,
    merge_and_qc_data,
    aggregate_data,
    save_to_fasta
)

if __name__ == "__main__":
    # Get library path and genome identifier with argparse
    parser = argparse.ArgumentParser(description="Pipeline for creating a featurized library object.")
    parser.add_argument("--library_path", type=str, required=True, help="Path to the library files.")
    parser.add_argument("--tesseract_library", action=argparse.BooleanOptionalAction, help="Bool to handle multi-genome tesseract library.")
    args = parser.parse_args()

    library_path = args.library_path
    tesseract_library = args.tesseract_library

    # If tesseract library type, 
    if tesseract_library:
        print('Handling tesseract library type...')
        # Return subdirectories in library_path
        samples = [f for f in listdir(library_path) if not f.startswith('.') and 
                      not f.endswith('.DS_Store') and 
                      os.path.isdir(join(library_path, f))]
        library_paths = [join(library_path, sample) for sample in samples]
        
    else:
        library_paths = [library_path]
        samples = [library_path.split('/')[-1]]  # Use the last part of the path as the sample name
    
    # Process for each sample
    aggregated_merge = []
    for sample, library_path_sample in tqdm(zip(samples, library_paths), desc="Processing samples", total=len(samples)):
        print(f'Processing sample: {sample}')
        
        # Load data
        print('Loading data...')
        barcodes, mapped_inserts, = load_data(library_path_sample)
        print('Data loaded.')

        # Calculate raw read lengths
        print('Calculating raw read lengths...')
        raw_read_lengths = calculate_raw_read_lengths(library_path_sample)
        print('Raw read lengths calculated.')
        
        # Merge and QC data
        print('Merging and QCing data...')
        merge = merge_and_qc_data(barcodes, mapped_inserts, raw_read_lengths)
        print('Data merged and QCed.')

        # Aggregate data
        print('Aggregating data...')
        merge = aggregate_data(merge).assign(library_sample=sample)
        print('Data aggregated.')

        # Append to aggregated merge
        aggregated_merge.append(merge)

    # Concatenate all aggregated dataframes
    aggregated_merge = pd.concat(aggregated_merge)
    print('All samples processed and aggregated.')

    # Remove duplicate barcodes
    n_pre_drop = len(aggregated_merge)
    aggregated_merge.drop_duplicates(subset='bc_sequence', keep=False, inplace=True)
    n_post_drop = len(aggregated_merge)
    print(f'Dropped {n_pre_drop - n_post_drop} degenerate barcodes.')
    print(f'Library barcode number after dropping duplicates: {n_post_drop}')

    # Export dataframe object
    out_path = library_path
    filename = 'library_dataframe.parquet'
    
    aggregated_merge.to_parquet(join(out_path, filename), index=False)
    print('Object exported.')

    # Save barcodes to a separate FASTA file for error correction
    barcodes_list = list(aggregated_merge['bc_sequence'])
    save_to_fasta(barcodes_list, join(out_path, 'library_barcodes.fasta'))
    print('Barcodes saved to FASTA file.')
