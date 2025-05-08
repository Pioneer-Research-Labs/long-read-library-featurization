import argparse
from os.path import join
from create_library_dataframe import (load_data, 
    calculate_raw_read_lengths,
    merge_and_qc_data,
    aggregate_data
)

if __name__ == "__main__":
    # Get library path and genome identifier with argparse
    parser = argparse.ArgumentParser(description="Pipeline for creating a featurized library object.")
    parser.add_argument("--library_path", type=str, required=True, help="Path to the library files.")
    args = parser.parse_args()

    library_path = args.library_path

    # Load data
    print('Loading data...')
    barcodes, mapped_inserts, = load_data(library_path)
    print('Data loaded.')

    # Calculate raw read lengths
    print('Calculating raw read lengths...')
    raw_read_lengths = calculate_raw_read_lengths(library_path)
    print('Raw read lengths calculated.')
    
    # Merge and QC data
    print('Merging and QCing data...')
    merge = merge_and_qc_data(barcodes, mapped_inserts, raw_read_lengths)
    print('Data merged and QCed.')

    # Aggregate data
    print('Aggregating data...')
    merge = aggregate_data(merge)
    print('Data aggregated.')

    # Export dataframe object
    out_path = library_path
    filename = 'library_dataframe.parquet'
    
    merge.to_parquet(join(out_path, filename), index=False)
    print('Object exported.')