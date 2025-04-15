import argparse
from create_featurized_object import (load_data, 
    merge_and_qc_data, 
    create_features_and_metadata,
    create_anndata_object,
    export_anndata_object
)

if __name__ == "__main__":
    # Get library path and genome identifier with argparse
    parser = argparse.ArgumentParser(description="Pipeline for creating a featurized library object.")
    parser.add_argument("--library_path", type=str, required=True, help="Path to the library.")
    parser.add_argument("--genome", type=str, required=True, help="Genome identifier.")
    args = parser.parse_args()

    library_path = args.library_path
    genome = args.genome

    # Load data
    barcodes, mapped_inserts, inserts_only = load_data(library_path, genome)

    print('Data loaded.')

    # Merge and QC data
    merge = merge_and_qc_data(barcodes, inserts_only)

    # Create features and metadata
    features_full, features_partial, bc_meta_avg, features_meta, empty_barcodes = create_features_and_metadata(
        merge,
        barcodes,
        mapped_inserts)

    # Create AnnData object
    obj = create_anndata_object(features_full, features_partial, bc_meta_avg, features_meta, empty_barcodes)

    # Export AnnData object
    out_path = library_path
    filename = 'featurized_object.h5ad'

    export_anndata_object(obj, out_path, filename)

    print('Object exported.')