### Integrate long read library data with short read fitness data

import pandas as pd
import numpy as np
from tqdm import tqdm
from os.path import join
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Levenshtein import distance
from io import StringIO
import pysam
import pybedtools as pbt
from scipy.sparse import lil_matrix, vstack, hstack, lil_array, csr_matrix, csr_array
from scipy.spatial.distance import hamming
import anndata as ad
import scanpy as sc
import ipyparallel as ipp
import itertools as it
from time import time
import sys

MIN_BC_LENGTH = 42
MAX_BC_LENGTH = 46
N_CORES = 32
DISTANCE_CUTOFF = 2 # This corresponds to ~0.0001 probability of distances >=3 for a 44bp barcode

def compute_distances_row(lr):    
    '''
    Compute Levenshtein distances between a long read and all short reads in the short read barcode list
    '''
    distances = csr_array((1,len(sr_bcs)))
    
    for i in range(len(sr_bcs)):
        d = distance(lr, sr_bcs[i], score_cutoff=distance_cutoff)
        if d <= distance_cutoff:
            distances[0,i] = d
        else:
            continue
    return distances

def load_data(long_read_object, fitness_data):
    '''
    Load AnnData object for long read library data as well as short read fitness data
    '''
    # Library data
    obj = ad.read_h5ad(long_read_object)

    # Create a column for barcode length
    obj.obs['single_bc_length'] = [len(x) for x in obj.obs.index]

    # Load short read data
    fitness = pd.read_csv(fitness_data, index_col=0)

    return obj, fitness

def get_filtered_barcodes_to_correct(obj, fitness):
    '''
    Get barcodes from long read library data and short read fitness data for error correction.
    We won't bother correcting the empty barcodes because there are a lot of them.
    '''
    # Get the barcodes
    lr_bcs = list(obj.obs.index) # Long read
    sr_bcs = list(fitness.barcode.unique()) # Short read
    df_sr_bcs = pd.DataFrame(index=sr_bcs, 
                            columns=['length'],
                            data=[len(x) for x in sr_bcs]) # dataframe for short read
    
    # For now let's only look at barcodes between a certain length (i.e. 2 Levenshtein distances away from 44)
    lr_filter = list(obj.obs[obj.obs.single_bc_length.between(MIN_BC_LENGTH, MAX_BC_LENGTH)].index)
    sr_filter = list(df_sr_bcs[df_sr_bcs.length.between(MIN_BC_LENGTH, MAX_BC_LENGTH)].index)

    # We only need to correct barcodes that are not exact matches between the two datasets
    intersection = set(lr_filter).intersection(set(sr_filter))

    lr_filter = list(set(lr_filter) - intersection)
    sr_filter = list(set(sr_filter) - intersection)
    
    return lr_filter, sr_filter, intersection

def calculate_distance_matrix(lr_filter, sr_filter):
    '''
    Calculate Levenshtein distance matrix between correctable long read and short read barcodes
    '''
    # Generate inputs
    inputs = lr_filter

    # Parallel computation
    n_cores = N_CORES
    cluster = ipp.Cluster(n=n_cores)
    cluster.start_cluster_sync()
    rc = cluster.connect_client_sync()
    rc.wait_for_engines(n_cores); rc.ids
    dview = rc[:] # Using direct view so we can push global variables

    # push global functions and definitions
    dview.push(dict(distance=distance,
                lil_array=lil_array,
                    csr_array=csr_array,
                sr_bcs=sr_filter,
                distance_cutoff=DISTANCE_CUTOFF))

    # submit the tasks
    asyncresult = dview.map_async(compute_distances_row, inputs)
    # wait interactively for results
    asyncresult.wait_interactive()

    # Shutdown the cluster
    try:
        await cluster.stop_cluster()
    except RuntimeError:
        pass
        
    # retrieve actual results
    results = asyncresult.get()

    # Stack sparse arrays together to get Levenshtein distance matrix
    distances = vstack(results) 

    return distances

def error_correct_barcodes(lr_filter, sr_filter, distances):
    '''
    Error correct long read barcodes based on Levenshtein distance matrix
    '''
    # Get indices of nonzero distances less than cutoff
    nz_indices = distances.nonzero()
    np_lr_sub = np.array(lr_filter) # Numpy array of long read barcodes
    np_sr_sub = np.array(sr_filter) # Numpy array of short read barcodes

    # Make dataframe of correctable barcodes
    correctable_bcs = pd.DataFrame(index=np.unique(np_lr_sub[nz_indices[0]]),
                                columns=['distance','possible_mapped_bc'],
                                data=[])

    # Pull out distances, mapped barcodes, and downstream stats
    for i, j in zip(nz_indices[0], nz_indices[1]):
        try:
            correctable_bcs.loc[np_lr_sub[i],'distance'].append(distances[i,j])
            correctable_bcs.loc[np_lr_sub[i],'possible_mapped_bc'].append(np_sr_sub[j])
        except AttributeError:
            correctable_bcs.loc[np_lr_sub[i],'distance'] = [distances[i,j]]
            correctable_bcs.loc[np_lr_sub[i],'possible_mapped_bc'] = [np_sr_sub[j]]
    correctable_bcs['n'] = [len(x) for x in correctable_bcs['distance']]
    correctable_bcs['min_dist'] = [min(x) for x in correctable_bcs['distance']]
    correctable_bcs['n_min_dist'] = [
        sum(np.array(x)==y) for x,y in zip(correctable_bcs['distance'],correctable_bcs['min_dist'])]

    # Filter out barcodes that have more than one mappable target based on minimum distance
    correctable_bcs = correctable_bcs[correctable_bcs.n_min_dist == 1]
    correctable_bcs['mapped_bc'] = [np.array(x)[np.array(y) == z][0] for x,y,z in zip(correctable_bcs['possible_mapped_bc'],
                                                                        correctable_bcs['distance'],
                                                                        correctable_bcs['min_dist'])]

    # Save the mapping dictionary
    correction_map = correctable_bcs['mapped_bc'].to_dict()

    # Verify that the mapping are the minimum distance
    assert np.array([distance(x,y)==z for x,y,z in zip(
        correctable_bcs.index, 
        correctable_bcs.mapped_bc,
        correctable_bcs.min_dist)]).all()

    return correction_map, correctable_bcs

def create_integrated_object(obj, fitness, intersection, correction_map, correctable_bcs):
    '''
    Create an integrated object with error corrected barcodes
    '''
    # Create a new corrected object
    obj_corr = obj[(obj.obs.index.isin(list(intersection))) | (obj.obs.index.isin(correctable_bcs.index))]

    # Create a new obs column for correction status
    obj_corr.obs.reset_index(inplace=True)
    obj_corr.obs.loc[:,'corrected'] = False
    obj_corr.obs['corrected_bc_sequence'] = obj_corr.obs['bc_sequence'].copy()
    obj_corr.obs.loc[obj_corr.obs.bc_sequence.isin(correctable_bcs.index),'corrected'] = True
    obj_corr.obs.loc[obj_corr.obs.bc_sequence.isin(correctable_bcs.index),'bc_sequence'] = obj_corr.obs.loc[
        obj_corr.obs.bc_sequence.isin(correctable_bcs.index),'corrected_bc_sequence'].map(correction_map)
    obj_corr.obs.set_index('corrected_bc_sequence', inplace=True)

    # Create integrated object by getting intersection of barcodes and using fitness data
    fitness_pivot = pd.pivot_table(fitness[fitness['Day #'] != 0],
                values='fitness',
                index=['Library Info','Media','Replicate #','barcode'],
                columns='Day #').reset_index()
    fitness_pivot.columns.name = None
    fitness_pivot = fitness_pivot.rename(columns={1:'fitness_day_1',
                                                2:'fitness_day_2',
                                                3:'fitness_day_3',
                                                4:'fitness_day_4',
                                                5:'fitness_day_5',
                                                'barcode':'corrected_bc_sequence'})

    # New metadata table
    int_meta = pd.merge(obj_corr.obs,
                        fitness_pivot,
                        on='corrected_bc_sequence',
                        how='left').dropna()
    bc_to_ind = {v: k for k, v in obj_corr.obs.reset_index()['corrected_bc_sequence'].to_dict().items()}

    # Generate new features matrix
    inds = int_meta.corrected_bc_sequence.map(bc_to_ind)
    int_features = obj_corr.X[inds,:]

    # Layers
    int_layers = {}
    for l in obj_corr.layers:
        int_layers[l] = obj_corr.layers[l][inds,:]

    # Create integrated object
    obj_int = ad.AnnData(X=int_features,
                        obs=int_meta,
                        var=obj_corr.var,
                        uns=obj_corr.uns)

    for l in obj_corr.layers:
        obj_int.layers[l] = int_layers[l]