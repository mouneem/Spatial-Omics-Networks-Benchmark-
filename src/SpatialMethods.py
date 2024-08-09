import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from PIL import Image
import json
import random
# import scanit
import scanpy as sc
import torch
import os
import sys

path_src = '.'
sys.path.append(path_src)
import importlib
import SpatialAnalysis as SpAn


path_src = '../mosna/'
sys.path.append(path_src)
import mosna2 as mosna

cluster_params = {
    'reducer_type': 'umap', 
    'n_neighbors': 30, 
    'metric': 'euclidean', #'manhattan', # or 'euclidean',
    'min_dist': 0.0,
    'clusterer_type': 'leiden', 
    'dim_clust': 2, 
    'k_cluster': 30, 
    'resolution_parameter': 0.005,
}

import tempfile

def run_NichesDetection(coords, types, expression, thr = 5, mosna_output_path = '../data/intermediate/MOSNA/', output_path = '../' , data_info = None, plot_network = False,  cluster_params = cluster_params):
    """
    out_path : path to export output files (edges)
    data_info : patient or sample or other...
    """
    edges = SpAn.delaunay_edges(coords)
    edges = SpAn.trim_edges_by_distance(coords, edges, thr )

    if plot_network:
        SpAn.plot_network(coords, edges, types = types, linewidth = .5)
    marker_col = expression.columns

    nodes = pd.concat([coords, expression],axis=1)


    # If no grouping info, add column (sample)
    if data_info == None :
        data_info =['sample']
        nodes[data_info] = '1'

    # exporting results
    os.makedirs(mosna_output_path, exist_ok=True)
    temp_dir = tempfile.TemporaryDirectory(dir=mosna_output_path)
    temp_dir_path = temp_dir.name
    csv_file_path = temp_dir_path + f"/nodes_patient-{data_info[0]}.csv"
    nodes.to_csv(csv_file_path, index=False)

    csv_file_path = temp_dir_path + f"/edges_patient-{data_info[0]}.csv"
    pd.DataFrame(edges).to_csv(csv_file_path, index=False)

    attributes_col = nodes.columns

    print(os.getcwd() )
    var_aggreg = mosna.compute_spatial_omic_features_single_network(
    method = 'NAS',
    net_dir = temp_dir_path,  
    data_info = data_info,
    extension = 'csv',
    read_fct = pd.read_csv,
    attributes_col=attributes_col, 
    use_attributes=marker_col, # use all attributes 
    make_onehot=False, 
    stat_funcs = 'default', 
    stat_names = 'default', 
    id_level_1='patient',
    save_intermediate_results=False, 
    dir_save_interm=None,
    add_sample_info = False,
    verbose=1,
    )

    # plt.gcf().set_size_inches(30, 30)

    cluster_labels, cluster_dir, nb_clust, _ = mosna.get_clusterer(var_aggreg, output_path, **cluster_params)

    return cluster_labels
    