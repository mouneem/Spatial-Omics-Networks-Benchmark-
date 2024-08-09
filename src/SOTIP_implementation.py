from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import issparse
from PIL import Image
import json
import random
import scanpy as sc
import sys
path_src = '../../src'
sys.path.append(path_src)

from sotip import *
import os
from scipy.spatial import KDTree
from scipy.sparse import csr_matrix

def get_anndata(coordinates, matrix, features, observations, neighbors, dim_red , k = 6):
    observations = pd.read_csv(observations, index_col=0)
    features = pd.read_csv(features, index_col=0)
    matrix = pd.read_csv(matrix, index_col=0)

    coordinates = (
        pd.read_csv(coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )
    features = features['phenotype'].unique()
    adata = ad.AnnData(obs=observations, var=features, obsm={"spatial": coordinates})

    if matrix.shape[0] > 0:
        X = matrix
        
        if issparse(X):
            X = X.tocsr()
        adata.X = X

    if neighbors != None:
        adata.obsp["spatial_connectivities"] = pd.read_csv(neighbors, index_col=0).values
    else:
        # Generate adjacency matrix if not provided
        tree = KDTree(coordinates)
        distances, indices = tree.query(coordinates, k=k + 1)  # k+1 because the closest point is itself
        indptr = []
        indices_list = []
        for i, neighbors in enumerate(indices):
            indptr.extend([i] * k)
            indices_list.extend(neighbors[1:])  # Skip the point itself

        data = [1] * len(indices_list)  # All connections are 1 (binary adjacency matrix)
        adjacency_matrix = csr_matrix((data, (indptr, indices_list)), shape=(coordinates.shape[0], coordinates.shape[0]))
        adata.obsp["spatial_connectivities"] = adjacency_matrix

    return adata

def run_sotip(
    coordinates,
    matrix,
    features,
    observations,
    neighbors,
    out_dir,
    dim_red,
    n_clusters,
    technology,
    seed,
    config
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)


    # Output files
    label_file = out_dir / "domains.tsv"
    embedding_file = out_dir / "embedding.tsv"


    adata = get_anndata(coordinates, matrix, features, observations, neighbors, dim_red)

    random.seed(seed)
    
    sc.pp.pca(adata)
    # Process the data with scanpy routine    
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Find resolution for the given n_clusters
    res_recom = search_res(adata, n_clusters, start=0.1, step=0.1, tol=5e-3, max_run=10)
    sc.tl.leiden(adata, resolution=res_recom)

    # ME size 
    knn = config.get('knn')
    n_neighbors = config.get('n_neighbours')

    # Order of cluster label for ME representation (adata.obsm['ME'])
    ME_var_names_np_unique = np.array(adata.obs['leiden'].cat.categories) 

    # Add a ME obsm for adata
    MED(adata, use_cls='leiden', nn=knn, copy=False, ME_var_names_np_unique=ME_var_names_np_unique, spatial_var='spatial') 


    return (adata)

# def compute_paga(adata, config):
#     # Compute the topological structure with paga

#     connectivities = adata.obsp['connectivities']
#     print(len(connectivities.data))

#     sc.tl.paga(adata, groups='leiden')
#     sc.pl.paga(adata, add_pos=True, show=False)
#     knn = config.get('knn')
#     n_neighbors = config.get('n_neighbors')
#     # Use the connectivities between cell clusters to guide the graph distance computation
#     gd = get_ground_distance(adata, method='paga_guided_umap', cls_key='leiden')

#     # Add a X_ME_EMD_mat obsm to adata_phEMD
#     adata_phEMD = MED_phEMD_mp(
#         adata.copy(),         
#         GD_method='paga_guided_umap',  
#         MED_knn=knn,          
#         CT_obs='leiden',       
#         ifspatialplot=False,  
#         OT_method='pyemd',    
#         ME_precompyted=True,  
#         GD_precomputed=True,  
#     )

#     adata.obsp['ME_EMD_mat'] = adata_phEMD.obsm['X_ME_EMD_mat']

#     # Compute the MEG, each node is a ME, edge is the connectivity
#     sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="reduced_dimensions")
#     knn_indices, knn_dists, forest = sc.neighbors.compute_neighbors_umap(adata_phEMD.obsm['X_ME_EMD_mat'], n_neighbors=n_neighbors, metric='precomputed')
#     adata.obsp['distances'], adata.obsp['connectivities'] = sc.neighbors._compute_connectivities_umap(
#         knn_indices,
#         knn_dists,
#         adata.shape[0],
#         n_neighbors, 
#     )

#     return (adata)

# def compute_MEgraph(adata, config):
#     out_dir = config.get('out_dir')
#     label_file = out_dir / "domains.tsv"
#     embedding_file = out_dir / "embedding.tsv"
    
#     # Set the ME graph's associated information (connectivity matrix, distance matrix) to neighbors_EMD
#     adata.uns['neighbors_EMD'] = adata.uns['neighbors'].copy()

#     # Use the computed MEG as input of umap and leiden clustering
#     sc.tl.umap(adata, neighbors_key='neighbors_EMD')
#     sc.tl.leiden(adata, neighbors_key='neighbors_EMD', key_added='leiden_EMD')

#     # Merge regions according to MEG connectivities
#     sc.tl.paga(adata, groups='leiden_EMD', neighbors_key='neighbors_EMD')
#     merge_cls_paga(adata, thresh=0, min_cls=n_clusters, paga_plot=False)

#     embedding_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names) # optional, DataFrame with index (cell-id/barcode) and n columns
#     label_df = adata.obs[['leiden_EMD']]  # DataFrame with index (cell-id/barcode) and 1 column (label)

#     ## Write output
#     label_df.columns = ["label"]
#     label_df.to_csv(label_file, sep="\t", index_label="")

#     if embedding_df is not None:
#         embedding_df.to_csv(embedding_file, sep="\t", index_label="")
#     return adata