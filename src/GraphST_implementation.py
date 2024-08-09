import json
from pathlib import Path
import random
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn import metrics
import anndata as ad
import scipy as sp
from GraphST import GraphST
from GraphST.utils import clustering

def map_technology(input_technology):
    technology_mapping = {
        "Visium": "10X",
        "Stereo-seq": "Stereo",
        "Slide-seq": "Slide"
    }
    return technology_mapping.get(input_technology, None)

from scipy.sparse import issparse
from scipy.spatial import KDTree
from scipy.sparse import csr_matrix


def get_anndata(coordinates, matrix, features, observations, neighbors = None, k = 6):
    # Load data
    observations_df = pd.read_csv(observations, index_col=0)  # cells
    features_df = pd.read_csv(features, index_col=0)  # genes
    matrix_df = pd.read_csv(matrix, index_col=0)

    # Ensure the matrix has the correct shape (cells x genes)
    matrix_data = matrix_df.values

    # Align the observations and matrix
    observations_df = observations_df.loc[matrix_df.index]


    # Align the features and matrix
    features_df = matrix_df

    coordinates_array = pd.read_csv(coordinates, index_col=0).loc[observations_df.index].to_numpy()

    # Check the shape of coordinates
    if coordinates_array.shape[0] != observations_df.shape[0]:
        raise ValueError(f"Shape mismatch: coordinates have {coordinates_array.shape[0]} rows, but observations have {observations_df.shape[0]} rows.")


    # Create AnnData object
    adata = ad.AnnData(X=matrix_data, obs=observations_df, var=features_df.T, obsm={"spatial": coordinates_array})

    adata.var['highly_variable'] = True

    if matrix_df.shape[0] > 0:
        X = matrix_df
        
        if issparse(X):
            X = X.tocsr()
        adata.X = X
        
    if neighbors != None:
        adata.obsp["spatial_connectivities"] = pd.read_csv(neighbors, index_col=0).values
    else:
        # Generate adjacency matrix if not provided
        tree = KDTree(coordinates_array)
        distances, indices = tree.query(coordinates_array, k=k + 1)  # k+1 because the closest point is itself
        indptr = []
        indices_list = []
        for i, neighbors in enumerate(indices):
            indptr.extend([i] * k)
            indices_list.extend(neighbors[1:])  # Skip the point itself

        data = [1] * len(indices_list)  # All connections are 1 (binary adjacency matrix)
        adjacency_matrix = csr_matrix((data, (indptr, indices_list)), shape=(coordinates_array.shape[0], coordinates_array.shape[0]))
        adata.obsp["spatial_connectivities"] = adjacency_matrix

    return adata

def run_graphST(coordinates, matrix, features, observations, neighbors, config , out_dir):
    label_file = Path(out_dir) / "domains.tsv"
    embedding_file = Path(out_dir) / "embedding.tsv"

    n_clusters = config['n_clusters']
    seed = config['seed']
    technology = map_technology(config['technology'])
    
    if technology is None:
        raise Exception(
            f"Invalid technology. GraphST only supports 10X Visium, Stereo-seq, and Slide-seq/Slide-seqV2 not {config['technology']}."
        )

    adata = get_anndata(
        coordinates = coordinates, 
        matrix = matrix, 
        features = features, 
        observations = observations, 
        neighbors = neighbors
    )
    adata.var_names_make_unique()
    sc.pp.scale(adata, zero_center=False, max_value=10)

    random.seed(seed)
    torch.manual_seed(seed)
    np.random.seed(seed)

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    model = GraphST.GraphST(adata, datatype=technology, device=device)
    adata = model.train()

    radius = config.get("radius", 50)
    clustering(adata, 
               n_clusters=n_clusters, 
               radius=radius, 
               method=config["method"], 
               refinement=config["refine"]
               )

    label_df = adata.obs[["domain"]]
    out_dir.mkdir(parents=True, exist_ok=True)
    label_df.columns = ["label"]
    label_df.to_csv(label_file, sep="\t", index_label="")
    return label_df

