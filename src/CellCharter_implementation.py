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
import anndata as ad
import cellcharter as cc
import squidpy as sq


def make_minimal_adata(nodes):
    X = pd.get_dummies(nodes['phenotype']).astype(np.float64)
    obs = nodes[['x', 'y']]
    obs['roi'] = 'ROI1'
    # Create 'var' dataframe with the one-hot encoded types columns
    var = pd.DataFrame(index=X.columns)

    # Create AnnData object
    minimal_adata = ad.AnnData(X=X.values, obs=obs, var=var)
    minimal_adata.obsm['spatial'] = nodes[['x', 'y']].values
    return minimal_adata

def run_CellCharter(nodes, n_clusters=2):
    adata = make_minimal_adata(nodes)
    sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)

    cc.gr.remove_long_links(adata)

    cc.gr.aggregate_neighbors(adata, n_layers=3)
    gmm = cc.tl.Cluster(
        n_clusters=2, 
        trainer_params=dict(accelerator='gpu', devices=1)
    )
    gmm.fit(adata, use_rep='X_cellcharter')
    adata.obs['spatial_cluster'] = gmm.predict(adata, use_rep='X_cellcharter')
    return adata
    