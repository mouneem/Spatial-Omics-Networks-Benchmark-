
import anndata as ad
import numpy as np
import pandas as pd
import scipy as sp
import scanpy as sc
import torch
import SEDR
import os
import warnings
import random
import json
from pathlib import Path
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from SEDR.graph_func import preprocess_graph


warnings.filterwarnings('ignore')
rpy2.robjects.numpy2ri.activate()


def get_anndata(observations_path, features_path, coordinates_path):
    observations = pd.read_csv(observations_path, index_col=0)
    features = pd.read_csv(features_path, index_col=0)

    coordinates = (
        pd.read_csv(coordinates_path, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(obs=observations, var=features, obsm={"spatial": coordinates})

    return adata

def res_search_fixed_clus_leiden(adata, n_clusters, increment=0.01, random_seed=2023):
    for res in np.arange(0.2, 2, increment):
        sc.tl.leiden(adata, random_state=random_seed, resolution=res)
        if len(adata.obs['leiden'].unique()) > n_clusters:
            break
    return res-increment

def leiden(adata, n_clusters, use_rep='SEDR', key_added='SEDR', random_seed=2023):
    sc.pp.neighbors(adata, use_rep=use_rep)
    res = res_search_fixed_clus_leiden(adata, n_clusters, increment=0.01, random_seed=random_seed)
    sc.tl.leiden(adata, random_state=random_seed, resolution=res)

    adata.obs[key_added] = adata.obs['leiden']
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    return adata

def res_search_fixed_clus_louvain(adata, n_clusters, increment=0.01, random_seed=2023):
    for res in np.arange(0.2, 2, increment):
        sc.tl.louvain(adata, random_state=random_seed, resolution=res)
        if len(adata.obs['louvain'].unique()) > n_clusters:
            break
    return res-increment

def louvain(adata, n_clusters, use_rep='SEDR', key_added='SEDR', random_seed=2023):
    sc.pp.neighbors(adata, use_rep=use_rep)
    res = res_search_fixed_clus_louvain(adata, n_clusters, increment=0.01, random_seed=random_seed)
    sc.tl.louvain(adata, random_state=random_seed, resolution=res)

    adata.obs[key_added] = adata.obs['louvain']
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    return adata

def mclust_R(adata, n_clusters, use_rep='SEDR', key_added='SEDR', modelNames="EEE", random_seed=2023):
    np.random.seed(random_seed)
    robjects.r.library("mclust")

    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm[use_rep]), n_clusters, modelNames)
    mclust_res = np.array(res[-2])

    adata.obs[key_added] = mclust_res
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    return adata

def run_sedr(config, out_dir_path, n_clusters, technology, seed, observations_path, features_path, coordinates_path, image_path=None, cluster_method="leiden"):
    # Set up output files
    out_dir = Path(out_dir_path)
    label_file = out_dir / "domains.tsv"
    embedding_file = out_dir / "embedding.tsv"

    # Set SEED for SEDR
    random.seed(seed)
    SEDR.fix_seed(seed)

    # Load data
    adata = get_anndata(observations_path, features_path, coordinates_path)
    adata.var_names_make_unique()

    # Generate neighborhood graph if not provided
    sc.pp.neighbors(adata, use_rep='spatial')

    # Copy from source code in order for customization
    adj_m1 = adata.obsp["connectivities"]
    adj_m1 = sp.sparse.coo_matrix(adj_m1)

    # Store original adjacency matrix (without diagonal entries) for later
    adj_m1 = adj_m1 - sp.sparse.dia_matrix((adj_m1.diagonal()[np.newaxis, :], [0]), shape=adj_m1.shape)
    adj_m1.eliminate_zeros()

    # Some preprocessing
    adj_norm_m1 = preprocess_graph(adj_m1)
    adj_m1 = adj_m1 + sp.sparse.eye(adj_m1.shape[0])

    adj_m1 = adj_m1.tocoo()
    shape = adj_m1.shape
    values = adj_m1.data
    indices = np.stack([adj_m1.row, adj_m1.col])
    adj_label_m1 = torch.sparse_coo_tensor(indices, values, shape)

    norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - adj_m1.sum()) * 2)

    graph_dict = {
        "adj_norm": adj_norm_m1,
        "adj_label": adj_label_m1.coalesce(),
        "norm_value": norm_m1
    }

    # Training SEDR
    # device: using cpu or gpu (if available)
    # using_dec: boolean, whether to use the unsupervised deep embedded clustering (DEC) method to improve clustering results 
    sedr_net = SEDR.Sedr(adata.obsm['spatial'], 
                         graph_dict, 
                         mode='clustering', 
                         device=config["device"])

    if config["using_dec"]:
        sedr_net.train_with_dec(N=1)
    else:
        sedr_net.train_without_dec(N=1)
    sedr_feat, _, _, _ = sedr_net.process()
    # latent embedding
    adata.obsm['SEDR'] = sedr_feat

    # Clustering 
    if cluster_method == "mclust":
        adata = mclust_R(adata, n_clusters=n_clusters, use_rep='SEDR', key_added='SEDR', random_seed=seed)
    elif cluster_method == "louvain":
        adata = louvain(adata, n_clusters=n_clusters, use_rep='SEDR', key_added='SEDR', random_seed=seed)
    elif cluster_method == "leiden":
        adata = leiden(adata, n_clusters=n_clusters, use_rep='SEDR', key_added='SEDR', random_seed=seed)
        
    # Output dataframes
    label_df = adata.obs[["SEDR"]]
    embedding_df = pd.DataFrame(adata.obsm['SEDR'])
    embedding_df.index = adata.obs_names

    # Write output
    out_dir.mkdir(parents=True, exist_ok=True)
    label_df.columns = ["label"]
    label_df.to_csv(label_file, sep="\t", index_label="")

    if embedding_df is not None:
        embedding_df.to_csv(embedding_file, sep="\t", index_label="")
    return label_df
