import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from PIL import Image
import tempfile
import os, sys
import scanpy as sc
import torch
import random
import importlib.util

def run_deepst(
    coordinates, features, observations, out_dir, n_clusters, technology, seed,
    matrix=None, neighbors=None, dim_red=None, image=None, config=None
):
    out_dir = Path(out_dir)

    # Output files
    label_file = out_dir / "domains.tsv"
    embedding_file = out_dir / "embedding.tsv"

    config = config

    def get_anndata(coordinates, features, observations, matrix, neighbors, dim_red, image):
        observations = pd.read_csv(observations, index_col=0)
        features = pd.read_csv(features, index_col=0)

        coordinates = (
            pd.read_csv(coordinates, index_col=0)
            .loc[observations.index, :]
            .to_numpy()
        )

        adata = ad.AnnData(obs=observations, var=features.T, obsm={"spatial": coordinates})

        if matrix is not None:
            X = pd.read_csv(matrix , index_col=0)
            adata.X = X


        from sklearn.neighbors import NearestNeighbors

        if neighbors is not None:
            adata.obsp["spatial_connectivities"] = sp.io.mmread(neighbors).T.tocsr()
        else:
            # Generate spatial connectivities if neighbors file is not provided
            nbrs = NearestNeighbors(n_neighbors=6).fit(coordinates)
            distances, indices = nbrs.kneighbors(coordinates)
            row_indices = np.repeat(np.arange(len(coordinates)), indices.shape[1])
            col_indices = indices.flatten()
            data = distances.flatten()
            spatial_connectivities = sp.coo_matrix((data, (row_indices, col_indices)), shape=(len(coordinates), len(coordinates)))
            adata.obsp["spatial_connectivities"] = spatial_connectivities.tocsr()


        # Filter by selected samples/features
        if "selected" in adata.obs.columns:
            adata = adata[observations["selected"].astype(bool), :]
        if "selected" in adata.var.columns:
            adata = adata[:, features["selected"].astype(bool)]

        if dim_red is not None:
            adata.obsm["reduced_dimensions"] = (
                pd.read_csv(dim_red, index_col=0).loc[adata.obs_names].to_numpy()
            )

        if image is not None:
            adata.uns["image"] = np.array(Image.open(image))
        else:
            adata.uns["image"] = None

        return adata
    
    adata = get_anndata(coordinates, features, observations, matrix, neighbors, dim_red, image)
    adata.var_names_make_unique()

    # Seeding for controlling randomization
    random.seed(seed)
    torch.manual_seed(seed)
    np.random.seed(seed)

    with tempfile.TemporaryDirectory() as tmpdir:
        gitdir = '../../../src/DeepST/deepst'

        sys.path.append(gitdir)
        print(os.getcwd())
        import adj as adj
        
        spec = importlib.util.spec_from_file_location("deepST", f"{gitdir}/DeepST.py")
        deepST = importlib.util.module_from_spec(spec)
        sys.modules["deepST"] = deepST
        spec.loader.exec_module(deepST)
        
        deepen = deepST.run(
            save_path = None,
            task = "Identify_Domain", 
            pre_epochs = config["pre_epochs"],  
            epochs = config["epochs"],  
            use_gpu = config["use_gpu"],
        )
        
        adj_pre = adata.obsp["spatial_connectivities"]
        adj_pre = sp.coo_matrix(adj_pre)
        adj_pre = adj_pre - sp.dia_matrix((adj_pre.diagonal()[np.newaxis, :], [0]), shape=adj_pre.shape)
        adj_pre.eliminate_zeros()
        
        graph = adj.graph(data = adata, k = 1)
        adj_norm = graph.pre_graph(adj_pre)
        adj_label = adj_pre + sp.eye(adj_pre.shape[0])
        adj_label = torch.FloatTensor(adj_label.toarray())
        norm = adj_pre.shape[0] * adj_pre.shape[0] / float((adj_pre.shape[0] * adj_pre.shape[0] - adj_pre.sum()) * 2)

        graph_dict = {
                     "adj_norm": adj_norm,
                     "adj_label": adj_label,
                     "norm_value": norm 
        }

        data = adata.X.astype(np.float64)
        
        deepst_embed = deepen._fit(data, graph_dict=graph_dict, Conv_type=config["Conv_type"])
        adata.obsm["DeepST_embed"] = deepst_embed
        adata = deepen._get_cluster_data(adata, n_domains=n_clusters, priori=True)
        
        label_df = adata.obs[["DeepST_refine_domain"]]
        embedding_df = pd.DataFrame(adata.obsm['DeepST_embed'])
        embedding_df.index = adata.obs_names
        
        out_dir.mkdir(parents=True, exist_ok=True)
        label_df.columns = ["label"]
        label_df.to_csv(label_file, sep="\t", index_label="")
        
        if embedding_df is not None:
            embedding_df.to_csv(embedding_file, sep="\t", index_label="")

    return label_df