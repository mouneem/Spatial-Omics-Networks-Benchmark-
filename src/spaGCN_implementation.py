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


from PIL import Image
import scipy as sp
import anndata as ad

import json
import random
import numpy as np
import pandas as pd
import SpaGCN as spg
import torch
from pathlib import Path


def run_spagcn(coordinates, features, observations, out_dir, n_clusters, technology, seed, matrix=None, neighbors=None, dim_red=None, image=None, config=None):

    # Load configuration

    if config["refine"] and technology not in ["Visium", "ST"]:
        raise Exception(
            f"Invalid parameter combination. Refinement only works with Visium and ST not {technology}"
        )

    # Set beta (TODO: determine this dynamically if needed)
    beta = 49

    # Set seeds for reproducibility
    random.seed(seed)
    torch.manual_seed(seed)
    np.random.seed(seed)

    # Create output directory
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Output files
    label_file = out_dir / "domains.tsv"
    embedding_file = out_dir / "embedding.tsv"

    observations_df = pd.read_csv(observations, index_col=0)


    features_df = pd.read_csv(features, index_col=0)
    X = features_df
    if sp.sparse.issparse(X):
        X = X.tocsr()

    coordinates_df = pd.read_csv(coordinates, index_col=0).loc[observations_df.index, :].to_numpy()
    print(coordinates_df.shape)
    print(features_df.shape)

    adata = ad.AnnData(
        X=X, obs=observations_df, var=features_df.T, obsm={"spatial_pixel": coordinates_df}
    )
    adata.uns["image"] = None

    if image is not None:
        adata.uns["image"] = np.array(Image.open(image))
    else:
        adata.uns["image"] = None

    # Ensure 'col' and 'row' exist in observations
    if 'col' not in adata.obs.columns or 'row' not in adata.obs.columns:
        adata.obs['col'] = coordinates_df[:, 0]
        adata.obs['row'] = coordinates_df[:, 1]
    

    # Calculate adjacency matrix
    if technology in ["Visium", "ST"]:
        if adata.uns["image"] is not None:
            adj = spg.calculate_adj_matrix(
                adata.obs["col"],
                adata.obs["row"],
                adata.obsm["spatial_pixel"][:, 0],
                adata.obsm["spatial_pixel"][:, 1],
                image=adata.uns["image"],
                alpha=config["alpha"],
                beta=beta,
                histology=True,
            )
        else:
            adj = spg.calculate_adj_matrix(
                adata.obs["col"], adata.obs["row"], histology=False
            )
    else:
        adj = spg.calculate_adj_matrix(
            adata.obsm["spatial_pixel"][:, 0],
            adata.obsm["spatial_pixel"][:, 1],
            histology=False,
        )

    clf = spg.SpaGCN()

    # Find the l value given p
    l = spg.search_l(config["p"], adj)
    clf.set_l(l)

    # Determine number of principal components
    num_samples, num_features = adata.shape
    n_pcs = min(config["n_pcs"], num_samples, num_features)

    # Train model

    clf.train(
        adata,
        adj,
        init_spa=True,
        init=config["method"],
        n_clusters=n_clusters,
        num_pcs=n_pcs,
    )
    
    y_pred, prob = clf.predict()
    adata.obs["cluster"] = pd.Series(y_pred, index=adata.obs_names, dtype="category")

    if technology in ["Visium", "ST"] and config["refine"]:
        adj_2d = spg.calculate_adj_matrix(
            adata.obs["col"], adata.obs["row"], histology=False
        )
        shape = "hexagon" if technology == "Visium" else "square"
        refined_pred = spg.refine(
            sample_id=adata.obs_names.tolist(),
            pred=adata.obs["cluster"].tolist(),
            dis=adj_2d,
            shape=shape,
        )
        adata.obs["refined_cluster"] = pd.Series(refined_pred, index=adata.obs_names, dtype="category")
        label_df = adata.obs[["refined_cluster"]]
    else:
        label_df = adata.obs[["cluster"]]

    # Write output
    label_df.columns = ["label"]
    label_df.to_csv(label_file, sep="\t", index_label="")
    return label_df