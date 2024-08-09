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

import json
import tempfile
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import torch
from scipy.io import mmread
from scipy.sparse import issparse
from pathlib import Path
from popari.components import PopariDataset
from popari.io import save_anndata
from popari.model import SpiceMix
from popari import tl

def get_anndata(coordinates, matrix, features, observations, neighbors=None, config=None):
    """Convert data input into SpiceMix anndata format."""
    counts_matrix = mmread(matrix)
    if issparse(counts_matrix):
        counts_matrix = counts_matrix.tocsr()

    observations = pd.read_csv(observations, delimiter="\t", index_col=0)
    features = pd.read_table(features, index_col=0)

    coordinates = (
        pd.read_table(coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adjacency_matrix = mmread(neighbors).T.tocsr() if neighbors else None

    adata = ad.AnnData(
        X=counts_matrix,
        obs=observations,
        obsm={"spatial": coordinates},
        obsp={"spatial_connectivities": adjacency_matrix} if adjacency_matrix is not None else None,
    )

    preprocess_parameters = config.pop("preprocess", None)
    if preprocess_parameters:
        if adjacency_matrix is not None:
            del adata.obsp["spatial_connectivities"]
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=preprocess_parameters["hvgs"])

        adata = PopariDataset(adata[:, adata.var["highly_variable"]], "processed")
        adata.compute_spatial_neighbors()

    else:

        if adjacency_matrix is not None:
            transpose_condition = adjacency_matrix.T > adjacency_matrix
            adjacency_matrix = (
                adjacency_matrix
                + adjacency_matrix.T.multiply(transpose_condition)
                - adjacency_matrix.multiply(transpose_condition)
            )

            num_cells = adjacency_matrix.shape[0]
            adjacency_list = [[] for _ in range(num_cells)]
            for x, y in zip(*adjacency_matrix.nonzero()):
                adjacency_list[x].append(y)

            adata.obsp["adjacency_matrix"] = adjacency_matrix
            adata.obs["adjacency_list"] = adjacency_list

            del adata.obsp["spatial_connectivities"]

        adata = PopariDataset(adata, "processed")

    return adata


def run_spicemix(
    coordinates,
    matrix,
    features,
    observations,
    out_dir,
    n_clusters,
    technology,
    seed,
    neighbors=None,
    dim_red=None,
    image=None,
    config=None,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    label_file = out_dir / "domains.tsv"

    if config is not None:
        with open(config, "r") as f:
            config = json.load(f)
    else:
        config = {}

    preprocess_parameters = config.pop("preprocess", None)

    adata = get_anndata(coordinates, matrix, features, observations, neighbors)

    num_preiterations = config.pop("num_preiterations", 5)
    num_iterations = config.pop("num_iterations", 200)

    device = config.pop("device", "cpu") if torch.cuda.is_available() else 'cpu'
    dtype = config.pop("dtype", "float32")

    if dtype == "float32":
        dtype = torch.float32
    elif dtype == "float64":
        dtype = torch.float64

    torch_context = {
        "device": device,
        "dtype": dtype,
    }

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        save_anndata(temp_path / "input.h5ad", [adata])

        model = SpiceMix(
            dataset_path=temp_path / "input.h5ad",
            random_state=seed,
            initial_context=torch_context,
            torch_context=torch_context,
            **config,
        )

        # Run
        model.train(num_preiterations=num_preiterations, num_iterations=num_iterations)

        tl.preprocess_embeddings(model, normalized_key="normalized_X")
        tl.leiden(
            model, joint=True, target_clusters=n_clusters, use_rep="normalized_X"
        )

        # TODO: add optional smoothing step
        label_df = model.datasets[0].obs[["leiden"]]

    ## Write output
    label_df.columns = ["label"]
    label_df.to_csv(label_file, sep="\t", index_label="")

# Example usage
run_spicemix(
    coordinates="path/to/coordinates.tsv",
    matrix="path/to/matrix.mtx",
    features="path/to/features.tsv",
    observations="path/to/observations.tsv",
    out_dir="path/to/output/directory",
    n_clusters=10,
    technology="Visium",
    seed=42,
    neighbors="path/to/neighbors.mtx",
    config="path/to/config.json"
)
