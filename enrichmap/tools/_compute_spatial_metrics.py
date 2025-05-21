from __future__ import annotations

import squidpy as sq
import numpy as np

from anndata import AnnData
from libpysal.weights import KNN
from esda.moran import Moran, Moran_Local
from esda.geary import Geary
from esda.getisord import G_Local


def compute_spatial_metrics(
    adata: AnnData,
    score_key: str = "enrichment_score",
    n_neighbors: int = 6,
):
    """
    Compute global and local spatial autocorrelation metrics for a spatial score.

    This function calculates several spatial statistics to assess spatial autocorrelation in gene set
    enrichment or similar scores stored in `adata.obs[score_key]`. It returns Moran's I (global),
    Geary's C (local variance), Local Moran's I (spatial clusters), and Getis-Ord G (hotspots).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm['spatial']`.
    score_key : str, default "enrichment_score"
        Column in `adata.obs` containing the score vector for which spatial metrics should be computed.
    n_neighbors : int, default 6
        Number of spatial neighbours to use for spatial weights calculation.

    Returns
    -------
    metrics : dict
        Dictionary containing:
            - "Moran's I": float
            - "Geary's C": float
            - "Local Moran I": np.ndarray of local Moran’s I values

    W : libpysal.weights.W
        Spatial weights matrix used for computation.

    local_moran : esda.Moran_Local
        Full Local Moran’s I result object.
    """
    # Compute spatial neighbors and weights
    sq.gr.spatial_neighbors(adata, n_neighs=n_neighbors, coord_type="generic", key_added="spatial")

    # Extract spatial weights from coordinates
    W = KNN.from_array(adata.obsm["spatial"])
    W.transform = "r"

    # Extract score vector
    score = adata.obs[score_key].values.astype(np.float64)

    # Compute global Moran's I
    moran = Moran(score, W)

    # Compute Geary's C
    geary = Geary(score, W)

    # Compute Local Moran's I
    local_moran = Moran_Local(score, W)

    metrics = {
        "Moran's I": moran.I,
        "Geary's C": geary.C,
        "Local Moran I": local_moran.Is,
    }

    return metrics, W
