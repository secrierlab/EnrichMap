from __future__ import annotations

import squidpy as sq
import numpy as np

from anndata import AnnData
from libpysal.weights import KNN
from esda.moran import Moran, Moran_Local
from esda.geary import Geary
from esda.getisord import G_Local
from libpysal.weights import KNN

def compute_spatial_metrics(
    adata: AnnData,
    score_key: str = "enrichment_score",
    n_neighs: int = 6,
):
    """
    Compute global and local spatial autocorrelation metrics for a spatial score.

    This function calculates several spatial statistics to assess spatial autocorrelation in gene set
    enrichment or similar scores stored in `adata.obs[score_key]`. It returns Moran's I (global),
    Geary's C (local variance), Local Moran's I (spatial clusters), and Getis-Ord G (hotspots).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm['spatial']`. If spatial neighbours
        are not precomputed, they will be computed using Squidpy.

    score_key : str, default "enrichment_score"
        Column in `adata.obs` containing the score vector for which spatial metrics should be computed.

    Returns
    -------
    metrics : dict
        Dictionary containing:
            - "Moran's I": float
            - "Geary's C": float
            - "Local Moran I": np.ndarray of local Moran’s I values
            - "Local G": np.ndarray of Getis-Ord G Z-scores

    W : libpysal.weights.W
        Spatial weights matrix used for computation.

    local_moran : esda.Moran_Local
        Full Local Moran’s I result object.

    local_g : esda.G_Local
        Full Getis-Ord G result object.

    Notes
    -----
    - Spatial neighbours are computed using 6 nearest neighbours by default if not already present.
    - The function expects coordinates in `adata.obsm["spatial"]` and computes weights using libpysal.
    - All values are computed using the `esda` and `libpysal` libraries.
    """
    # Ensure spatial neighbors are computed
    if "spatial_neighbors" not in adata.uns:
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type="generic", key_added="spatial")
    
    # Extract spatial weights
    W = KNN.from_array(adata.obsm["spatial"])
    W.transform = "r"

    # Extract score
    score = adata.obs[score_key].values

    # Compute Moran’s I (global spatial autocorrelation)
    moran = Moran(score, W)
    
    # Compute Geary’s C (measuring local variations)
    geary = Geary(score, W)

    # Compute Local Moran’s I (spatial clusters/outliers)
    score = score.astype(np.float64)
    local_moran = Moran_Local(score, W)

    # Compute Getis-Ord G (hotspot analysis)
    local_g = G_Local(score, W)

    # Collect results
    metrics = {
        "Moran's I": moran.I,
        "Geary's C": geary.C,
        "Local Moran I": local_moran.Is,
        "Local G": local_g.Zs
    }

    return metrics, W, local_moran, local_g