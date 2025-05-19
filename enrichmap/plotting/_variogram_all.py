from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from anndata import AnnData
from typing import List
from scipy.spatial.distance import pdist
from skgstat import Variogram

plt.rcParams["pdf.fonttype"] = "truetype"

def variogram_all(
    adata: AnnData,
    score_keys: List[str],
    save: None | str = None,
    max_lag: float | None = None,
    lag_percentile: float = 95
) -> None:
    """
    Compute and plot empirical variograms for multiple score keys on the same plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm["spatial"]`.
    score_keys : list of str
        List of keys in `adata.obs` to compute variograms for.
    save : str or None, optional
        If provided, path to save the figure as a PDF file.
    max_lag : float or None, optional
        If set, limits the x-axis of the variogram plot to this value.
        If None, computed from pairwise distances using `lag_percentile`.
    lag_percentile : float, optional
        Percentile of pairwise distances to set as max_lag when max_lag is None (default: 95).
    """
    coords = adata.obsm["spatial"]

    if max_lag is None:
        dists = pdist(coords)
        max_lag = np.percentile(dists, lag_percentile)

    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)

    colours = plt.cm.tab10(np.linspace(0, 1, len(score_keys)))
    variograms = []

    for key, colour in zip(score_keys, colours):
        values = adata.obs[key].values
        V = Variogram(coords, values, method="cressie", model="gaussian", verbose=False)
        variograms.append(V)

        ax.plot(V.bins, V.experimental, "o-", label=key, color=colour)
        ax.axhline(np.var(values), color=colour, linestyle="--", alpha=0.5)

    ax.set_xlabel("Spatial lag")
    ax.set_ylabel("Semivariance")
    ax.set_title("Empirical variograms")
    ax.legend(title="Signatures")
    ax.grid(False)
    ax.set_xlim(0, max_lag)

    if save:
        # Ensure 'figures/' directory exists
        os.makedirs("figures", exist_ok=True)

        # If 'save' has no directory path, prepend 'figures/'
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()