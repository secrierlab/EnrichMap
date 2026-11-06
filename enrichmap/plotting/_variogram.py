from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from anndata import AnnData
from typing import List

from scipy.spatial.distance import pdist
from skgstat import Variogram

plt.rcParams["pdf.fonttype"] = "truetype"


def variogram(
    adata: AnnData,
    score_key: List[str],
    save: None | str = None,
    max_lag: float | None = None,
    lag_percentile: float = 95,
) -> None:
    """
    Compute empirical variograms to assess spatial dependence for any score key.

    Semivariance is shown divided by the variance of the score, so scores on different scales share
    one axis. The dashed line at 1 is the variance: values below it mean that spots at that lag are
    more alike than spots picked at random.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm["spatial"]`.
    score_key : list of str
        List of keys in `adata.obs` to compute variograms for.
    save : str or None, optional
        If provided, path to save the figure as a PDF file.
    max_lag : float or None, optional
        If set, limits the x-axis of the variogram plot to this value.
        If None, computed from pairwise distances using `lag_percentile`.
    lag_percentile : float, optional
        Percentile of pairwise distances to set as max_lag when max_lag is None (default: 95).

    Returns
    -------
    variograms : list of Variogram
        List of computed Variogram objects for each score key.
    """
    coords = adata.obsm["spatial"]

    if max_lag is None:
        dists = pdist(coords)
        max_lag = np.percentile(dists, lag_percentile)

    n = len(score_key)
    fig, axes = plt.subplots(
        1, n, figsize=(2.5 * n, 2.1), sharey=True, constrained_layout=True
    )
    if n == 1:
        axes = [axes]

    variograms = []
    for ax, key in zip(axes, score_key):
        values = adata.obs[key].values
        V = Variogram(coords, values, method="cressie", model="gaussian", verbose=False)
        variograms.append(V)

        inside = V.bins <= max_lag
        ax.plot(
            V.bins[inside],
            V.experimental[inside] / np.var(values),
            "o-",
            ms=2,
            lw=0.8,
            color="#2f0bdf",
        )
        ax.axhline(1, color="0.5", ls="--", lw=0.5)
        ax.set_xlabel("Spatial lag", fontsize=7)
        ax.set_title(f"{key}", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.grid(False)
        ax.set_xlim(0, max_lag)
        sns.despine(ax=ax)

    axes[0].set_ylabel("Semivariance / variance", fontsize=7)
    axes[0].set_ylim(0, None)

    if save:
        # Ensure 'figures/' directory exists
        os.makedirs("figures", exist_ok=True)

        # If 'save' has no directory path, prepend 'figures/'
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()
