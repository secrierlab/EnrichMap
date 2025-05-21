from __future__ import annotations

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from scipy.stats import pearsonr
from libpysal.weights import KNN
from libpysal.weights.spatial_lag import lag_spatial

plt.rcParams["pdf.fonttype"] = "truetype"

def cross_moran_scatter(
    adata: AnnData,
    score_x: str,
    score_y: str,
    batch_key: str | None = None,
    n_neighbours: int = 6,
    save: str | None = None
):
    """
    Plot cross-Moran scatterplots: score_x vs spatial lag of score_y, per sample or globally.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm["spatial"]`.
    score_x : str
        Column in `adata.obs` for X-axis variable.
    score_y : str
        Column in `adata.obs` for which the spatial lag is computed (Y-axis).
    batch_key : str or None
        Key in `adata.obs` indicating batch/sample/library identifiers.
        If None, computes a single plot using the entire dataset.
    n_neighbours : int
        Number of nearest neighbours for spatial weights.
    save : str or None
        If provided, saves the figure to the given file name.
    """
    if batch_key is None:
        fig, ax = plt.subplots(figsize=(3, 3))
        batches = [adata]
        titles = ["Cross-Moran scatter plot"]
    else:
        batches = [adata[adata.obs[batch_key] == b].copy() for b in adata.obs[batch_key].unique()]
        titles = adata.obs[batch_key].unique().tolist()
        n = len(batches)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), constrained_layout=True)
        axes = np.ravel(axes)

    for i, (ad, title) in enumerate(zip(batches, titles)):
        coords = ad.obsm["spatial"]
        x = ad.obs[score_x].values
        y = ad.obs[score_y].values

        W = KNN.from_array(coords, k=n_neighbours)
        W.transform = "r"
        y_lag = lag_spatial(W, y)

        mask = ~np.isnan(x) & ~np.isnan(y_lag)
        if np.sum(mask) > 1:
            r, p = pearsonr(x[mask], y_lag[mask])
        else:
            r, p = np.nan, np.nan

        ax = ax if batch_key is None else axes[i]
        ax.scatter(x, y_lag, s=10, alpha=0.3, color="lightblue")
        sns.regplot(x=x, y=y_lag, scatter=False, ax=ax, color="black", line_kws={"lw": 1})
        ax.axhline(0, color="grey", lw=1)
        ax.axvline(0, color="grey", lw=1)
        ax.set_title(f"Slide {title}\nr = {r:.2f}, p = {p:.2g}", fontsize=10)
        ax.set_xlabel(score_x)
        ax.set_ylabel(f"Spatial lag of {score_y}")
        ax.grid(False)

    if batch_key is not None:
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")

    if save:
        os.makedirs("figures", exist_ok=True)
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")

    plt.show()
