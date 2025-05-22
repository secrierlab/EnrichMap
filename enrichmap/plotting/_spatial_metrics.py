from __future__ import annotations

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from ..tools._compute_spatial_metrics import compute_spatial_metrics
from anndata import AnnData
from collections.abc import Sequence

plt.rcParams["pdf.fonttype"] = "truetype"


def spatial_metrics(
    adata: AnnData,
    score_keys: Sequence[str],
    metric: str = "Moran's I" or "Geary's C",
    n_neighbors: int = 6,
    figsize: tuple[int, int] = (4, 4),
    save=None
) -> None:
    """
    Compute and visualise spatial metrics for different scoring methods in a given dataset.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing spatial and gene expression information.
    score_keys : sequence of str
        A list of method names for which to compute the spatial metric. These methods should correspond to
        columns in `adata.obs`.
    metric : str, optional
        The spatial metric to compute, e.g., "Moran's I" or "Geary's C". Defaults to "Moran's I".
    n_neighs : int, optional
        Number of neighbours to use when computing spatial weights. Defaults to 6.
    figsize : tuple of int, optional
        The size of the figure to be generated for the bar plot. Defaults to (4, 4).
    save : str or None, optional
        If specified, the plot will be saved to the file with the given filename. If None, the plot will not be saved.

    Returns
    -------
    None
    """
    results = []

    for method in score_keys:
        if method in adata.obs.columns:
            print(f"Computing {metric} for {method} with {n_neighbors} neighbours...")
            metrics, _ = compute_spatial_metrics(adata, score_key=method, n_neighbors=n_neighbors)
            results.append((method, metrics[metric]))

    df = pd.DataFrame(results, columns=["Method", metric])

    plt.figure(figsize=figsize)
    sns.barplot(data=df, x="Method", y=metric, palette="muted", edgecolor=None)

    plt.axhline(0, linestyle="--", color="gray", linewidth=1)
    plt.title(f"{metric}", fontsize=8)
    if metric == "Moran's I":
        plt.ylabel("Spatial coherence", fontsize=6)
    elif metric == "Geary's C":
        plt.ylabel("Local dissimilarity", fontsize=6)
    plt.xlabel("Scoring methods", fontsize=6)
    plt.yticks(fontsize=6)
    plt.xticks(rotation=90, fontsize=6)
    plt.grid(False)

    if save:
        os.makedirs("figures", exist_ok=True)
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()
