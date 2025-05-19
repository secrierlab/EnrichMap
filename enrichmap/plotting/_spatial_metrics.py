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
    metric: str = "Moran's I",
    figsize: tuple[int, int] = (4, 4),
    save=None
) -> None:
    """
    Compute and visualise spatial metrics for different scoring methods in a given dataset.

    This function computes a specified spatial metric (e.g., Moran's I) for different scoring methods
    stored in the `adata.obs` and visualises the results using a bar plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing spatial and gene expression information. 
        The `methods` specified must be columns in `adata.obs`.

    score_keys : sequence of str
        A list of method names for which to compute the spatial metric. These methods should correspond to
        columns in `adata.obs`.

    metric : str, optional
        The spatial metric to compute, e.g., "Moran's I" or "Geary's C". Defaults to "Moran's I".

    figsize : tuple of int, optional
        The size of the figure to be generated for the bar plot. Defaults to (4, 4).

    save : str or None, optional
        If specified, the plot will be saved to the file with the given filename (including extension). 
        If None, the plot will not be saved.

    Returns
    -------
    None
        Computes the spatial metric for each method and displays a bar plot with the results.
    """
    results = []

    for method in score_keys:
        if method in adata.obs.columns:
            print(f"Computing {metric} for {method}...")
            metrics, _, _, _ = compute_spatial_metrics(adata, method)
            results.append((method, metrics[metric]))

    # Convert results into a DataFrame
    df = pd.DataFrame(results, columns=["Method", metric])

    # Plotting
    plt.figure(figsize=figsize)
    sns.barplot(data=df, x="Method", y=metric, palette="muted", edgecolor=None)
    
    plt.axhline(0, linestyle="--", color="gray", linewidth=1)
    plt.title(f"{metric}", fontsize=8)
    plt.ylabel("Continuity score", fontsize=6)
    plt.xlabel("Scoring methods", fontsize=6)
    plt.yticks(fontsize=6)
    plt.xticks(rotation=45, fontsize=6)
    plt.grid(False)

    if save:
        # Ensure 'figures/' directory exists
        os.makedirs("figures", exist_ok=True)

        # If 'save' has no directory path, prepend 'figures/'
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()