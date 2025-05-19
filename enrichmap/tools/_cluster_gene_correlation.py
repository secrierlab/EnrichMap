from __future__ import annotations
from typing import TYPE_CHECKING

from anndata import AnnData
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd

from scipy.cluster.hierarchy import linkage

if TYPE_CHECKING:
    from typing import Literal

plt.rcParams["pdf.fonttype"] = "truetype"

def cluster_gene_correlation(
    adata: AnnData,
    signature_name: str = "enrichmap",
    top_n_genes: int = 10,
    method: Literal["single", "complete", "average", "weighted", "centroid", "median", "ward"] = 'ward',
    metric: Literal["euclidean", "correlation", "cosine"] = 'euclidean',
    figsize: tuple = (10, 10),
    cmap: str = "seismic",
    save_fig: str | None = None
) -> None:
    """
    Visualise clustered gene-gene correlation heatmap based on gene contributions to a signature.

    Computes pairwise gene correlations, constructs a distance matrix, performs hierarchical clustering, 
    and plots a clustered heatmap of the top contributing genes.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must contain `gene_contributions` for the given `signature_name` in `adata.uns`.

    signature_name : str, default "enrichmap"
        Name of the gene signature whose gene contributions are used to compute the correlation matrix.

    top_n_genes : int, default 10
        Number of top genes (based on overall correlation strength) to include in the heatmap. 
        Set to None to use all genes.

    method : str, default 'ward'
        Linkage method for hierarchical clustering. One of 'single', 'complete', 'average', 'weighted', 
        'centroid', 'median', 'ward'.

    metric : str, default 'euclidean'
        Distance metric for linkage computation. One of 'euclidean', 'correlation', 'cosine'.

    figsize : tuple, default (10, 10)
        Figure size of the clustered heatmap.

    cmap : str, default "seismic"
        Colormap used for the heatmap.

    save_fig : str, optional
        File path to save the figure. If None, the plot is only shown.

    Returns
    -------
    None
        Displays the clustered heatmap.

    """
    contribution_matrix = adata.uns['gene_contributions'][signature_name]
    correlation_matrix = pd.DataFrame(contribution_matrix).corr()

    if top_n_genes is not None:
        top_genes = correlation_matrix.abs().sum(axis=1).nlargest(top_n_genes).index
        correlation_matrix = correlation_matrix.loc[top_genes, top_genes]

    # Ensure the correlation matrix is symmetric and fill NaNs with 0
    correlation_matrix = correlation_matrix.fillna(0)
    correlation_matrix = (correlation_matrix + correlation_matrix.T) / 2
    np.fill_diagonal(correlation_matrix.values, 1)

    # Convert correlation to a valid distance matrix
    dist_matrix = (1 - correlation_matrix) / 2  # Ensures values are in [0, 1]

    # Ensure no NaN or Inf values
    if np.any(np.isnan(dist_matrix)) or np.any(np.isinf(dist_matrix)):
        raise ValueError("Distance matrix contains NaN or infinite values after transformation.")

    # Convert to condensed format for linkage
    condensed_dist = ssd.pdist(dist_matrix)

    # Compute hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=method, metric=metric)

    # Generate cluster map
    sns.clustermap(correlation_matrix, row_linkage=linkage_matrix, col_linkage=linkage_matrix, figsize=figsize, cmap=cmap, square=True, vmin=0, vmax=1, dendrogram_ratio=(0.2, 0.2), cbar_kws={"shrink": 0.5})
    plt.title(f"Clustered gene correlation heatmap - {signature_name}")
    plt.grid(False)
    if save_fig is not None:
        plt.savefig(save_fig, dpi=300)
    plt.show()