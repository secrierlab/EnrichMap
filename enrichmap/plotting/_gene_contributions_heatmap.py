from __future__ import annotations

import os
from typing import Literal

from anndata import AnnData
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
import squidpy as sq

plt.rcParams["pdf.fonttype"] = "truetype"


def gene_contributions_heatmap(
    adata: AnnData,
    score_key: str,
    top_n_genes: int = 5,
    bottom_n_genes: int = 5,
    cluster_genes: bool = True,
    order_spots: Literal["spatial", "cluster", "coherence"] | None = "coherence",
    spatial_key: str = "spatial",
    cmap: str = "RdBu_r",
    fontsize: int = 8,
    center: float = 0.0,
    batch_key: str | None = None,
    library_id: str | list | None = None,
    ncols: int = 2,
    figsize: tuple = (12, 6),
    save: str | None = None,
    n_neighbors: int = 6,
) -> None:
    """
    Spatial gene contribution heatmap with expression scaling and spatial ordering.

    This function visualises gene contribution matrices derived from spatial
    gene set scoring. It combines expression magnitude, spatial autocorrelation,
    and graph-based spatial coherence into a single interpretable heatmap.

    Core design:
    - Rows = genes (selected by extreme contribution: high + low)
    - Columns = spots
    - Values = z-scored gene contributions (per gene)

    The pipeline explicitly separates three ordering layers:
    1. Gene selection by mean absolute contribution
    2. Gene ordering by expression + Moran’s I spatial autocorrelation
    3. Spot ordering by spatial structure or spatial coherence

    Parameters
    ----------
    adata : AnnData
        Input AnnData containing gene contribution matrices in:
        adata.uns["gene_contributions"][score_key]

    score_key : str
        Key identifying the gene set or signature.

    top_n_genes : int
        Number of highest contributing genes to include.

    bottom_n_genes : int
        Number of lowest contributing genes to include.

    cluster_genes : bool
        Whether to hierarchically cluster genes by similarity of spatial
        contribution profiles.

    order_spots : {"spatial", "cluster", "coherence", None}
        Strategy for ordering spots:
        - spatial: PCA projection of coordinates
        - cluster: hierarchical clustering of spot profiles
        - coherence: graph-based spatial smoothness (local signal consistency)
        - None: original order

    spatial_key : str
        Key in adata.obsm containing spatial coordinates.

    cmap : str
        Colormap used for heatmap. Diverging recommended due to z-scoring.

    fontsize : int
        Font size for labels and titles.

    center : float
        Center value for colormap scaling (typically 0 for z-scores).

    batch_key : str or None
        Optional grouping variable for per-batch visualisation.

    library_id : str or list or None
        Subset of batches to plot.

    ncols : int
        Number of columns in multi-panel layout.

    figsize : tuple
        Figure size per panel.

    save : str or None
        File path to save output figure.

    n_neighbors : int
        Number of neighbours used in spatial graph construction.

    Returns
    -------
    None
        Displays heatmap figure.

    Notes
    -----
    - Expression values are z-scored per gene before plotting.
    - Diverging colormap is meaningful only after scaling.
    - Spatial coherence is computed using graph-weighted signal smoothing.
    - Moran’s I is used to encode spatial autocorrelation at gene level.
    """

    if "gene_contributions" not in adata.uns:
        raise ValueError(
            "Gene contributions not found in adata.uns. "
            "Run enrichmap.tools.score first."
        )

    contribution_matrix = adata.uns["gene_contributions"][score_key]

    def _select_extreme_genes(contributions, n_top, n_bottom):
        mean_abs = {
            gene: np.mean(np.abs(scores)) for gene, scores in contributions.items()
        }
        sorted_genes = sorted(mean_abs, key=mean_abs.get, reverse=True)
        top = sorted_genes[:n_top]
        bottom = sorted_genes[-n_bottom:] if n_bottom > 0 else []
        bottom = [g for g in bottom if g not in top]
        return top + bottom, mean_abs

    def _cluster_order(matrix):
        if matrix.shape[0] < 2:
            return np.arange(matrix.shape[0])
        dist = pdist(matrix, metric="correlation")
        dist = np.nan_to_num(dist, nan=0.0)
        Z = linkage(dist, method="average")
        return leaves_list(Z)

    def _spatial_order(coords):
        if coords is None or coords.shape[0] < 2:
            return np.arange(coords.shape[0])
        pc1 = PCA(n_components=1).fit_transform(coords).ravel()
        return np.argsort(pc1)

    def _scale_rows(X):
        """
        Z-score per gene (row-wise scaling).
        This is what makes RdBu_r meaningful.
        """
        mean = X.mean(axis=1, keepdims=True)
        std = X.std(axis=1, keepdims=True)
        std[std == 0] = 1.0
        return (X - mean) / std

    def _compute_spot_coherence(data, coords):
        """
        Local spatial coherence via graph-weighted signal smoothness.
        """
        from scipy.sparse import csr_matrix
        from anndata import AnnData as _AnnData

        signal = np.mean(np.abs(data), axis=0)

        if coords is None or len(signal) < 3:
            return np.zeros_like(signal)

        ad = _AnnData(X=np.zeros((coords.shape[0], 1)))
        ad.obsm[spatial_key] = coords

        sq.gr.spatial_neighbors(ad, coord_type="generic")
        W = ad.obsp["spatial_connectivities"]

        if not isinstance(W, csr_matrix):
            W = csr_matrix(W)

        row_sum = np.array(W.sum(axis=1)).flatten()
        row_sum[row_sum == 0] = 1.0
        W = W.multiply(1 / row_sum[:, None])

        smoothed = W.dot(signal)
        return smoothed

    def _order_genes(contributions, genes, mean_abs, coords):
        moran = _compute_gene_morans(contributions, genes, coords)
        return sorted(
            genes,
            key=lambda g: (mean_abs[g], moran.get(g, 0.0)),
            reverse=True,
        )

    def _compute_gene_morans(contributions, genes, coords):
        from anndata import AnnData as _AnnData

        X = np.array([contributions[g] for g in genes]).T
        ad = _AnnData(X=X)
        ad.obsm[spatial_key] = coords
        ad.var_names = genes

        sq.gr.spatial_neighbors(ad, coord_type="generic")
        sq.gr.spatial_autocorr(ad, mode="moran")

        return ad.uns["moranI"]["I"].to_dict()

    def _order_spots(data, coords):
        if order_spots == "spatial":
            return _spatial_order(coords)
        elif order_spots == "cluster":
            return _cluster_order(data.T)
        elif order_spots == "coherence":
            score = _compute_spot_coherence(data, coords)
            return np.argsort(score)[::-1]
        return np.arange(data.shape[1])

    def _build_heatmap(contributions, genes, mean_abs, coords):
        ordered_genes = _order_genes(contributions, genes, mean_abs, coords)

        X = np.array([contributions[g] for g in ordered_genes])

        X = _scale_rows(X)

        spot_order = _order_spots(X, coords)
        X = X[:, spot_order]

        return X, ordered_genes

    def _draw_heatmap(data, genes, ax, title):
        vmax = np.nanpercentile(np.abs(data), 99)
        vmin = -vmax

        sns.heatmap(
            data,
            yticklabels=genes,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            xticklabels=False,
            ax=ax,
        )
        ax.set_xlabel("Spots (spatial coherence)", fontsize=fontsize)
        ax.set_ylabel("Genes (z-scored)", fontsize=fontsize)
        ax.set_title(title, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        ax.grid(False)

    # Single panel
    if batch_key is None:
        coords = adata.obsm.get(spatial_key)

        genes, mean_abs = _select_extreme_genes(
            contribution_matrix, top_n_genes, bottom_n_genes
        )

        data, ordered_genes = _build_heatmap(
            contribution_matrix,
            genes,
            mean_abs,
            coords,
        )

        fig, ax = plt.subplots(figsize=figsize)
        _draw_heatmap(data, ordered_genes, ax, score_key)

    # Multi-panel
    else:
        if batch_key not in adata.obs.columns:
            raise ValueError(f"{batch_key} not found in adata.obs")

        if library_id is not None:
            if isinstance(library_id, str):
                library_id = [library_id]
            unique_batches = library_id
        else:
            unique_batches = list(adata.obs[batch_key].unique())

        n_batches = len(unique_batches)
        n_rows = (n_batches + ncols - 1) // ncols

        fig, axes = plt.subplots(
            n_rows,
            ncols,
            figsize=(ncols * figsize[0] / 2, n_rows * figsize[1] / 2),
            constrained_layout=True,
        )
        axes = np.atleast_2d(axes).flatten()

        for i, batch in enumerate(unique_batches):
            mask = (adata.obs[batch_key] == batch).values
            coords = adata.obsm[spatial_key][mask]

            batch_contributions = {g: v[mask] for g, v in contribution_matrix.items()}

            genes, mean_abs = _select_extreme_genes(
                batch_contributions, top_n_genes, bottom_n_genes
            )

            data, ordered_genes = _build_heatmap(
                batch_contributions,
                genes,
                mean_abs,
                coords,
            )

            _draw_heatmap(data, ordered_genes, axes[i], f"{score_key} ({batch})")

        for j in range(n_batches, len(axes)):
            axes[j].axis("off")

    if save:
        os.makedirs("figures", exist_ok=True)
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")

    plt.show()
