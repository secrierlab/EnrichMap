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

plt.rcParams["pdf.fonttype"] = "truetype"


def gene_contributions_heatmap(
    adata: AnnData,
    score_key: str,
    top_n_genes: int = 10,
    cluster_genes: bool = True,
    order_spots: Literal["spatial", "cluster"] | None = "spatial",
    spatial_key: str = "spatial",
    cmap: str = "Reds",
    fontsize: int = 8,
    center: int = 0,
    batch_key: str | None = None,
    library_id: str | list | None = None,
    ncols: int = 2,
    figsize: tuple = (12, 6),
    save: str | None = None,
) -> None:
    """
    Visualise gene contributions from spatial gene set scoring as heatmaps,
    optionally clustered by contribution profile similarity and ordered by
    spatial proximity.

    The top contributing genes (by mean absolute contribution) are selected
    first. Their ordering and the ordering of spots along the x-axis are
    then determined by the ``cluster_genes`` and ``order_spots`` parameters:

    - **Gene axis (rows):** when ``cluster_genes=True``, genes are
      hierarchically clustered by the similarity of their spot-level
      contribution profiles, so genes that co-localise in the same tissue
      regions appear adjacent. This captures both expression magnitude and
      spatial pattern in a single ordering. When ``False``, genes are
      sorted by descending mean contribution (original behaviour).

    - **Spot axis (columns):** when ``order_spots="spatial"``, spots are
      ordered by their position along the first principal component of the
      spatial coordinates, giving a linear sweep across the tissue that
      makes spatial domains visible as contiguous colour blocks. When
      ``"cluster"``, spots are hierarchically clustered by their gene
      contribution profiles, grouping spots with similar enrichment
      patterns regardless of physical location. When ``None``, spots
      retain their original order.

    If ``batch_key`` is provided, one subplot is generated per batch (or
    per specified ``library_id``), enabling batch-wise comparison.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Gene contributions must be stored in
        ``adata.uns["gene_contributions"]`` (populated by
        ``enrichmap.tools.score``).

    score_key : str
        Name of the gene signature as used in ``enrichmap.tools.score``.

    top_n_genes : int, default 10
        Number of top contributing genes to display, selected by mean
        absolute contribution.

    cluster_genes : bool, default True
        Whether to hierarchically cluster genes (rows) by the similarity
        of their contribution profiles. When ``True``, genes that share
        spatial co-localisation patterns appear adjacent in the heatmap.
        When ``False``, genes are sorted by descending mean contribution.

    order_spots : ``"spatial"`` | ``"cluster"`` | None, default ``"spatial"``
        How to order spots (columns) in the heatmap:

        - ``"spatial"``: order by the first principal component of the
          spatial coordinates, producing a linear sweep across the tissue.
        - ``"cluster"``: hierarchical clustering of spots by their gene
          contribution profiles.
        - ``None``: retain the original spot order.

    spatial_key : str, default ``"spatial"``
        Key in ``adata.obsm`` containing the 2D spatial coordinates. Only
        used when ``order_spots="spatial"``.

    cmap : str, default ``"Reds"``
        Colourmap for the heatmap.

    fontsize : int, default 8
        Font size for labels and titles.

    center : int, default 0
        Value at which to centre the colourmap.

    batch_key : str or None, optional
        Column in ``adata.obs`` identifying batches or libraries. When
        provided, one subplot is generated per batch.

    library_id : str, list of str, or None, optional
        Specific library IDs to include. If ``None``, all batches are
        shown.

    ncols : int, default 2
        Number of subplot columns when ``batch_key`` is set.

    figsize : tuple, default (12, 6)
        Figure size for the single-panel case. For multi-panel, the size
        is computed automatically from the number of batches and ``ncols``.

    save : str or None, optional
        Path to save the figure (e.g. ``"figures/heatmap.pdf"``).

    Returns
    -------
    None
        Displays heatmaps and optionally saves the plot.

    Examples
    --------
    Default: clustered genes, spatially ordered spots:

    >>> gene_contributions_heatmap(adata, score_key="EMT")

    Mean-sorted genes, no spot reordering (original behaviour):

    >>> gene_contributions_heatmap(
    ...     adata, score_key="EMT",
    ...     cluster_genes=False, order_spots=None,
    ... )

    Per-batch comparison with spot clustering:

    >>> gene_contributions_heatmap(
    ...     adata, score_key="EMT",
    ...     batch_key="library_id",
    ...     order_spots="cluster",
    ... )
    """
    if "gene_contributions" not in adata.uns:
        raise ValueError(
            "Gene contributions not found in adata.uns. "
            "Run enrichmap.tools.score first."
        )

    contribution_matrix = adata.uns["gene_contributions"][score_key]

    # Helpers

    def _select_top_genes(contributions, n):
        """Select top genes by mean absolute contribution."""
        mean_abs = {
            gene: np.mean(np.abs(scores)) for gene, scores in contributions.items()
        }
        return sorted(mean_abs, key=mean_abs.get, reverse=True)[:n]

    def _cluster_order(matrix):
        """Hierarchical clustering order for rows of a matrix."""
        if matrix.shape[0] < 2:
            return np.arange(matrix.shape[0])
        dist = pdist(matrix, metric="correlation")
        # Guard against NaN distances (constant rows)
        dist = np.nan_to_num(dist, nan=0.0)
        Z = linkage(dist, method="average")
        return leaves_list(Z)

    def _spatial_order(coords):
        """Order spots by first PC of spatial coordinates."""
        if coords.shape[0] < 2:
            return np.arange(coords.shape[0])
        pc1 = PCA(n_components=1).fit_transform(coords).ravel()
        return np.argsort(pc1)

    def _build_heatmap(contributions, genes, spot_mask, coords=None):
        """
        Build the heatmap matrix with gene and spot ordering applied.

        Returns the (n_genes, n_spots) matrix and the ordered gene labels.
        """
        data = np.array([contributions[g] for g in genes])

        # Gene ordering
        if cluster_genes:
            gene_order = _cluster_order(data)
        else:
            gene_order = np.arange(len(genes))
        ordered_genes = [genes[i] for i in gene_order]
        data = data[gene_order]

        # Spot ordering
        if order_spots == "spatial" and coords is not None:
            spot_order = _spatial_order(
                coords[spot_mask] if spot_mask is not None else coords
            )
            data = data[:, spot_order]
        elif order_spots == "cluster":
            spot_order = _cluster_order(data.T)
            data = data[:, spot_order]

        return data, ordered_genes

    def _draw_heatmap(data, genes, ax, title):
        """Render a single heatmap panel."""
        sns.heatmap(
            data,
            yticklabels=genes,
            cmap=cmap,
            center=center,
            annot=False,
            xticklabels=False,
            ax=ax,
        )
        ax.set_xlabel("Spots", fontsize=fontsize)
        ax.set_ylabel("Top contributing genes", fontsize=fontsize)
        ax.tick_params(axis="both", labelsize=fontsize)
        ax.set_title(title, fontsize=fontsize)
        ax.grid(False)

    # Single panel

    if batch_key is None:
        genes = _select_top_genes(contribution_matrix, top_n_genes)
        data, ordered_genes = _build_heatmap(
            contribution_matrix,
            genes,
            spot_mask=None,
            coords=adata.obsm.get(spatial_key),
        )

        fig, ax = plt.subplots(figsize=figsize)
        _draw_heatmap(data, ordered_genes, ax, title=score_key)

    # Multi-panel by batch

    else:
        if batch_key not in adata.obs.columns:
            raise ValueError(f"Batch key '{batch_key}' not found in adata.obs.")

        if library_id is not None:
            if isinstance(library_id, str):
                library_id = [library_id]
            missing = [
                lid for lid in library_id if lid not in adata.obs[batch_key].values
            ]
            if missing:
                raise ValueError(
                    f"Library IDs {missing} not found in adata.obs['{batch_key}']."
                )
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
            batch_mask = (adata.obs[batch_key] == batch).values

            batch_contributions = {
                gene: scores[batch_mask] for gene, scores in contribution_matrix.items()
            }

            genes = _select_top_genes(batch_contributions, top_n_genes)
            data, ordered_genes = _build_heatmap(
                batch_contributions,
                genes,
                spot_mask=batch_mask,
                coords=adata.obsm.get(spatial_key),
            )

            _draw_heatmap(
                data,
                ordered_genes,
                axes[i],
                title=f"{score_key} ({batch})",
            )

        for j in range(n_batches, len(axes)):
            axes[j].axis("off")

    # Save

    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")

    plt.show()
