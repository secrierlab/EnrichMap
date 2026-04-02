from __future__ import annotations

import os
import numpy as np
import pandas as pd
import seaborn as sns

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from anndata import AnnData
from scipy.stats import spearmanr, pearsonr
from pathlib import Path


def compute_corr_and_pval(df, method="spearman"):
    n = df.shape[1]
    corr = np.zeros((n, n))
    pvals = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            if method == "spearman":
                r, p = spearmanr(df.iloc[:, i], df.iloc[:, j], nan_policy="omit")
            elif method == "pearson":
                r, p = pearsonr(df.iloc[:, i], df.iloc[:, j])
            else:
                raise ValueError("Unsupported method: choose 'spearman' or 'pearson'")
            corr[i, j] = r
            pvals[i, j] = p
    return (
        pd.DataFrame(corr, index=df.columns, columns=df.columns),
        pd.DataFrame(pvals, index=df.columns, columns=df.columns),
    )


def get_star_annot(pvals):
    annot = pd.DataFrame("", index=pvals.index, columns=pvals.columns)
    for i in range(pvals.shape[0]):
        for j in range(pvals.shape[1]):
            if i == j:
                continue
            p = pvals.iloc[i, j]
            if p < 0.001:
                annot.iloc[i, j] = "***"
            elif p < 0.01:
                annot.iloc[i, j] = "**"
            elif p < 0.05:
                annot.iloc[i, j] = "*"
    return annot


def signature_correlation_heatmap(
    adata: AnnData,
    score_keys: list[str],
    library_key: str | None = None,
    library_id: str | list[str] | None = None,
    method: str = "spearman",
    tile_size: float = 0.25,
    cmap: str | None = "coolwarm",
    save: str | Path | None = None,
    **kwargs: dict[str, any],
):
    """
    Plot a heatmap of correlations between gene set scores in `adata.obs`.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with gene signature scores in `adata.obs`.

    score_keys : list of str
        Column names in `adata.obs` corresponding to gene set signature scores.

    library_key : str or None
        Key in `adata.obs` for library identifiers (e.g. patient ID).

    library_id : str or None
        If set, filter `adata` for this library before plotting.

    method : str
        Correlation method: "spearman" or "pearson".

    tile_size : float
        Size of each heatmap tile in inches. The figure dimensions are
        computed from this value and the number of signatures, so tiles
        remain a consistent size regardless of how many panels are drawn.

    cmap : str or None
        Colormap for the heatmap. Defaults to "coolwarm".

    save : str or Path or None
        Path to save the figure.

    **kwargs : dict
        Additional keyword arguments passed to `sns.heatmap`.
    """
    # Optional library filtering
    if library_key and library_id is not None:
        mask = (
            adata.obs[library_key] == library_id
            if isinstance(library_id, str)
            else adata.obs[library_key].isin(library_id)
        )
        adata = adata[mask]

    n_sigs = len(score_keys)

    # Panel dimensions derived from tile size
    heatmap_size = tile_size * n_sigs
    label_margin = 1.0  # space for rotated tick labels
    title_margin = 0.35  # space for the panel title
    panel_w = heatmap_size + label_margin
    panel_h = heatmap_size + label_margin + title_margin

    def _clean_label(text):
        return text.replace("_score", "")

    def _cluster_order(corr):
        """Hierarchical clustering order from a correlation matrix."""
        dist = 1 - corr.values
        np.fill_diagonal(dist, 0)
        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        return leaves_list(Z)

    def _plot_heatmap(corr, pvals, title, ax):
        """Draw a single correlation heatmap on the given axes."""
        annot = get_star_annot(pvals)
        order = _cluster_order(corr)
        corr = corr.iloc[order, order]
        annot = annot.iloc[order, order]

        sns.heatmap(
            corr,
            annot=annot,
            fmt="",
            cmap=cmap,
            square=True,
            ax=ax,
            vmin=-1,
            vmax=1,
            cbar=False,
            annot_kws={"size": 8},
            linewidths=0,
            **kwargs,
        )

        ax.set_title(title, fontsize=10, pad=8)
        ax.set_xticklabels(
            [_clean_label(t.get_text()) for t in ax.get_xticklabels()],
            fontsize=6,
            rotation=90,
            ha="right",
        )
        ax.set_yticklabels(
            [_clean_label(t.get_text()) for t in ax.get_yticklabels()],
            fontsize=6,
        )
        ax.grid(False)
        ax.tick_params(length=0)

    def _add_colorbar(fig, axes_list, cmap, pad=0.4):
        """
        Add a shared colorbar aligned to the actual heatmap extent.

        Reads the rendered positions of the top-right and bottom-right
        axes so the colorbar matches the heatmap height exactly.
        """
        fig.canvas.draw()

        # Find the bounding box of all heatmap axes
        top = max(ax.get_position().y1 for ax in axes_list)
        bottom = min(ax.get_position().y0 for ax in axes_list)
        right = max(ax.get_position().x1 for ax in axes_list)

        cbar_width = 0.012
        cbar_left = right + pad * cbar_width + 0.01
        cbar_ax = fig.add_axes([cbar_left, bottom, cbar_width, top - bottom])

        sm = ScalarMappable(cmap=cmap, norm=Normalize(vmin=-1, vmax=1))
        sm.set_array([])
        fig.colorbar(sm, cax=cbar_ax)
        cbar_ax.tick_params(labelsize=7, length=2)
        cbar_ax.grid(False)

    # Single panel
    if library_key is None:
        df = adata.obs[score_keys].dropna()
        corr, pvals = compute_corr_and_pval(df, method=method)

        fig_w = panel_w + 0.6  # extra for colorbar
        fig_h = panel_h
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        _plot_heatmap(corr, pvals, "Correlation heatmap", ax)

        fig.tight_layout()
        _add_colorbar(fig, [ax], cmap)

    # Multi-panel by batch
    else:
        batches = sorted(adata.obs[library_key].dropna().unique())
        n_batches = len(batches)
        ncols = int(np.ceil(np.sqrt(n_batches)))
        nrows = int(np.ceil(n_batches / ncols))

        gap_w = 0.15  # inches between columns
        gap_h = 0.35  # inches between rows
        fig_w = panel_w * ncols + gap_w * (ncols - 1) + 0.6  # +cbar
        fig_h = panel_h * nrows + gap_h * (nrows - 1)

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(fig_w, fig_h),
            squeeze=False,
        )

        active_axes = []
        for i, batch in enumerate(batches):
            row, col = divmod(i, ncols)
            ax = axes[row, col]

            df = adata[adata.obs[library_key] == batch].obs[score_keys].dropna()
            corr, pvals = compute_corr_and_pval(df, method=method)
            _plot_heatmap(corr, pvals, batch, ax)
            active_axes.append(ax)

        # Blank unused cells
        for j in range(n_batches, nrows * ncols):
            row, col = divmod(j, ncols)
            axes[row, col].axis("off")

        fig.tight_layout()
        _add_colorbar(fig, active_axes, cmap)

    # Save
    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")

    plt.show()
