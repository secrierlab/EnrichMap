from __future__ import annotations

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as gridspec

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib import pyplot as plt
from anndata import AnnData
from scipy.stats import spearmanr, pearsonr
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
    annot = pvals.copy()
    annot[:] = ""
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
    figsize: tuple[int, int] = (5, 5),
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

    figsize : tuple of int
        Size of the figure (width, height).

    cmap : str or None
        Colormap for the heatmap. Defaults to "coolwarm".

    save : str or Path or None
        Path to save the figure.

    **kwargs : dict
        Additional keyword arguments passed to `sns.heatmap`.
    """
    # optional library filtering
    if library_key and library_id is not None:
        mask = (
            adata.obs[library_key] == library_id
            if isinstance(library_id, str)
            else adata.obs[library_key].isin(library_id)
        )
        adata = adata[mask]

    # internal helper
    def plot_single_heatmap(corr, pvals, title, fig, position):
        annot = get_star_annot(pvals)

        # Compute distance matrix for clustering (1 - abs(corr))
        dist = 1 - corr
        condensed_dist = squareform(dist, checks=False)  # convert to 1D condensed form

        # Compute linkage and order
        linkage_matrix = linkage(condensed_dist, method="average")
        order = leaves_list(linkage_matrix)

        # Reorder rows and columns
        corr = corr.iloc[order, order]
        pvals = pvals.iloc[order, order]
        annot = annot.iloc[order, order]

        # Create main axis and colorbar axis
        ax = fig.add_subplot(position)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)

        sns.heatmap(
            corr,
            annot=annot,
            fmt="",
            cmap=cmap,
            square=True,
            ax=ax,
            vmax=1,
            vmin=-1,
            cbar=True,
            cbar_ax=cax,
            annot_kws={"size": 8},
            **kwargs,
        )

        ax.set_title(title, fontsize=10, pad=8)

        xticks = [
            label.get_text().replace("_score", "") for label in ax.get_xticklabels()
        ]
        yticks = [
            label.get_text().replace("_score", "") for label in ax.get_yticklabels()
        ]

        ax.set_xticklabels(xticks, fontsize=6, rotation=90, ha="right")
        ax.set_yticklabels(yticks, fontsize=6)
        ax.grid(False)
        cax.tick_params(labelsize=6)

    # single panel
    if library_key is None:
        df = adata.obs[score_keys].dropna()
        corr, pvals = compute_corr_and_pval(df, method=method)

        fig = plt.figure(figsize=figsize, constrained_layout=True)
        gs = gridspec.GridSpec(1, 1, figure=fig)
        plot_single_heatmap(corr, pvals, "Correlation heatmap", fig, gs[0])

    # multi-panel by batch
    else:
        batches = sorted(adata.obs[library_key].dropna().unique())
        n_batches = len(batches)
        ncols = int(np.ceil(np.sqrt(n_batches)))
        nrows = int(np.ceil(n_batches / ncols))

        fig = plt.figure(
            figsize=(figsize[0] * ncols, figsize[1] * nrows), constrained_layout=True
        )
        gs_outer = gridspec.GridSpec(nrows, ncols, figure=fig, wspace=0.4, hspace=0.4)

        for i, batch in enumerate(batches):
            row, col = divmod(i, ncols)
            idx = gs_outer[row, col]

            df = adata[adata.obs[library_key] == batch].obs[score_keys].dropna()
            corr, pvals = compute_corr_and_pval(df, method=method)
            plot_single_heatmap(
                corr, pvals, f"Correlation heatmap for {batch}", fig, idx
            )

        # blank unused cells
        for j in range(i + 1, nrows * ncols):
            fig.add_subplot(gs_outer[j // ncols, j % ncols]).axis("off")

    # optional save
    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")

    plt.show()
