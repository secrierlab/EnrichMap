from __future__ import annotations

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from anndata import AnnData
from scipy.stats import spearmanr, pearsonr

plt.rcParams["pdf.fonttype"] = "truetype"


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
    return pd.DataFrame(corr, index=df.columns, columns=df.columns), pd.DataFrame(
        pvals, index=df.columns, columns=df.columns
    )


def get_star_annot(pvals):
    annot = pvals.copy()
    annot[:] = ""
    for i in range(pvals.shape[0]):
        for j in range(pvals.shape[1]):
            p = pvals.iloc[i, j]
            if p < 0.001:
                annot.iloc[i, j] = "***"
            elif p < 0.01:
                annot.iloc[i, j] = "**"
            elif p < 0.05:
                annot.iloc[i, j] = "*"
            else:
                annot.iloc[i, j] = ""
    return annot


def signature_correlation_heatmap(
    adata: AnnData,
    score_keys: list[str],
    batch_key: str | None = None,
    method: str = "spearman",
    figsize: tuple[int, int] = (5, 4),
    save: str | Path | None = None,
):
    """
    Plot a heatmap of correlations between gene set scores in `adata.obs`.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, where gene signature scores are stored in `adata.obs`.

    score_keys : list of str
        List of column names in `adata.obs` corresponding to gene set signature scores.

    batch_key : str or None, optional (default: None)
        Key in `adata.obs` specifying batch labels. If provided, separate heatmaps
        will be plotted for each batch group.

    method : str, optional (default: "spearman")
        Correlation method to use. Can be `"spearman"` or `"pearson"`.

    save : str or Path or None, optional (default: None)
        Path to save the figure.
    """

    def plot_heatmap(corr, pvals, title=None, ax=None, cbar=False, cbar_ax=None):
        annot = get_star_annot(pvals)
        hm = sns.heatmap(
            corr,
            annot=annot,
            fmt="",
            cmap="seismic",
            center=0,
            vmin=-1,
            vmax=1,
            ax=ax,
            cbar=cbar,
            cbar_ax=cbar_ax,
            annot_kws={"size": 10},
        )
        if title:
            ax.set_title(title, fontsize=10)
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha="right")
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)
        ax.grid(False)
        return hm

    star_legend = {
        "***": "p < 0.001",
        "**": "p < 0.01",
        "*": "p < 0.05",
        "n.s.": "n.s.",
    }

    handles = [
        plt.Line2D(
            [],
            [],
            linestyle="None",
            marker="",
            markersize=10,
            markerfacecolor="black",
            markeredgewidth=0,
            label=f"{star} ({label})",
        )
        for star, label in star_legend.items()
    ]

    if batch_key is None:
        df = adata.obs[score_keys].dropna()
        corr, pvals = compute_corr_and_pval(df, method=method)

        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 2, width_ratios=[20, 1])
        ax = fig.add_subplot(gs[0])
        cbar_ax = fig.add_subplot(gs[1])

        plot_heatmap(
            corr,
            pvals,
            title="Correlation of gene set scores",
            ax=ax,
            cbar=True,
            cbar_ax=cbar_ax,
        )
        fig.legend(handles=handles, loc="lower left", fontsize=8)

        if save:
            os.makedirs("figures", exist_ok=True)
            if not os.path.dirname(save):
                save = os.path.join("figures", save)
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()

    else:
        batch_values = adata.obs[batch_key].unique()
        n = len(batch_values)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))

        fig_width = 3 * ncols + 1.5  # add space for colorbar
        fig_height = 3 * nrows
        fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
        axes = np.ravel(axes)

        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])

        mesh = None
        for i, batch in enumerate(batch_values):
            df = adata[adata.obs[batch_key] == batch].obs[score_keys].dropna()
            corr, pvals = compute_corr_and_pval(df, method=method)
            hm = plot_heatmap(
                corr, pvals, title=f"{batch_key}: {batch}", ax=axes[i], cbar=False
            )
            if mesh is None:
                mesh = hm.get_children()[0]

        for j in range(i + 1, len(axes)):
            axes[j].axis("off")

        if mesh:
            fig.colorbar(mesh, cax=cbar_ax)

        fig.legend(handles=handles, loc="lower left", fontsize=8)

        if save:
            os.makedirs("figures", exist_ok=True)
            if not os.path.dirname(save):
                save = os.path.join("figures", save)
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
