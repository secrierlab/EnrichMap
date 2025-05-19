from __future__ import annotations

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from anndata import AnnData

plt.rcParams["pdf.fonttype"] = "truetype"

def signature_correlation_heatmap(
    adata: AnnData,
    signatures: list[str],
    batch_key: str | None = None,
    method: str = "spearman",
    save: str | Path | None = None
):
    """
    Plot correlation heatmaps between gene signature scores, optionally per batch.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with signature scores in `adata.obs`.
    signatures : list of str
        List of gene set score keys to correlate.
    batch_key : str or None
        Column in `adata.obs` indicating batches/libraries.
        If None, one global heatmap is plotted.
    method : str
        Correlation method ('spearman', 'pearson', etc.).
    save : str, Path or None
        If provided, base filename to save figure(s).
    """
    if batch_key is None:
        df = adata.obs[signatures].dropna()
        corr = df.corr(method=method)
        plt.figure(figsize=(6, 6), constrained_layout=True)
        sns.heatmap(corr, annot=True, cmap="coolwarm", center=0, vmin=-1, vmax=1, cbar=False, annot_kws={"size": 8})
        plt.title("Correlation of gene set scores")
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.grid(False)
        if save:
            # Ensure 'figures/' directory exists
            os.makedirs("figures", exist_ok=True)

            # If 'save' has no directory path, prepend 'figures/'
            if not os.path.dirname(save):
                save = os.path.join("figures", save)
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
    else:
        batch_values = adata.obs[batch_key].unique()
        n = len(batch_values)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows), constrained_layout=True)
        axes = np.ravel(axes)

        for i, batch in enumerate(batch_values):
            df = adata[adata.obs[batch_key] == batch].obs[signatures].dropna()
            corr = df.corr(method=method)
            sns.heatmap(corr, annot=True, cmap="coolwarm", center=0, vmin=-1, vmax=1, ax=axes[i], cbar=False, annot_kws={"size": 8})
            axes[i].set_title(f"{batch_key}: {batch}")
            axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=8)
            axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=8)
            axes[i].grid(False)

        # Hide unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")

        if save:
            # Ensure 'figures/' directory exists
            os.makedirs("figures", exist_ok=True)

            # If 'save' has no directory path, prepend 'figures/'
            if not os.path.dirname(save):
                save = os.path.join("figures", save)
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()