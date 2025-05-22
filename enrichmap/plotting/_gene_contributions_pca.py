from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt

from anndata import AnnData
from sklearn.decomposition import PCA
from adjustText import adjust_text

plt.rcParams["pdf.fonttype"] = "truetype"

def gene_contributions_pca(
    adata: AnnData,
    score_key: str,
    top_n_genes: int = 5,
    highlight_genes: list[str] | None = None,
    figsize: tuple = (4, 4),
    fontsize: int = 8,
    save: str | None = None
) -> None:
    """
    Visualise gene contributions to a spatial gene signature using PCA-reduced scatter plot.

    Each point represents a gene's contribution vector projected into 2D via PCA. 
    The top contributing genes and any specified genes of interest can be highlighted.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must contain `gene_contributions` for the given `signature_name` in `adata.uns`.

    score_key : str
        Name of the gene signature whose gene contributions are to be visualised.
        Should match a key in `adata.uns["gene_contributions"]`.

    top_n_genes : int, default 5
        Number of top genes (by vector norm) to highlight in the plot.

    highlight_genes : list of str, optional
        Specific gene names to highlight, regardless of their contribution magnitude.

    figsize : tuple, default (4, 4)
        Size of the resulting plot.

    fontsize : int, default 8
        Font size used for gene labels in the plot.

    save : str, optional
        Path to save the plot. If None, the plot is only shown and not saved.

    Returns
    -------
    None
        Displays the PCA plot and optionally saves it to a file.

    Raises
    ------
    ValueError
        If `gene_contributions` for the specified signature is not found in `adata.uns`.
    """
    if "gene_contributions" not in adata.uns:
        raise ValueError("Gene contributions not found in adata.uns.")

    if score_key.endswith("_score"):
        score_key = score_key[:-6]

    contrib_dict = adata.uns["gene_contributions"][score_key]

    # Stack into matrix: genes x spatial spots
    gene_names = list(contrib_dict.keys())
    X = np.vstack([contrib_dict[gene] for gene in gene_names])
    
    # PCA to reduce to 2D
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(X)

    # Select top N genes by norm of their contribution vector
    norms = np.linalg.norm(X, axis=1)
    top_idx = np.argsort(norms)[-top_n_genes:]

    # Prepare for highlighting
    highlight_genes = highlight_genes or []
    highlight_set = set(highlight_genes)

    plt.figure(figsize=figsize)
    plt.scatter(coords[:, 0], coords[:, 1], alpha=0.2, s=10, color="grey")

    texts = []

    for idx in top_idx:
        gene = gene_names[idx]
        if gene in highlight_set:
            continue
        x, y = coords[idx]
        plt.scatter(x, y, s=20, color="lightblue")
        texts.append(plt.text(x, y, gene, fontsize=fontsize, ha="center", va="center"))

    for gene in highlight_genes:
        if gene in gene_names:
            idx = gene_names.index(gene)
            x, y = coords[idx]
            plt.scatter(x, y, s=20, color="grey")
            texts.append(plt.text(x, y, gene, fontsize=fontsize, ha="center", va="center"))

    # Add dashed lines at x = 0 and y = 0
    plt.axhline(0, color="grey", linestyle="--", linewidth=0.5)
    plt.axvline(0, color="grey", linestyle="--", linewidth=0.5)

    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(f"Gene contributions for {score_key}")
    plt.grid(False)
    if save:
        # Ensure 'figures/' directory exists
        os.makedirs("figures", exist_ok=True)

        # If 'save' has no directory path, prepend 'figures/'
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()