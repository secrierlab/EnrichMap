from __future__ import annotations

from sklearn.metrics import f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import Sequence
from anndata import AnnData

def compute_f1_scores(
    adata: AnnData,
    score_keys: Sequence[str],
    binary_labels: np.ndarray,
    figsize: tuple[int, int] = (3, 3),
    save: str | None = None
) -> dict[str, float]:
    """
    Compute and plot F1 scores for multiple gene set scoring methods.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with scoring results in `adata.obs`.
    score_keys : list of str
        List of method score keys to benchmark.
    binary_labels : np.ndarray
        Ground truth binary labels.
    figsize : tuple, optional
        Size of the figure.
    save : str or None, optional
        If given, path to save the figure.

    Returns
    -------
    f1_scores : dict
        Dictionary of F1 scores per `score_keys` provided.
    """
    f1_scores = {}

    for method in score_keys:
        if method in adata.obs.columns:
            preds = adata.obs[method] > np.median(adata.obs[method])
            f1_scores[method] = f1_score(binary_labels, preds)

    # Plot
    plt.figure(figsize=figsize)
    sns.barplot(x=list(f1_scores.keys()), y=list(f1_scores.values()), palette="muted", edgecolor=None)
    plt.title("F1 scores", fontsize=10)
    plt.ylabel("F1 score", fontsize=8)
    plt.xlabel("Method", fontsize=8)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=8)
    plt.grid(False)

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()

    return f1_scores
