from __future__ import annotations

import numpy as np

from anndata import AnnData
from scipy.sparse import issparse
from scipy.stats import spearmanr

from ._infer_gene_weights import infer_gene_weights


def build_signed_weights(
    adata: AnnData,
    up_genes: list[str],
    down_genes: list[str],
    score_key: str = "enrichmap",
) -> dict[str, dict[str, float]]:
    """
    Build a signed gene weight dictionary from upregulated and downregulated gene lists.

    CV-inferred weights capture spatial variability (magnitude). The sign
    is determined by differential correlation: for each gene, the Spearman
    correlation with the mean expression of upregulated genes minus the
    correlation with the mean expression of downregulated genes. This
    difference cancels shared technical inflation (e.g. library size) and
    produces genuinely signed weights. Each group is normalised by its own
    size so that upregulated and downregulated genes contribute equally
    regardless of gene count.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    up_genes : list of str
        Upregulated gene symbols.
    down_genes : list of str
        Downregulated gene symbols.
    score_key : str, default "enrichmap"
        Name for the signature in the returned dictionary.

    Returns
    -------
    dict
        Gene weights dictionary ready for ``enrichmap.tl.score(gene_weights=...)``.

    Examples
    --------
    >>> from enrichmap.tools import build_signed_weights
    >>> gw = build_signed_weights(adata, up_genes=up, down_genes=down, score_key="G0")
    >>> em.tl.score(adata, gene_weights=gw)
    """
    all_genes = list(set(up_genes + down_genes))
    cv_weights = infer_gene_weights(adata, all_genes)

    # Filter to genes present in data
    up_present = [g for g in up_genes if g in cv_weights]
    down_present = [g for g in down_genes if g in cv_weights]

    n_up = len(up_present)
    n_down = len(down_present)

    if n_up == 0:
        raise ValueError("No upregulated genes found in the dataset.")
    if n_down == 0:
        raise ValueError("No downregulated genes found in the dataset.")

    # Compute both reference patterns
    up_expr = adata[:, up_present].X
    if issparse(up_expr):
        up_expr = up_expr.toarray()
    up_reference = np.mean(up_expr, axis=1)

    down_expr = adata[:, down_present].X
    if issparse(down_expr):
        down_expr = down_expr.toarray()
    down_reference = np.mean(down_expr, axis=1)

    # Differential correlation determines sign
    all_present = up_present + down_present
    weights = {}

    for g in all_present:
        expr = adata[:, g].X
        expr = expr.toarray().flatten() if issparse(expr) else expr.flatten()

        r_up, _ = spearmanr(expr, up_reference)
        r_down, _ = spearmanr(expr, down_reference)
        r_up = r_up if not np.isnan(r_up) else 0
        r_down = r_down if not np.isnan(r_down) else 0

        # Sign from differential correlation, magnitude from CV
        delta_r = r_up - r_down

        if g in up_present:
            weights[g] = cv_weights[g] * delta_r / n_up
        else:
            weights[g] = cv_weights[g] * delta_r / n_down

    return {score_key: weights}
