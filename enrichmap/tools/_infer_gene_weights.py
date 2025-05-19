from __future__ import annotations

import numpy as np

from anndata import AnnData
from scipy.sparse import issparse

def infer_gene_weights(
    adata: AnnData,
    gene_set: list,
) -> dict | None:
    """
    Infer gene weights based on the coefficient of variation (CV) of gene expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, containing gene expression values.

    gene_set : list
        List of gene names for which weights (CV values) will be inferred.
        Only genes present in `adata.var_names` will be considered.

    Returns
    -------
    dict
        A dictionary mapping gene names to their inferred weights based on the coefficient of variation.
        The weight reflects the relative variability of each gene's expression across all cells.

    """
    # Ensure genes are in the dataset
    common_genes = list(set(gene_set).intersection(set(adata.var_names)))
    if len(common_genes) == 0:
        raise ValueError("No common genes found between gene set and dataset")
    
    # Subset the data to only include the common genes
    expr_matrix = adata[:, common_genes].X
    if issparse(expr_matrix):
        expr_matrix = expr_matrix.toarray()  # Convert sparse matrix to dense if needed
    
    # Compute the mean and standard deviation for each gene
    mean_expr = np.mean(expr_matrix, axis=0)
    std_expr = np.std(expr_matrix, axis=0)
    
    # Avoid division by zero and calculate CV (coefficient of variation)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.divide(std_expr, mean_expr)
        cv[mean_expr == 0] = 0  # Set CV to 0 where the mean expression is zero

    # Create a dictionary mapping gene names to their CV (gene weights)
    gene_weights = dict(zip(common_genes, cv))
    
    return gene_weights