from __future__ import annotations

from anndata import AnnData
import squidpy as sq
import scanpy as sc
import numpy as np

from tqdm import tqdm
from scipy.sparse import issparse
from statsmodels.gam.api import GLMGam, BSplines
from scipy.stats import median_abs_deviation

from enrichmap.tools._infer_gene_weights import infer_gene_weights

sc.settings.verbosity = 0


def score(
    adata: AnnData,
    gene_set: list | dict | None = None,
    gene_weights: dict | None = None,
    score_key: str | list | None = None,
    spatial_key: str | None = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: str | None = None,
) -> None:
    """
    Compute spatially smoothed and spatially corrected gene set enrichment scores for one or more gene signatures.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, containing expression values and spatial coordinates in `obsm`.

    gene_set : list or dict or None
        Gene set(s) to be scored. If a list is provided, it is interpreted as a single gene signature.
        If a dict is provided, keys are signature names and values are lists of gene symbols.
        If None, `gene_weights` must be provided and gene sets will be inferred from the keys of `gene_weights`.

    gene_weights : dict, optional
        Dictionary mapping signature names to dictionaries of gene weights (default is None).
        If None, gene weights are inferred automatically. If provided, `gene_set` is overridden to match the keys.

    score_key : str, list, or None, optional
        Name or list of names to assign to the gene signature(s) if `gene_set` is provided as a list.
        Ignored if `gene_set` is already a dictionary.

    spatial_key : str
        Key in `adata.obsm` containing spatial coordinates used for spatial covariate correction. By default, it is set to "spatial".

    n_neighbours : int, default 6
        Number of nearest spatial neighbours used for smoothing.

    smoothing : bool, default True
        Whether to perform spatial smoothing of signature scores.

    correct_spatial_covariates : bool, default True
        Whether to correct scores for spatial covariates using a GAM.

    batch_key : str or None, optional
        Column in `adata.obs` indicating batch labels.

    Returns
    -------
    None
        Scores are stored in `adata.obs` and gene contributions in `adata.uns["gene_contributions"]`.
    """

    if gene_set is None:
        if gene_weights is not None:
            gene_set = {
                sig: list(weights.keys()) for sig, weights in gene_weights.items()
            }
        else:
            raise ValueError("Either gene_set or gene_weights must be provided.")

    if isinstance(gene_set, list):
        gene_set = {score_key or "enrichmap": gene_set}

    inferred_gene_weights = {}
    gene_weights = gene_weights or {}

    if "gene_contributions" not in adata.uns:
        adata.uns["gene_contributions"] = {}

    for sig_name, genes in tqdm(gene_set.items(), desc="Scoring signatures"):
        common_genes = list(set(genes).intersection(set(adata.var_names)))
        if len(common_genes) < 2:
            raise ValueError(
                f"Signature '{sig_name}' has fewer than two genes in the dataset"
            )

        # Determine gene weights
        if sig_name in gene_weights and gene_weights[sig_name] is not None:
            current_gene_weights = {
                g: gene_weights[sig_name].get(g, 1.0) for g in common_genes
            }
        else:
            inferred_gene_weights[sig_name] = infer_gene_weights(adata, common_genes)
            current_gene_weights = {
                g: inferred_gene_weights[sig_name].get(g, 1.0) for g in common_genes
            }

        # Compute weighted expression
        weighted_matrix = np.zeros(adata.n_obs)
        contribution_matrix = {}

        for gene in common_genes:
            expr = adata[:, gene].X
            expr = expr.toarray().flatten() if issparse(expr) else expr.flatten()
            weighted_expr = expr * current_gene_weights[gene]
            weighted_matrix += weighted_expr
            contribution_matrix[gene] = weighted_expr

        raw_scores = weighted_matrix.copy()

        # Spatial smoothing
        smoothed_scores = raw_scores.copy()
        if smoothing:
            batches = adata.obs[batch_key].unique() if batch_key else [None]
            for batch in batches:
                mask = (
                    adata.obs[batch_key] == batch
                    if batch_key
                    else np.ones(adata.n_obs, bool)
                )
                adata_batch = adata[mask].copy()
                sq.gr.spatial_neighbors(
                    adata_batch,
                    n_neighs=n_neighbors,
                    coord_type="generic",
                    key_added="spatial",
                )
                conn = adata_batch.obsp["spatial_connectivities"]
                smoothed_scores[mask] = conn.dot(raw_scores[mask]) / np.maximum(
                    conn.sum(axis=1).A1, 1e-10
                )

        # Spatial covariate correction
        corrected_scores = smoothed_scores.copy()
        if correct_spatial_covariates:
            batches = adata.obs[batch_key].unique() if batch_key else [None]
            for batch in batches:
                mask = (
                    adata.obs[batch_key] == batch
                    if batch_key
                    else np.ones(adata.n_obs, bool)
                )
                coords = adata.obsm[spatial_key][mask]
                bs = BSplines(
                    coords, df=[10] * coords.shape[1], degree=[3] * coords.shape[1]
                )
                gam = GLMGam.from_formula(
                    "y ~ 1", data={"y": smoothed_scores[mask]}, smoother=bs
                )
                result = gam.fit()
                corrected_scores[mask] = smoothed_scores[mask] - result.fittedvalues

        # Robust scaling per batch
        scaled_scores = np.zeros_like(corrected_scores)
        batches = adata.obs[batch_key].unique() if batch_key else [None]
        for batch in batches:
            mask = (
                adata.obs[batch_key] == batch
                if batch_key
                else np.ones(adata.n_obs, bool)
            )
            batch_scores = corrected_scores[mask]
            median = np.median(batch_scores)
            mad = median_abs_deviation(batch_scores, scale="normal")
            scaled_scores[mask] = (batch_scores - median) / mad

        adata.obs[f"{sig_name}_score"] = scaled_scores
        adata.uns["gene_contributions"][sig_name] = contribution_matrix

    return None
