from __future__ import annotations

from anndata import AnnData
import squidpy as sq
import scanpy as sc
import numpy as np
import logging
from tqdm import tqdm
from scipy.sparse import issparse
from pygam import LinearGAM, te

from enrichmap.tools._infer_gene_weights import infer_gene_weights
from spatialdata._logging import logger

logger.setLevel(logging.ERROR)
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

    n_neighbors : int, default 6
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

        # Step 1: Compute weighted average
        # Z_j = (1 / sum(w_i)) * sum(w_i * x_ij)
        weighted_matrix = np.zeros(adata.n_obs)
        contribution_matrix = {}
        weight_sum = sum(current_gene_weights[g] for g in common_genes)

        for gene in common_genes:
            expr = adata[:, gene].X
            expr = expr.toarray().flatten() if issparse(expr) else expr.flatten()
            weighted_expr = expr * current_gene_weights[gene]
            weighted_matrix += weighted_expr
            contribution_matrix[gene] = weighted_expr

        raw_scores = weighted_matrix / weight_sum

        # Step 2: Batch z-score normalisation
        # Z_j_tilde^(b) = (Z_j - mu_b) / sigma_b
        if batch_key is not None:
            for batch in adata.obs[batch_key].unique():
                mask = adata.obs[batch_key] == batch
                batch_scores = raw_scores[mask]
                mu_b = np.mean(batch_scores)
                sigma_b = np.std(batch_scores)
                if sigma_b > 0:
                    raw_scores[mask] = (batch_scores - mu_b) / sigma_b

        # Step 3: Spatial smoothing
        # Z_j_tilde = (1 / sum(A_jk)) * sum(A_jk * Z_k)
        scores = raw_scores.copy()
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
                scores[mask] = conn.dot(scores[mask]) / np.maximum(
                    conn.sum(axis=1).A1, 1e-10
                )

        # Step 4: Spatial covariate correction
        # Z_j_tilde = alpha + f(x_j, y_j) + epsilon_j
        # Corrected_j = Z_j_tilde - Z_j_tilde_hat
        if correct_spatial_covariates:
            batches = adata.obs[batch_key].unique() if batch_key else [None]
            for batch in batches:
                mask = (
                    adata.obs[batch_key] == batch
                    if batch_key
                    else np.ones(adata.n_obs, bool)
                )
                coords = adata.obsm[spatial_key][mask]
                gam = LinearGAM(te(0, 1, n_splines=[10, 10])).fit(coords, scores[mask])
                scores[mask] = scores[mask] - gam.predict(coords)

        adata.obs[f"{sig_name}_score"] = scores
        adata.uns["gene_contributions"][sig_name] = contribution_matrix

    return None
