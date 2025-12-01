from __future__ import annotations

from anndata import AnnData
import squidpy as sq
import scanpy as sc
import numpy as np

from tqdm import tqdm
from scipy.sparse import issparse
from statsmodels.gam.api import GLMGam, BSplines
from scipy.stats import median_abs_deviation

from ._infer_gene_weights import infer_gene_weights

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
    clip_min: float = -5.0,
    clip_max: float = 5.0,
) -> AnnData | None:
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
        If the coordinates are stored in a different key, specify that key here.

    n_neighbours : int, default 6
        Number of nearest spatial neighbours used for smoothing. This is passed to `squidpy.gr.spatial_neighbors`.

    smoothing : bool, default True
        Whether to perform spatial smoothing of signature scores using neighbour connectivity.

    correct_spatial_covariates : bool, default True
        Whether to correct scores for spatial covariates using a Generalised Additive Model.

    batch_key : str or None, optional
        Column in `adata.obs` indicating batch labels for batch-wise z-score normalisation and smoothing.
        If None, all cells are treated as a single batch.

    clip_min : float, default -5.0
        Minimum value to clip the final scores to limit the influence of outliers.
    clip_max : float, default 5.0
        Maximum value to clip the final scores to limit the influence of outliers.

    Returns
    -------
    AnnData
        The original `AnnData` object with additional scores stored in `adata.obs` under the key `{signature_name}_score`
        and gene contributions stored in `adata.uns["gene_contributions"]`.
    """

    # Determine gene_set if not provided
    if gene_set is None:
        if gene_weights is not None:
            # Use genes from gene_weights
            gene_set = {
                sig: list(weights.keys()) for sig, weights in gene_weights.items()
            }
        else:
            raise ValueError("Either gene_set or gene_weights must be provided.")

    # Convert list → dict if needed
    if isinstance(gene_set, list):
        gene_set = {score_key or "enrichmap": gene_set}

    inferred_gene_weights = {}
    gene_weights = gene_weights or {}

    if "gene_contributions" not in adata.uns:
        adata.uns["gene_contributions"] = {}

    # Iterate over signatures
    for sig_name, genes in tqdm(gene_set.items(), desc="Scoring signatures"):
        # Restrict to genes present in adata
        common_genes = list(set(genes).intersection(set(adata.var_names)))
        if len(common_genes) == 0:
            raise ValueError(f"No common genes found for signature {sig_name}")
        if len(common_genes) < 2:
            raise ValueError(
                f"Signature '{sig_name}' has fewer than two genes in the dataset"
            )

        # Determine gene weights
        if sig_name in gene_weights and gene_weights[sig_name] is not None:
            # Use provided weights, restrict to common genes
            current_gene_weights = {
                g: gene_weights[sig_name][g]
                for g in common_genes
                if g in gene_weights[sig_name]
            }
            # Fill missing genes with weight 1.0
            for g in common_genes:
                if g not in current_gene_weights:
                    current_gene_weights[g] = 1.0
        else:
            inferred_gene_weights[sig_name] = infer_gene_weights(adata, common_genes)
            current_gene_weights = inferred_gene_weights[sig_name]

        # Align dictionary of weights to common genes
        gene_weight_dict = {g: current_gene_weights.get(g, 1.0) for g in common_genes}

        # Compute weighted expression
        all_weighted = np.zeros(adata.n_obs)
        contribution_matrix = {}

        for gene in common_genes:
            expr = adata[:, gene].X
            expr = expr.toarray().flatten() if issparse(expr) else expr.flatten()
            weighted_expr = expr * gene_weight_dict[gene]
            all_weighted += weighted_expr
            contribution_matrix[gene] = weighted_expr

        # Normalise by total weight
        all_weighted /= np.sum(list(gene_weight_dict.values()))
        raw_scores = all_weighted.copy()

        # Spatial smoothing
        if smoothing:
            smoothed_scores = np.zeros_like(raw_scores)
            batch_values = adata.obs[batch_key].unique() if batch_key else [None]

            for batch in batch_values:
                mask = (
                    adata.obs[batch_key] == batch
                    if batch_key
                    else np.ones(adata.n_obs, dtype=bool)
                )
                adata_batch = adata[mask].copy()

                sq.gr.spatial_neighbors(
                    adata_batch,
                    n_neighs=n_neighbors,
                    coord_type="generic",
                    key_added="spatial",
                )

                conn = adata_batch.obsp["spatial_connectivities"]
                smoothed = conn.dot(raw_scores[mask]) / np.maximum(
                    conn.sum(axis=1).A1, 1e-10
                )
                smoothed_scores[mask] = smoothed

        else:
            smoothed_scores = raw_scores

        # Spatial covariate correction
        if correct_spatial_covariates:
            corrected_scores = np.zeros_like(smoothed_scores)
            batch_values = adata.obs[batch_key].unique() if batch_key else [None]

            for batch in batch_values:
                mask = (
                    adata.obs[batch_key] == batch
                    if batch_key
                    else np.ones(adata.n_obs, dtype=bool)
                )
                coords = adata.obsm[spatial_key][mask]

                bs = BSplines(coords, df=[10, 10], degree=[3, 3])
                gam = GLMGam.from_formula(
                    "y ~ 1",
                    data={"y": smoothed_scores[mask]},
                    smoother=bs,
                )
                result = gam.fit()

                corrected_scores[mask] = smoothed_scores[mask] - result.fittedvalues

        else:
            corrected_scores = smoothed_scores

        # Global robust scaling
        median = np.median(corrected_scores)
        mad = median_abs_deviation(corrected_scores, scale="normal")
        corrected_scores = (corrected_scores - median) / mad

        # Clipping
        corrected_scores = np.clip(corrected_scores, clip_min, clip_max)

        # Store outputs
        adata.obs[f"{sig_name}_score"] = corrected_scores
        adata.uns["gene_contributions"][sig_name] = contribution_matrix

    return None
