from __future__ import annotations

from anndata import AnnData
import squidpy as sq
import scanpy as sc
import numpy as np

from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from scipy.sparse import issparse
from statsmodels.gam.api import GLMGam, BSplines

from ._infer_gene_weights import infer_gene_weights

sc.settings.verbosity = 0

def score(
    adata: AnnData,
    gene_set: list | dict,
    gene_weights: dict | None = None,
    score_key: str | list | None = None,
    spatial_key: str | None = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: str | None = None
) -> AnnData | None:
    """
    Compute spatially smoothed and spatially corrected gene set enrichment scores for one or more gene signatures.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, containing expression values and spatial coordinates in `obsm`.

    gene_set : list or dict
        Gene set(s) to be scored. If a list is provided, it is interpreted as a single gene signature.
        If a dict is provided, keys are signature names and values are lists of gene symbols.

    gene_weights : dict, optional
        Dictionary mapping signature names to dictionaries of gene weights (default is None). 
        If None, gene weights are inferred automatically.

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

    Returns
    -------
    AnnData
        The original `AnnData` object with additional scores stored in `adata.obs` under the key `{signature_name}_score`
        and gene contributions stored in `adata.uns["gene_contributions"]`.
    """

    if isinstance(gene_set, list):
        gene_set = {score_key or "enrichmap": gene_set}
    
    inferred_gene_weights = {}
    gene_weights = gene_weights or {}

    # Initialise or extend gene_contributions dict
    if "gene_contributions" not in adata.uns:
        adata.uns["gene_contributions"] = {}

    for sig_name, genes in tqdm(gene_set.items(), desc="Scoring signatures"):
        common_genes = list(set(genes).intersection(set(adata.var_names)))
        if len(common_genes) == 0:
            raise ValueError(f"No common genes found for signature {sig_name}")
        if len(common_genes) < 2:
            raise ValueError(f"Signature '{sig_name}' has fewer than two genes present in the dataset. At least two genes are required.")

        if sig_name not in gene_weights:
            inferred_gene_weights[sig_name] = infer_gene_weights(adata, common_genes)
        
        current_gene_weights = inferred_gene_weights.get(sig_name, gene_weights.get(sig_name, {}))
        if len(current_gene_weights) != len(common_genes):
            raise ValueError(f"The number of gene weights does not match the number of genes in the gene set for {sig_name}.")
        
        gene_weight_dict = dict(zip(common_genes, [current_gene_weights.get(gene, 1) for gene in common_genes]))
        all_z_scores = np.zeros(adata.n_obs)
        contribution_matrix = {}

        for idx, gene in enumerate(common_genes):
            gene_expr = adata[:, gene].X.toarray() if issparse(adata[:, gene].X) else adata[:, gene].X
            weighted_expr = gene_expr.flatten() * gene_weight_dict.get(gene, 1)
            all_z_scores += weighted_expr
            contribution_matrix[gene] = weighted_expr

        all_z_scores /= np.sum(list(gene_weight_dict.values()))
        z_scores = all_z_scores

        if batch_key is not None and batch_key in adata.obs.columns:
            batch_labels = adata.obs[batch_key].values
            scaler = StandardScaler()
            for batch in np.unique(batch_labels):
                mask = batch_labels == batch
                z_scores[mask] = scaler.fit_transform(z_scores[mask].reshape(-1, 1)).flatten()

        if smoothing:
            smoothed_scores = np.zeros_like(z_scores)
            batch_values = adata.obs[batch_key].unique() if batch_key else [None]
            for batch in batch_values:
                mask = adata.obs[batch_key] == batch if batch_key else np.ones(adata.n_obs, dtype=bool)
                adata_batch = adata[mask].copy()
                sq.gr.spatial_neighbors(adata_batch, n_neighs=n_neighbors, coord_type="generic", key_added="spatial")
                conn = adata_batch.obsp["spatial_connectivities"]
                smoothed = conn.dot(z_scores[mask]) / conn.sum(axis=1).A1
                smoothed_scores[mask] = smoothed
        else:
            smoothed_scores = z_scores

        if correct_spatial_covariates:
            coords = adata.obsm[spatial_key]
            bs = BSplines(coords, df=[10, 10], degree=[3, 3])
            gam = GLMGam.from_formula("smoothed_scores ~ 1", data=adata.obs.assign(smoothed_scores=smoothed_scores), smoother=bs)
            result = gam.fit()
            corrected_scores = smoothed_scores - result.fittedvalues
        else:
            corrected_scores = smoothed_scores

        adata.obs[f"{sig_name}_score"] = corrected_scores

        # Append to gene_contributions without overwriting previous signatures
        adata.uns["gene_contributions"][sig_name] = contribution_matrix

    return None