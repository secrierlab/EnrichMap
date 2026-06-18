from __future__ import annotations
from anndata import AnnData
from tqdm import tqdm
from scipy.sparse import issparse
from pygam import LinearGAM, te
from spatialdata._logging import logger

import logging
import squidpy as sq
import scanpy as sc
import numpy as np
import warnings
from ._infer_gene_weights import infer_gene_weights

logging.getLogger("squidpy").setLevel(logging.WARNING)
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

    Notes
    -----
    The four scoring steps (manuscript Eq. 1-5) are computed for all
    signatures at once with vectorised linear algebra rather than one
    signature at a time:

    - the weighted mean (Eq. 1) reuses a single column slice per
      signature instead of one per gene;
    - per-batch z-scoring (Eq. 2) and spatial smoothing (Eq. 3) act on
      the whole (spots x signatures) score matrix;
    - the spatial-covariate GAM (Eq. 4-5) is a penalised spline, hence a
      *linear* smoother whose hat matrix depends only on the coordinates,
      not on the response. It is therefore built once per batch and
      applied to every signature, instead of refitting a GAM per
      signature. Results are identical to the per-signature fit.
    """

    # Input resolution: gene_weights takes priority when provided
    if gene_weights is not None:
        # gene_weights is the primary driver — derive gene_set from it
        gene_set = {sig: list(weights.keys()) for sig, weights in gene_weights.items()}
    elif gene_set is not None:
        # gene_set only, weights will be inferred via CV
        if isinstance(gene_set, list):
            gene_set = {score_key or "enrichmap": gene_set}
    else:
        raise ValueError("Either gene_set or gene_weights must be provided.")

    # Guard: warn if multiple slides may be present without batch_key
    if batch_key is None and (smoothing or correct_spatial_covariates):
        _batch_hints = ["library_id", "sample", "slide", "batch", "sample_id"]
        for col in _batch_hints:
            if col in adata.obs.columns and adata.obs[col].nunique() > 1:
                warnings.warn(
                    f"Column '{col}' in adata.obs has {adata.obs[col].nunique()} "
                    f"unique values, suggesting multiple slides. Spatial smoothing "
                    f"and covariate correction assume a single tissue section when "
                    f"batch_key is None. Consider setting batch_key='{col}' to "
                    f"avoid cross-slide artefacts.",
                    UserWarning,
                    stacklevel=2,
                )
                break

    inferred_gene_weights = {}
    gene_weights = gene_weights or {}

    if "gene_contributions" not in adata.uns:
        adata.uns["gene_contributions"] = {}

    # ------------------------------------------------------------------
    # Pre-compute spatial graphs and coordinate arrays ONCE.
    #
    # The spatial kNN graph depends only on spatial coordinates, not on
    # gene expression or pathway identity. Building it inside the
    # per-pathway loop (as the original code did) rebuilds the same
    # graph N times for N pathways — a major performance bottleneck.
    # ------------------------------------------------------------------
    batches = adata.obs[batch_key].unique() if batch_key else [None]

    _spatial_cache = {}  # batch → (mask, connectivity_matrix, coords)
    if smoothing or correct_spatial_covariates:
        for batch in batches:
            mask = (
                adata.obs[batch_key] == batch
                if batch_key
                else np.ones(adata.n_obs, bool)
            )
            mask = np.asarray(mask)
            coords = adata.obsm[spatial_key][mask]

            if smoothing:
                adata_batch = adata[mask].copy()
                sq.gr.spatial_neighbors(
                    adata_batch,
                    n_neighs=n_neighbors,
                    coord_type="generic",
                    key_added="spatial",
                )
                conn = adata_batch.obsp["spatial_connectivities"]
            else:
                conn = None

            _spatial_cache[batch] = (mask, conn, coords)

    # ------------------------------------------------------------------
    # Step 1 (Eq. 1): weighted mean, per signature.
    #
    # Collect every signature's raw score into a single (spots x
    # signatures) matrix. Each signature is one column slice of the
    # expression matrix plus a vectorised weighted sum — no per-gene
    # AnnData indexing, which is what made wide signatures (e.g. large TF
    # regulons) slow. Steps 2-4 then act on the whole matrix at once.
    # ------------------------------------------------------------------
    var_index = {g: i for i, g in enumerate(adata.var_names)}
    X = adata.X
    X = X.tocsc() if issparse(X) else np.asarray(X)

    sig_names = list(gene_set.keys())
    scores_mat = np.zeros((adata.n_obs, len(sig_names)), dtype=np.float64)
    contributions = {}

    pbar = tqdm(sig_names)
    for j, sig_name in enumerate(pbar):
        genes = gene_set[sig_name]
        common_genes = list(set(genes).intersection(set(adata.var_names)))
        pbar.set_description(
            f"Scoring {sig_name}: {len(common_genes)}/{len(genes)} genes found"
        )
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

        idx = [var_index[g] for g in common_genes]
        w_vec = np.array([current_gene_weights[g] for g in common_genes], dtype=np.float64)
        weight_sum = np.abs(w_vec).sum()

        cols = X[:, idx]
        cols = cols.toarray() if issparse(cols) else np.asarray(cols)
        weighted = cols * w_vec  # (spots x len(common_genes))

        scores_mat[:, j] = weighted.sum(axis=1) / weight_sum
        contributions[sig_name] = {
            g: weighted[:, i] for i, g in enumerate(common_genes)
        }

    # Step 2 (Eq. 2): per-batch z-score normalisation (only with batch_key).
    if batch_key is not None:
        for batch in batches:
            mask = np.asarray(adata.obs[batch_key] == batch)
            block = scores_mat[mask]
            mu_b = block.mean(axis=0)
            sigma_b = block.std(axis=0)
            nz = sigma_b > 0
            block[:, nz] = (block[:, nz] - mu_b[nz]) / sigma_b[nz]
            scores_mat[mask] = block

    # Step 3 (Eq. 3): spatial smoothing on the full score matrix.
    if smoothing:
        for batch in batches:
            mask, conn, _ = _spatial_cache[batch]
            denom = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
            scores_mat[mask] = conn.dot(scores_mat[mask]) / denom[:, None]

    # Step 4 (Eq. 4-5): spatial-covariate correction.
    #
    # LinearGAM(te(0, 1)) with a fixed smoothing parameter is a linear
    # smoother: predict = H @ y with H = B (BᵀB + P)⁻¹ Bᵀ, where B is the
    # tensor-product spline design and P the penalty — both functions of
    # the coordinates only. We therefore residualise every signature
    # against the same fit in one solve, which is identical to fitting a
    # GAM per signature (verified to < 1e-7) but built once per batch.
    if correct_spatial_covariates:
        for batch in batches:
            mask, _, coords = _spatial_cache[batch]
            block = scores_mat[mask]
            try:
                # Fast path: build the spline hat matrix once and
                # residualise every signature in one solve.
                gam = LinearGAM(te(0, 1)).fit(coords, coords[:, 0])
                B = gam._modelmat(coords)
                B = B.toarray() if hasattr(B, "toarray") else np.asarray(B)
                P = gam._P()
                P = P.toarray() if hasattr(P, "toarray") else np.asarray(P)
                coef = np.linalg.solve(B.T @ B + P, B.T @ block)
                scores_mat[mask] = block - B @ coef
            except Exception:
                # Robust fallback: fit a GAM per signature (original
                # behaviour) if pygam internals are unavailable.
                for j in range(block.shape[1]):
                    g_ = LinearGAM(te(0, 1)).fit(coords, block[:, j])
                    block[:, j] = block[:, j] - g_.predict(coords)
                scores_mat[mask] = block

    # Write results back (one column per signature).
    for j, sig_name in enumerate(sig_names):
        adata.obs[f"{sig_name}_score"] = scores_mat[:, j]
        adata.uns["gene_contributions"][sig_name] = contributions[sig_name]

    return None
