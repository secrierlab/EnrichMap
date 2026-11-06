from __future__ import annotations
from anndata import AnnData
from tqdm import tqdm
from scipy.sparse import issparse, csc_matrix, csr_matrix
from scipy.stats import rankdata
from pygam import LinearGAM, te
from spatialdata._logging import logger

import logging
import zlib
import squidpy as sq
import scanpy as sc
import numpy as np
import warnings
from ._infer_gene_weights import infer_gene_weights

logging.getLogger("squidpy").setLevel(logging.WARNING)
logger.setLevel(logging.ERROR)


sc.settings.verbosity = 0

# Null model: each signature gene is replaced by one of the
# `_N_MATCHED_CANDIDATES` genes closest to it in CV among the `_N_MEAN_WINDOW`
# non-signature genes closest to it in mean expression (both widened when the
# pools are too small to give distinct control sets), and the per-spot null SD
# is moderated towards its median with `_SD_PRIOR_DF` prior degrees of freedom.
_N_MEAN_WINDOW = 100
_N_MATCHED_CANDIDATES = 30
_SD_PRIOR_DF = 4
# Control sets sharing more genes than this (mean Jaccard) widen their pools / are flagged; offsets larger than
# this many null SDs mean the signature genes have no expression peers.
_MAX_CONTROL_OVERLAP = 0.1
_MAX_OFFSET_RATIO = 5.0
# Signed control sets: spots x sets x genes handled at once when refitting weights (about 200 MB per float64
# array), and the most spots x genes ranks kept for reuse across the signatures of one call (about 800 MB).
_RANK_BLOCK_ELEMENTS = 50_000_000
_RANK_CACHE_ELEMENTS = 100_000_000


def _matched_candidates(
    mean_rank: np.ndarray,
    cv_rank: np.ndarray,
    idx: np.ndarray,
    free: np.ndarray,
    n_mean: int,
    n_cv: int,
    tiebreak: np.ndarray,
) -> np.ndarray:
    """
    Candidate control genes for each signature gene.

    Parameters
    ----------
    mean_rank, cv_rank : np.ndarray
        Rank of every gene's mean expression and coefficient of variation.
    idx : np.ndarray
        Indices of the signature genes.
    free : np.ndarray
        Indices of the genes that may be drawn as controls.
    n_mean, n_cv : int
        Size of the mean-expression window and number of candidates kept from it by CV.
    tiebreak : np.ndarray
        Rank of every gene's name, used to break ties and order rows.

    Returns
    -------
    np.ndarray
        (len(idx) x n_cv) indices, each row ordered by gene name: per signature gene, the free genes
        closest in CV among those closest in mean expression. Ties are broken by gene name, so the result
        depends neither on the numpy build nor on the column order of the data.
    """
    n_mean = min(n_mean, len(free))
    n_cv = min(n_cv, n_mean)
    order = free[np.lexsort((tiebreak[free], mean_rank[free]))]  # free genes by (mean rank, name)
    sorted_rank = mean_rank[order]
    out = np.empty((len(idx), n_cv), dtype=np.int64)
    for i, g in enumerate(idx):
        p = np.searchsorted(sorted_rank, mean_rank[g])
        near = order[max(0, p - n_mean): min(len(order), p + n_mean)]
        window = near[np.lexsort((tiebreak[near], np.abs(mean_rank[near] - mean_rank[g])))][:n_mean]
        chosen = window[np.lexsort((tiebreak[window], np.abs(cv_rank[window] - cv_rank[g])))][:n_cv]
        out[i] = chosen[np.argsort(tiebreak[chosen])]
    return out


def _expand_rings(conn, n_rings: int):
    """
    Neighbourhood graph extended to every node within `n_rings` hops.

    Parameters
    ----------
    conn : sparse matrix
        Spatial connectivity (direct neighbours, focal node excluded).
    n_rings : int
        Number of hops to include (1 returns the direct neighbours).

    Returns
    -------
    scipy.sparse.csr_matrix
        Binary connectivity without the focal node.
    """
    A = csr_matrix((conn > 0).astype(np.float64))
    reach, step = A.copy(), A.copy()
    for _ in range(n_rings - 1):
        step = csr_matrix((step @ A > 0).astype(np.float64))
        reach = csr_matrix(((reach + step) > 0).astype(np.float64))
    reach.setdiag(0)
    reach.eliminate_zeros()
    return reach



def _centred_ranks(cols: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Ranks of each column over spots, centred, and their Euclidean norms."""
    ranks = rankdata(cols, axis=0)
    ranks -= ranks.mean(axis=0)
    return ranks, np.sqrt((ranks**2).sum(axis=0))


def _signed_weights(cols: np.ndarray, is_up: np.ndarray, cv: np.ndarray) -> np.ndarray:
    """
    Signed gene weights by the `build_signed_weights` rule, vectorised.

    Parameters
    ----------
    cols : np.ndarray
        (spots x genes) expression of the gene set.
    is_up : np.ndarray
        Boolean per gene: True for the up genes, False for the down genes.
    cv : np.ndarray
        Coefficient of variation per gene (the weight magnitude).

    Returns
    -------
    np.ndarray
        CV x (Spearman correlation with the mean of the up genes - correlation with the mean of the
        down genes), divided by the size of the gene's own group.
    """
    ranks, norm = _centred_ranks(cols)

    def _spearman(reference):
        r = rankdata(reference)
        r -= r.mean()
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = (ranks * r[:, None]).sum(axis=0) / (norm * np.sqrt((r**2).sum()))
        return np.nan_to_num(rho)

    delta = _spearman(cols[:, is_up].mean(axis=1)) - _spearman(cols[:, ~is_up].mean(axis=1))
    return cv * delta / np.where(is_up, is_up.sum(), (~is_up).sum())


def _refit_signed_weights(X, ctrl_idx: np.ndarray, is_up: np.ndarray, cv: np.ndarray, rank_cache: dict) -> np.ndarray:
    """
    `_signed_weights` of every control set at once (equal up to floating-point rounding).

    A gene's ranks do not depend on the set, so they are computed once (and kept in `rank_cache`, shared by the
    signatures of one call, while it stays within _RANK_CACHE_ELEMENTS); per set only its two reference profiles
    (mean of its up and of its down genes) are ranked. Reference profiles and Spearman numerators are matrix
    products over all sets, taken in blocks of genes of bounded memory.

    Parameters
    ----------
    X : array or sparse matrix
        Expression (spots x genes), CSC if sparse.
    ctrl_idx : np.ndarray
        (sets x signature genes) gene indices of the control sets.
    is_up : np.ndarray
        Boolean per signature gene: up (True) or down (False).
    cv : np.ndarray
        CV of every gene (weight magnitude).
    rank_cache : dict
        Gene index -> (centred ranks, norm), filled as genes are ranked.

    Returns
    -------
    np.ndarray
        (sets x signature genes) refitted weights.
    """
    n_spots = X.shape[0]
    n_sets, k = ctrl_idx.shape
    genes, pos = np.unique(ctrl_idx, return_inverse=True)
    pos = pos.reshape(n_sets, k)
    set_of = np.repeat(np.arange(n_sets), k).reshape(n_sets, k)
    width = max(1, int(_RANK_BLOCK_ELEMENTS // n_spots))
    blocks = [np.arange(a, min(a + width, len(genes))) for a in range(0, len(genes), width)]

    def columns(b):
        c = X[:, genes[b]]
        return (c.toarray() if issparse(c) else np.asarray(c)).astype(np.float64)

    # reference profiles (sets x spots): mean of each set's up genes and of its down genes
    sides = (is_up, ~is_up)
    members = [
        csr_matrix((np.ones(n_sets * sel.sum()), (pos[:, sel].ravel(), set_of[:, sel].ravel())), shape=(len(genes), n_sets))
        for sel in sides
    ]
    sums = [np.zeros((n_spots, n_sets)) for _ in sides]
    kept = None
    for b in blocks:
        cb = columns(b)
        for s, m in zip(sums, members):
            s += np.asarray(cb @ m[b])
        kept = cb if len(blocks) == 1 else None
    both = np.empty((2 * n_sets, n_spots))  # centred ranks of the up references, then of the down references
    ref_norm = []
    for side, sel in enumerate(sides):
        r = both[side * n_sets : (side + 1) * n_sets]
        r[:] = rankdata((sums[side] / sel.sum()).T, axis=1)
        r -= r.mean(axis=1, keepdims=True)
        ref_norm.append(np.sqrt((r**2).sum(axis=1)))
    del sums

    # Spearman numerators of every gene against every reference, one gene block at a time
    num = np.empty((2, n_sets, k))
    norm = np.empty(len(genes))
    for b in blocks:
        new = [i for i in b if genes[i] not in rank_cache]
        cached = [i for i in b if genes[i] in rank_cache]
        ranks_b = np.empty((len(b), n_spots))
        for i in cached:
            ranks_b[i - b[0]], norm[i] = rank_cache[genes[i]]
        if new:
            cb = kept if kept is not None else columns(b)
            r_new, n_new = _centred_ranks(cb[:, np.asarray(new) - b[0]])
            r_new = np.ascontiguousarray(r_new.T)
            if (len(rank_cache) + len(new)) * n_spots <= _RANK_CACHE_ELEMENTS:
                rank_cache.update({genes[i]: (r_new[j], n_new[j]) for j, i in enumerate(new)})
            ranks_b[np.asarray(new) - b[0]] = r_new
            norm[new] = n_new
        prod = ranks_b @ both.T  # genes in block x (2 sets)
        in_b = (pos >= b[0]) & (pos <= b[-1])
        num[0][in_b] = prod[pos[in_b] - b[0], set_of[in_b]]
        num[1][in_b] = prod[pos[in_b] - b[0], n_sets + set_of[in_b]]

    delta = 0.0
    for side, sign in ((0, 1.0), (1, -1.0)):
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = num[side] / (norm[pos] * ref_norm[side][:, None])
        delta = delta + sign * np.nan_to_num(rho)
    return cv[ctrl_idx] * delta / np.where(is_up, is_up.sum(), (~is_up).sum())



def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    p = np.asarray(p, dtype=np.float64)
    n = len(p)
    order = np.argsort(p, kind="stable")
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def score(
    adata: AnnData,
    gene_set: list | dict | None = None,
    gene_weights: dict | None = None,
    score_key: str | list | None = None,
    spatial_key: str | None = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = False,
    batch_key: str | None = None,
    null_model: bool = True,
    n_control_sets: int | None = None,
    compute_q_vals: bool = False,
    gam_kwargs: dict | None = None,
    random_state: int | None = 0,
    delaunay: bool = False,
    n_rings: int = 1,
) -> None:
    """
    Compute spatially smoothed gene set enrichment scores, standardised against a null model of matched random gene sets, for one or more gene signatures.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix, containing expression values and spatial coordinates in `obsm`.

    gene_set : list or dict or None
        Gene set(s) to be scored. If a list is provided, it is interpreted as a single gene signature.
        If a dict is provided, keys are signature names and values are lists of gene symbols, or
        ``{"up": [...], "down": [...]}`` for a signed signature: its weights are then fitted to the data by
        the `build_signed_weights` rule and, with `null_model`, refitted by the same rule for every control
        set (pass signed programmes this way rather than as fitted `gene_weights`, whose signs the control
        sets would only inherit). If None, `gene_weights` must be provided and gene sets will be inferred
        from the keys of `gene_weights`.

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

    correct_spatial_covariates : bool, default False
        Whether to correct scores for spatial covariates using a GAM of the spatial coordinates.
        The GAM is fitted to each score itself, so it also removes broad biological patterns
        (large domains, gradients) and leaves negative halos around compact ones. With
        `null_model`, trends shared with random genes already cancel, so the GAM only removes
        the signature's own broad pattern.

    batch_key : str or None, optional
        Column in `adata.obs` indicating batch labels. Inferred gene weights (CV) then use the
        within-batch standard deviation, so a gene is not up-weighted for differing between
        batches; one weight per gene is still shared by all batches, which keeps scores
        comparable across them. The null model matches control genes on the same CV.

    null_model : bool, default True
        Whether to standardise each score against matched random gene sets. Each signature gene
        is replaced by a random gene of similar mean expression and coefficient of variation
        (signature genes excluded) that keeps its weight, `n_control_sets` times; any remaining
        difference in mean expression (e.g. for the most highly expressed genes, which have no
        exact peers) is removed by shifting each control set, and every control set is scored
        through the same steps. The stored score is z = (s - mean_null) / sd_null at
        each spot, with sd_null moderated towards its median. z > 0 means the signature scores
        above matched random genes at that spot. z is a competitive standardisation for ranking and
        display, not a calibrated p-value: the per-spot null is right-skewed and signature genes
        co-vary, so even random gene sets exceed z = 1.96 at roughly 3-5% of spots. For per-spot
        calls, use `compute_q_vals`.

    n_control_sets : int or None, default None
        Number of control gene sets per signature when `null_model` is True (at least 2). None uses 100,
        or 1000 with `compute_q_vals`. 100 gives z that agrees with 500 sets at Spearman 0.99; q-values
        need more, since the smallest attainable p is 1 / (n_control_sets + 1).

    compute_q_vals : bool, default False
        Whether to also compute per-spot empirical p-values and their Benjamini-Hochberg q-values
        against the control sets (requires `null_model`; uses 1000 control sets unless
        `n_control_sets` is given). The q-values are written to `adata.obs[f"{signature}_q"]`, and
        p- and q-values to `adata.uns["null_model"][signature]`. Call spots with q (e.g. q < 0.05),
        not z > 1.96. By default only the score is written to `adata.obs`.

    gam_kwargs : dict or None, optional
        Keyword arguments passed to `pygam.te(0, 1, ...)` when `correct_spatial_covariates` is True,
        e.g. ``{"n_splines": 4, "lam": 60}`` for a near-bilinear trend. Defaults to pyGAM's
        (10 splines per axis, lam = 0.6). Constrained terms fall back to a slower fit per column.

    random_state : int or None, default 0
        Non-negative seed for drawing control gene sets. Each signature's control sets depend only
        on this seed, its genes (and up/down labels), its weights and the data (not on its name, gene order,
        the column order of `adata`, the other signatures in the call, the storage format of `adata.X`
        or the numpy build).

    delaunay : bool, default False
        Whether to build the smoothing graph from a Delaunay triangulation of the coordinates instead of
        the `n_neighbors` nearest neighbours (`n_neighbors` is then ignored). Recommended for single-cell
        spatial data, where it adapts to local cell density.

    n_rings : int, default 1
        Number of graph hops included in the smoothing neighbourhood (1 = direct neighbours). With
        `delaunay=True`, `n_rings=2` gives about 20 neighbours per cell.

    Returns
    -------
    None
        Scores are stored in `adata.obs` and gene contributions in `adata.uns["gene_contributions"]`.
        With `null_model`, the per-spot null mean and moderated null SD are stored in
        `adata.uns["null_model"][signature]`. With `compute_q_vals`, so are the empirical p-value
        (1 + #{control >= signature}) / (n_control_sets + 1) and its Benjamini-Hochberg q-value (per
        batch, across spots), and the q-value is also written to `adata.obs[f"{signature}_q"]`.
        Per-spot calls should use q (e.g. q < 0.05), not z > 1.96: co-expressed gene sets exceed z = 1.96
        more often than random genes do, whereas q controls false discoveries against the decoys.

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

    With `null_model`, per-batch z-scoring (Eq. 2) is skipped: z is
    unchanged by any per-batch shift or scaling shared by a signature and
    its control sets, so the null model already puts batches on one
    scale, whereas z-scoring each column separately would divide the
    signature by its own (signal-bearing) spread.
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

    if not isinstance(n_rings, (int, np.integer)) or n_rings < 1:
        raise ValueError("n_rings must be an integer >= 1.")
    if batch_key is not None and adata.obs[batch_key].isna().any():
        raise ValueError(f"adata.obs['{batch_key}'] has missing batch labels.")
    if compute_q_vals and not null_model:
        raise ValueError("compute_q_vals=True needs null_model=True (q-values are computed against the control sets).")

    if null_model:
        if random_state is not None and (not isinstance(random_state, (int, np.integer)) or random_state < 0):
            raise ValueError("random_state must be a non-negative integer or None.")
        if not adata.var_names.is_unique:
            warnings.warn(
                "adata.var_names are not unique; columns sharing a name with a signature gene are "
                "excluded from the control pool.",
                UserWarning,
                stacklevel=2,
            )
        if n_control_sets is None:
            n_control_sets = 1000 if compute_q_vals else 100
        if not isinstance(n_control_sets, (int, np.integer)) or n_control_sets < 2:
            raise ValueError("n_control_sets must be an integer >= 2 when null_model=True.")
        if n_control_sets < 20:
            warnings.warn(
                "n_control_sets < 20: the per-spot null SD is noisy and z has heavy tails.",
                UserWarning,
                stacklevel=2,
            )
        if correct_spatial_covariates:
            warnings.warn(
                "With null_model=True, spatial trends shared with matched random genes already "
                "cancel, so the GAM only removes the signature's own broad pattern and leaves "
                "negative halos around domains that can pass z < -1.96. Consider "
                "correct_spatial_covariates=False.",
                UserWarning,
                stacklevel=2,
            )

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
    if null_model and any(
        w is not None and min(w.values(), default=0) < 0 < max(w.values(), default=0)
        for w in gene_weights.values()
    ):
        warnings.warn(
            "Signed gene_weights are inherited unchanged by the null model's control sets. If the signs "
            "were fitted to this data (e.g. build_signed_weights), the null is anti-conservative; pass "
            "the lists as gene_set={name: {'up': [...], 'down': [...]}} so control weights are refitted.",
            UserWarning,
            stacklevel=2,
        )

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
                    delaunay=delaunay,
                    key_added="spatial",
                )
                conn = adata_batch.obsp["spatial_connectivities"]
                if n_rings > 1:
                    conn = _expand_rings(conn, n_rings)
                # Spots without neighbours (e.g. cells with duplicated centroids, which the
                # Delaunay triangulation drops) keep their own value instead of becoming 0.
                isolated = np.flatnonzero(np.asarray(conn.sum(axis=1)).ravel() == 0)
                if len(isolated):
                    warnings.warn(
                        f"{len(isolated)} spot(s) have no spatial neighbours (e.g. duplicated "
                        f"coordinates); they keep their own unsmoothed score.",
                        UserWarning,
                        stacklevel=2,
                    )
                    conn = csr_matrix(conn) + csr_matrix(
                        (np.ones(len(isolated)), (isolated, isolated)), shape=conn.shape
                    )
            else:
                conn = None

            _spatial_cache[batch] = (mask, conn, coords)

    # The GAM hat matrix depends only on the coordinates: build the spline
    # design B and the penalised normal matrix BᵀB + P once per batch and
    # reuse them for every signature and control set (None if pygam
    # internals are unavailable or the terms are constrained, which makes
    # the fit non-linear; the per-column fallback is used then).
    _gam_cache = {}  # batch → (B, BᵀB + P) or None
    if correct_spatial_covariates:
        for batch in batches:
            _, _, coords = _spatial_cache[batch]
            try:
                gam = LinearGAM(te(0, 1, **(gam_kwargs or {}))).fit(coords, coords[:, 0])
                if getattr(gam.terms, "hasconstraint", False):
                    _gam_cache[batch] = None
                    continue
                B = gam._modelmat(coords)
                B = B.toarray() if hasattr(B, "toarray") else np.asarray(B)
                P = gam._P()
                P = P.toarray() if hasattr(P, "toarray") else np.asarray(P)
                _gam_cache[batch] = (B, B.T @ B + P)
            except Exception:
                _gam_cache[batch] = None

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
    sig_design = {}  # signature → (gene indices, weights, up/down flags or None), in input gene order

    pbar = tqdm(sig_names)
    for j, sig_name in enumerate(pbar):
        genes = gene_set[sig_name]
        is_up = None
        if isinstance(genes, dict):
            # signed signature: up genes, then down genes (a gene listed in both counts as up)
            up = [g for g in dict.fromkeys(genes.get("up", [])) if g in var_index]
            down = [g for g in dict.fromkeys(genes.get("down", [])) if g in var_index and g not in up]
            if not up or not down:
                raise ValueError(f"Signed signature '{sig_name}' needs up and down genes present in the dataset")
            common_genes = up + down
            is_up = np.array([True] * len(up) + [False] * len(down))
            genes = common_genes
        else:
            # input order (not set order) keeps results independent of PYTHONHASHSEED
            common_genes = [g for g in dict.fromkeys(genes) if g in var_index]
        pbar.set_description(
            f"Scoring {sig_name}: {len(common_genes)}/{len(genes)} genes found"
        )
        if len(common_genes) < 2:
            raise ValueError(
                f"Signature '{sig_name}' has fewer than two genes in the dataset"
            )

        # Determine gene weights
        if is_up is not None:
            cv = infer_gene_weights(adata, common_genes, batch_key=batch_key)
            cols = X[:, [var_index[g] for g in common_genes]]
            cols = cols.toarray() if issparse(cols) else np.asarray(cols, dtype=np.float64)
            w_signed = _signed_weights(cols.astype(np.float64), is_up, np.array([cv[g] for g in common_genes]))
            current_gene_weights = dict(zip(common_genes, w_signed))
        elif sig_name in gene_weights and gene_weights[sig_name] is not None:
            current_gene_weights = {
                g: gene_weights[sig_name].get(g, 1.0) for g in common_genes
            }
        else:
            inferred_gene_weights[sig_name] = infer_gene_weights(adata, common_genes, batch_key=batch_key)
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
        weighted = weighted.astype(np.float32)  # single precision is plenty for contributions and halves their memory
        contributions[sig_name] = {
            g: weighted[:, i] for i, g in enumerate(common_genes)
        }
        sig_design[sig_name] = (np.asarray(idx), w_vec, is_up)

    def _transform(mat: np.ndarray, standardise: bool = True) -> np.ndarray:
        """Steps 2-4 (Eq. 2-5) on a (spots x columns) score matrix, in place."""
        # Step 2 (Eq. 2): per-batch z-score normalisation (only with batch_key).
        if batch_key is not None and standardise:
            for batch in batches:
                mask = np.asarray(adata.obs[batch_key] == batch)
                block = mat[mask]
                mu_b = block.mean(axis=0)
                sigma_b = block.std(axis=0)
                nz = sigma_b > 0
                block[:, nz] = (block[:, nz] - mu_b[nz]) / sigma_b[nz]
                mat[mask] = block

        # Step 3 (Eq. 3): spatial smoothing on the full score matrix.
        if smoothing:
            for batch in batches:
                mask, conn, _ = _spatial_cache[batch]
                denom = np.maximum(np.asarray(conn.sum(axis=1)).ravel(), 1e-10)
                mat[mask] = conn.dot(mat[mask]) / denom[:, None]

        # Step 4 (Eq. 4-5): spatial-covariate correction.
        #
        # LinearGAM(te(0, 1)) with a fixed smoothing parameter is a linear
        # smoother: predict = H @ y with H = B (BᵀB + P)⁻¹ Bᵀ, where B is the
        # tensor-product spline design and P the penalty — both functions of
        # the coordinates only. We therefore residualise every column
        # against the same fit in one solve, which is identical to fitting a
        # GAM per signature (verified to < 1e-7) but built once per batch.
        if correct_spatial_covariates:
            for batch in batches:
                mask, _, coords = _spatial_cache[batch]
                block = mat[mask]
                try:
                    # Fast path: one solve for every column.
                    B, A = _gam_cache[batch]
                    coef = np.linalg.solve(A, B.T @ block)
                    mat[mask] = block - B @ coef
                except Exception:
                    # Robust fallback: fit a GAM per column (original
                    # behaviour) if pygam internals are unavailable.
                    for j in range(block.shape[1]):
                        g_ = LinearGAM(te(0, 1, **(gam_kwargs or {}))).fit(coords, block[:, j])
                        block[:, j] = block[:, j] - g_.predict(coords)
                    mat[mask] = block
        return mat

    _transform(scores_mat, standardise=not null_model)

    # ------------------------------------------------------------------
    # Null model: standardise against matched random gene sets.
    #
    # Each signature gene is replaced by a random non-signature gene of
    # similar mean expression and CV (see `_matched_candidates`), and the
    # replacement keeps the signature gene's weight, so
    # a control set is the same weighted mean over random genes of matched
    # expression level and variability. The control sets go through
    # exactly the same steps, so whatever the pipeline does to random genes
    # — smoothing, the GAM — is also done to the reference, and technical
    # patterns shared with random genes (library size, RNA quality
    # gradients) cancel. The score at each spot is then
    # z = (s - mean_null) / sd_null, with sd_null moderated towards its
    # per-batch median (limma-style), which keeps z continuous where the
    # control sets barely vary (sparse spots, low-plex panels).
    # ------------------------------------------------------------------
    if null_model:
        if "null_model" not in adata.uns:
            adata.uns["null_model"] = {}
        n_genes = adata.n_vars

        # Gene mean and CV over column chunks (no copy of the whole matrix), in
        # float64 with the same arithmetic for sparse and dense storage, rounded
        # so that rounding-level differences cannot reorder genes.
        batch_masks = [
            np.asarray(adata.obs[batch_key] == b) if batch_key else np.ones(adata.n_obs, bool)
            for b in batches
        ]
        # batch indicator (batches x spots) for per-batch gene means
        batch_ind = csr_matrix(np.vstack(batch_masks).astype(np.float64))
        batch_n = np.asarray(batch_ind.sum(axis=1)).ravel()
        gene_mean = np.empty(n_genes)
        gene_sd = np.empty(n_genes)
        gene_mean_batch = np.empty((len(batches), n_genes))
        for a0 in range(0, n_genes, 1024):
            sl = slice(a0, min(a0 + 1024, n_genes))
            if issparse(X):
                Xc = X[:, sl].astype(np.float64)
                m = np.asarray(Xc.sum(axis=0)).ravel() / adata.n_obs
                e2 = np.asarray(Xc.multiply(Xc).sum(axis=0)).ravel() / adata.n_obs
                mb = np.asarray((batch_ind @ Xc).todense()) / batch_n[:, None]
            else:
                Xc = np.asarray(X[:, sl], dtype=np.float64)
                m = Xc.sum(axis=0) / adata.n_obs
                e2 = np.square(Xc).sum(axis=0) / adata.n_obs
                mb = (batch_ind @ Xc) / batch_n[:, None]
            # with batches, the SD is within-batch (around each batch's mean), as for the gene weights
            var = e2 - (batch_n / adata.n_obs) @ np.square(mb) if batch_key else e2 - m**2
            gene_mean[sl], gene_sd[sl] = m, np.sqrt(np.maximum(var, 0))
            gene_mean_batch[:, sl] = mb
        gene_cv = np.divide(gene_sd, gene_mean, out=np.zeros(n_genes), where=gene_mean > 0)
        mean_rank = rankdata(np.round(gene_mean, 12))
        cv_rank = rankdata(np.round(gene_cv, 12))
        name_rank = rankdata(np.asarray(adata.var_names, dtype=str), method="ordinal")

        rank_cache = {}  # gene -> centred ranks, shared by the signed signatures of this call
        coarse = []  # signatures with too few candidate genes for a good null
        no_peers = []  # signatures whose genes have no peers of similar expression
        for j, sig_name in enumerate(sig_names):
            idx, w_vec, is_up = sig_design[sig_name]
            o = np.argsort(name_rank[idx], kind="stable")  # canonical order by gene name
            idx, w_vec = idx[o], w_vec[o]
            is_up = None if is_up is None else is_up[o]
            # Seed from the genes themselves (sorted names, plus up/down labels for signed
            # sets), not the signature's name: the same gene set always gets the same decoys.
            gene_key = ";".join(adata.var_names[idx].astype(str))
            if is_up is not None:
                gene_key += "|" + "".join("u" if f else "d" for f in is_up)
            rng = np.random.default_rng(
                None if random_state is None else [int(random_state), zlib.crc32(gene_key.encode())]
            )

            # Exclude every column named like a signature gene (duplicate var_names).
            free = np.flatnonzero(~adata.var_names.isin(adata.var_names[idx]))
            if len(free) == 0:
                raise ValueError(
                    f"Signature '{sig_name}' covers every gene, so no control genes are left; "
                    f"use null_model=False."
                )
            k = len(idx)
            # Draw the control sets; if they share many genes (pools of the most
            # highly expressed genes overlap), widen the pools and draw again.
            # Shared genes make the sets correlated, so the null SD would be
            # underestimated and |z| inflated.
            n_mean, n_cv = _N_MEAN_WINDOW, _N_MATCHED_CANDIDATES
            coarse_size = len(free) < max(5 * k, 50)
            while True:
                cand = _matched_candidates(mean_rank, cv_rank, idx, free, n_mean, n_cv, name_rank)
                n_cand = cand.shape[1]
                # One candidate per signature gene and set; redraw genes that
                # repeat within a set so each control set has distinct genes.
                ctrl_idx = cand[np.arange(k), rng.integers(0, n_cand, size=(n_control_sets, k))]
                dups_left = False
                for _ in range(100):
                    order = np.argsort(ctrl_idx, axis=1, kind="stable")
                    srt = np.take_along_axis(ctrl_idx, order, axis=1)
                    dup_sorted = np.zeros(srt.shape, dtype=bool)
                    dup_sorted[:, 1:] = srt[:, 1:] == srt[:, :-1]
                    if not dup_sorted.any():
                        break
                    dup = np.zeros_like(dup_sorted)
                    np.put_along_axis(dup, order, dup_sorted, axis=1)
                    rows, cols = np.nonzero(dup)
                    ctrl_idx[rows, cols] = cand[cols, rng.integers(0, n_cand, size=len(cols))]
                else:
                    dups_left = True  # no redraw round left every set free of repeats
                # mean Jaccard overlap between control sets
                used, inv = np.unique(ctrl_idx, return_inverse=True)
                inc = csr_matrix(
                    (np.ones(ctrl_idx.size), (np.repeat(np.arange(n_control_sets), k), inv.ravel())),
                    shape=(n_control_sets, len(used)),
                )
                inc.data[:] = 1.0
                shared = (inc @ inc.T).toarray()
                size = np.diag(shared)
                jaccard = shared / (size[:, None] + size[None, :] - shared)
                overlap = float(jaccard[~np.eye(n_control_sets, dtype=bool)].mean())
                if overlap <= _MAX_CONTROL_OVERLAP or coarse_size or n_mean >= len(free):
                    break
                n_mean, n_cv = 2 * n_mean, 2 * n_cv
            n_distinct = len({tuple(np.sort(r)) for r in ctrl_idx})
            if coarse_size or n_distinct < n_control_sets or dups_left or overlap > _MAX_CONTROL_OVERLAP:
                coarse.append(sig_name)

            # Control weights: inherited from the signature gene each control
            # gene replaces, or - for signed signatures, whose signs are fitted
            # to the data - refitted for every control set by the same rule, so
            # signature and controls go through an identical procedure.
            if is_up is None:
                w_ctrl = np.tile(w_vec, (n_control_sets, 1))
            else:
                w_ctrl = _refit_signed_weights(X, ctrl_idx, is_up, gene_cv, rank_cache)
            w_abs = np.abs(w_ctrl).sum(axis=1, keepdims=True)
            w_ctrl_norm = w_ctrl / np.where(w_abs > 0, w_abs, 1.0)

            # Step 1 for every control set at once: X[:, used genes] @ W,
            # W sparse (used genes x sets).
            W = csc_matrix(
                (w_ctrl_norm.ravel(), (inv.ravel(), np.repeat(np.arange(n_control_sets), k))),
                shape=(len(used), n_control_sets),
            )
            ctrl = X[:, used] @ W
            ctrl = ctrl.toarray() if issparse(ctrl) else np.asarray(ctrl, dtype=np.float64)
            # Shift each control set, within each batch, by its remaining mismatch
            # in mean expression, so every set has the signature's mean in that
            # batch and the null compares how expression is distributed in space,
            # not its overall level (per batch, so that shifts or scalings shared
            # by a batch still cancel). Smoothing and the GAM preserve or remove
            # constants, so the shift is applied first.
            w_norm = w_vec / np.abs(w_vec).sum()
            offsets = np.empty((len(batches), n_control_sets))
            for b, mask in enumerate(batch_masks):
                gm = gene_mean_batch[b]
                offsets[b] = (w_norm * gm[idx]).sum() - (w_ctrl_norm * gm[ctrl_idx]).sum(axis=1)
                ctrl[mask] += offsets[b][None, :]
            _transform(ctrl, standardise=False)
            # Unshifted controls (the GAM removes the constant shift itself): used where
            # neither signature nor controls are expressed, which must stay ties.
            ctrl0 = ctrl.copy()
            if not correct_spatial_covariates:
                for b, mask in enumerate(batch_masks):
                    ctrl0[mask] -= offsets[b][None, :]

            s = scores_mat[:, j].copy()
            mu = ctrl.mean(axis=1)
            sd = ctrl.std(axis=1, ddof=1)
            sd0 = ctrl0.std(axis=1, ddof=1)
            sd_mod = sd.copy()
            for mask in batch_masks:
                pos = sd0[mask] > 0
                s0_sq = np.median(sd[mask][pos] ** 2) if pos.any() else 1.0
                sd_mod[mask] = np.sqrt(
                    (_SD_PRIOR_DF * s0_sq + (n_control_sets - 1) * sd[mask] ** 2)
                    / (_SD_PRIOR_DF + n_control_sets - 1)
                )
            # "zero" up to the rounding left by subtracting the shift after smoothing
            ctrl0_zero = np.abs(ctrl0) <= 1e-9 * max(1.0, float(np.abs(offsets).max()))
            z = (s - mu) / sd_mod
            empty = (s == 0) & ctrl0_zero.all(axis=1)
            z[empty] = 0.0
            scores_mat[:, j] = z
            offset_ratio = float(np.abs(offsets.mean(axis=1)).max() / max(np.median(sd_mod), 1e-12))
            if offset_ratio > _MAX_OFFSET_RATIO:
                no_peers.append(sig_name)
            adata.uns["null_model"][sig_name] = {
                "null_mean": mu.astype(np.float32),
                "null_sd": sd_mod.astype(np.float32),
                "n_control_sets": n_control_sets,
                "mean_level_offset": offsets.mean(axis=1).astype(np.float32),
                "level_offset_ratio": offset_ratio,
                "control_overlap": overlap,
                "weights_refitted": is_up is not None,
            }
            if compute_q_vals:
                ge = (ctrl >= (s - 1e-12 * np.maximum(1, np.abs(s)))[:, None]) | ((s == 0)[:, None] & ctrl0_zero)
                p_emp = (1 + ge.sum(axis=1)) / (n_control_sets + 1)
                q_emp = np.empty_like(p_emp)
                for mask in batch_masks:
                    q_emp[mask] = _bh(p_emp[mask])
                adata.uns["null_model"][sig_name]["p_empirical"] = p_emp.astype(np.float32)
                adata.uns["null_model"][sig_name]["q_empirical"] = q_emp.astype(np.float32)

        if coarse:
            warnings.warn(
                f"Few distinct genes are available to build control sets for {len(coarse)} "
                f"signature(s) (e.g. '{coarse[0]}'), as on low-plex panels or for very large "
                f"signatures or sets of the most highly expressed genes (control sets overlapping "
                f"heavily); the null model is coarse there and z may be inflated.",
                UserWarning,
                stacklevel=2,
            )
        if no_peers:
            warnings.warn(
                f"{len(no_peers)} signature(s) (e.g. '{no_peers[0]}') are made of genes with no peers of "
                f"similar expression (e.g. the most highly expressed genes); the null model can match their "
                f"level but not their variability, so their z is descriptive rather than calibrated.",
                UserWarning,
                stacklevel=2,
            )

    # Write results back (one column per signature).
    for j, sig_name in enumerate(sig_names):
        adata.obs[f"{sig_name}_score"] = scores_mat[:, j]
        if compute_q_vals:
            adata.obs[f"{sig_name}_q"] = adata.uns["null_model"][sig_name]["q_empirical"]
        else:
            adata.obs.drop(columns=f"{sig_name}_q", errors="ignore", inplace=True)  # no stale q-values
        adata.uns["gene_contributions"][sig_name] = contributions[sig_name]
        if not null_model and "null_model" in adata.uns:
            adata.uns["null_model"].pop(sig_name, None)

    return None
