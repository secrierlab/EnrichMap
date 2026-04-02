from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from tqdm import tqdm

logging.getLogger("squidpy").setLevel(logging.WARNING)


def compare_wasserstein(
    adata: AnnData,
    score_key: str,
    batch_key: str,
    spatial_key: str = "spatial",
    spatial_weight: float = 1.0,
    score_weight: float = 1.0,
    n_subsample: int | None = 2000,
    n_permutations: int = 999,
    random_state: int = 0,
    group_key: str | None = None,
    plot: bool = True,
    figsize: tuple[float, float] = (7, 6),
    cmap: str = "magma_r",
    linkage_method: str = "average",
) -> pd.DataFrame:
    """
    Pairwise Wasserstein (earth mover's) distance between patients based on
    spatially embedded EnrichMap scores.

    Each patient's score field is represented as an empirical distribution
    over a joint (x, y, score) space. The Wasserstein-2 distance then
    quantifies how much "work" is needed to transform one patient's spatial
    score landscape into another's, capturing differences in both the
    spatial arrangement and the magnitude of scores simultaneously. Two
    patients can have identical marginal score distributions yet produce a
    large Wasserstein distance if their spatial patterns differ — for
    instance, a single coherent hotspot versus many scattered foci.

    The output is a patient-by-patient distance matrix that can be used
    directly for hierarchical clustering, multidimensional scaling or
    downstream statistical testing.

    Normalisation
    ~~~~~~~~~~~~~
    Spatial coordinates are min-max normalised per patient to [0, 1] so
    that varying slide extents (e.g. different tissue sizes or imaging
    resolutions) do not dominate the distance. Scores are likewise
    normalised per patient to [0, 1]. The ``spatial_weight`` and
    ``score_weight`` parameters allow tuning the relative contribution of
    location versus score amplitude to the final distance:

    - ``spatial_weight > score_weight``: emphasises where scores are
      located in tissue space; useful for detecting pattern rearrangements.
    - ``score_weight > spatial_weight``: emphasises score magnitudes;
      useful when amplitude differences are biologically meaningful.
    - Equal weights (default): balanced comparison.

    Statistical testing
    ~~~~~~~~~~~~~~~~~~~
    When exactly **two patients** are present, a spot-label permutation
    test is run automatically. All spots from both patients are pooled and
    randomly split into two groups of the original sizes. Coordinates are
    preserved — each spot keeps its original spatial position — so the test
    asks: "is the observed Wasserstein distance larger than expected if
    these score values were randomly distributed across these two tissue
    architectures?" The p-value is computed as the fraction of permuted
    distances that equal or exceed the observed distance.

    When **multiple patients** are present and ``group_key`` is provided,
    a PERMANOVA (permutational multivariate analysis of variance) is run
    on the distance matrix, analogous to ``vegan::adonis2`` in R. This
    tests whether within-group distances are smaller than between-group
    distances — i.e. whether the clinical grouping explains a significant
    proportion of variance in spatial score organisation.

    All test results are stored in ``df.attrs`` and displayed in the plot
    title when ``plot=True``.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must contain the EnrichMap score column in
        ``adata.obs`` and spatial coordinates in ``adata.obsm[spatial_key]``.

    score_key : str
        Column name in ``adata.obs`` holding the EnrichMap score to analyse,
        e.g. ``"enrichmap_score"`` or ``"EMT_score"``.

    batch_key : str
        Column name in ``adata.obs`` identifying individual patients or
        slides, e.g. ``"patient_id"`` or ``"library_id"``. Each unique
        value is treated as a separate sample.

    spatial_key : str, default ``"spatial"``
        Key in ``adata.obsm`` containing the 2D spatial coordinates as an
        (n_spots, 2) array.

    spatial_weight : float, default 1.0
        Multiplicative weight applied to the (normalised) spatial
        dimensions before computing pairwise distances. Increase to make
        the metric more sensitive to where scores are located; decrease
        to downweight spatial location relative to score amplitude.

    score_weight : float, default 1.0
        Multiplicative weight applied to the (normalised) score dimension.
        Increase to make the metric more sensitive to score amplitude
        differences; decrease to emphasise spatial arrangement.

    n_subsample : int or None, default 2000
        If set, subsample each patient to at most this many spots before
        computing distances. The optimal transport solver operates on an
        (n × m) cost matrix, so runtime is roughly O(n² · m) for each
        pair. Subsampling to 2000 keeps pairwise computation under ~1s
        for typical Visium data. Set to ``None`` to use all spots (may be
        slow for >5000 spots per sample).

    n_permutations : int, default 999
        Number of permutations for significance testing. Higher values
        give more precise p-values. For exploratory analysis 499 is
        sufficient; for publication-ready results use 999 or higher.

    random_state : int, default 0
        Seed for the random number generator, used for subsampling and
        permutation tests. Set for reproducibility.

    group_key : str or None, optional
        Column name in ``adata.obs`` for a higher-level clinical grouping,
        e.g. ``"subtype"``, ``"treatment_arm"`` or ``"response"``. When
        provided, the clustermap is annotated with a colour sidebar and
        a PERMANOVA test is run on the distance matrix. When ``None``,
        no group-level testing is performed.

    plot : bool, default True
        Whether to produce a hierarchically clustered heatmap of the
        distance matrix. For exactly two patients this is not very
        informative; consider setting ``plot=False``.

    figsize : tuple of float, default (7, 6)
        Figure size in inches (width, height) for the clustermap.

    cmap : str, default ``"magma_r"``
        Colourmap for the distance heatmap. ``"magma_r"`` gives a dark
        (low distance) to light (high distance) gradient.

    linkage_method : str, default ``"average"``
        Linkage method for hierarchical clustering of the distance matrix.
        ``"average"`` (UPGMA) is the standard choice for distance matrices.
        Other options include ``"ward"``, ``"complete"`` and ``"single"``.

    Returns
    -------
    pd.DataFrame
        Square distance matrix indexed and columned by patient name.
        Values are Wasserstein-2 distances in the joint
        (x_norm, y_norm, score_norm) space.

        Statistical test results, when applicable, are stored as
        dictionaries in ``df.attrs``:

        - ``df.attrs["pairwise_test"]``: spot-label permutation test
          result (two-patient case), containing keys ``"observed_distance"``,
          ``"null_mean"``, ``"null_std"``, ``"p_value"`` and
          ``"n_permutations"``.
        - ``df.attrs["permanova"]``: PERMANOVA result (multi-patient with
          ``group_key``), containing keys ``"pseudo_F"``, ``"p_value"``
          and ``"n_permutations"``.

    Examples
    --------
    Pairwise comparison of two patients:

    >>> dist = compare_wasserstein(
    ...     adata,
    ...     score_key="enrichmap_score",
    ...     batch_key="patient_id",
    ...     plot=False,
    ... )
    >>> print(dist)
                patient_01  patient_02
    patient_01      0.0000      0.1375
    patient_02      0.1375      0.0000
    >>> print(dist.attrs["pairwise_test"]["p_value"])
    0.025

    Multi-patient comparison with PERMANOVA:

    >>> dist = compare_wasserstein(
    ...     adata,
    ...     score_key="enrichmap_score",
    ...     batch_key="patient_id",
    ...     group_key="subtype",
    ... )
    >>> print(dist.attrs["permanova"])
    {'test': 'PERMANOVA', 'pseudo_F': 6.22, 'p_value': 0.088, ...}

    Emphasise spatial arrangement over score amplitude:

    >>> dist = compare_wasserstein(
    ...     adata,
    ...     score_key="enrichmap_score",
    ...     batch_key="patient_id",
    ...     spatial_weight=2.0,
    ...     score_weight=0.5,
    ... )

    See Also
    --------
    compare_morans_i : Spatial autocorrelation comparison via Moran's I.
    compare_variograms : Semivariogram-based comparison of spatial scale
        and structure.
    """
    try:
        import ot
    except ImportError:
        raise ImportError("The POT package is required. Install with: pip install POT")

    if score_key not in adata.obs.columns:
        raise KeyError(f"'{score_key}' not found in adata.obs")

    rng = np.random.default_rng(random_state)
    patients = adata.obs[batch_key].unique()

    def _embed(ad_sub, rng_local):
        """Normalise and embed a single patient's spots."""
        coords = ad_sub.obsm[spatial_key].copy().astype(np.float64)
        scores = ad_sub.obs[score_key].values.astype(np.float64)

        for d in range(coords.shape[1]):
            dmin, dmax = coords[:, d].min(), coords[:, d].max()
            rng_d = dmax - dmin
            coords[:, d] = (coords[:, d] - dmin) / rng_d if rng_d > 0 else 0.0

        smin, smax = scores.min(), scores.max()
        if smax - smin > 0:
            scores = (scores - smin) / (smax - smin)

        emb = np.column_stack([coords * spatial_weight, scores[:, None] * score_weight])
        if n_subsample is not None and emb.shape[0] > n_subsample:
            idx = rng_local.choice(emb.shape[0], n_subsample, replace=False)
            emb = emb[idx]
        return emb

    def _w2(emb_a, emb_b):
        """Wasserstein-2 distance between two embeddings."""
        a = np.ones(emb_a.shape[0]) / emb_a.shape[0]
        b = np.ones(emb_b.shape[0]) / emb_b.shape[0]
        M = ot.dist(emb_a, emb_b, metric="sqeuclidean")
        return np.sqrt(max(ot.emd2(a, b, M), 0))

    embeddings = {}
    for patient in patients:
        mask = adata.obs[batch_key] == patient
        embeddings[patient] = _embed(adata[mask].copy(), rng)

    # Pairwise distances
    n_pat = len(patients)
    dist_matrix = np.zeros((n_pat, n_pat))
    pairs = [(i, j) for i in range(n_pat) for j in range(i + 1, n_pat)]

    for i, j in tqdm(pairs, desc="Wasserstein distances"):
        d = _w2(embeddings[patients[i]], embeddings[patients[j]])
        dist_matrix[i, j] = d
        dist_matrix[j, i] = d

    result = pd.DataFrame(dist_matrix, index=patients, columns=patients)

    # Pairwise permutation test (2 patients)
    if n_pat == 2:
        result.attrs["pairwise_test"] = _pairwise_permutation_wasserstein(
            adata,
            score_key,
            batch_key,
            spatial_key,
            patients,
            spatial_weight,
            score_weight,
            n_subsample,
            n_permutations=n_permutations,
            random_state=random_state,
            observed_dist=dist_matrix[0, 1],
            _embed_fn=_embed,
            _w2_fn=_w2,
        )

    # PERMANOVA (multiple patients with group_key)
    if group_key is not None and group_key in adata.obs.columns:
        sample_groups = (
            adata.obs.groupby(batch_key)[group_key].first().reindex(patients)
        )
        if sample_groups.nunique() == 2:
            result.attrs["permanova"] = _permanova(
                dist_matrix,
                sample_groups.values,
                n_permutations=n_permutations,
                random_state=random_state,
            )

    if plot:
        _plot_wasserstein(
            result,
            adata=adata,
            batch_key=batch_key,
            group_key=group_key,
            figsize=figsize,
            cmap=cmap,
            linkage_method=linkage_method,
        )

    return result


def _pairwise_permutation_wasserstein(
    adata,
    score_key,
    batch_key,
    spatial_key,
    patients,
    spatial_weight,
    score_weight,
    n_subsample,
    n_permutations,
    random_state,
    observed_dist,
    _embed_fn,
    _w2_fn,
):
    """
    Spot-label permutation test for Wasserstein distance between two patients.

    Constructs a null distribution by pooling all spots from both patients
    and randomly reassigning them to two groups of the original sizes.
    Each spot retains its spatial coordinates, so the null tests whether
    the observed distance is larger than expected if score values were
    randomly distributed across the two tissue architectures.

    This is an internal function called automatically by
    :func:`compare_wasserstein` when exactly two patients are present.

    Parameters
    ----------
    adata : AnnData
        Full annotated data matrix (both patients).
    score_key, batch_key, spatial_key : str
        Keys passed through from the parent call.
    patients : array-like
        Two-element array of patient identifiers.
    spatial_weight, score_weight : float
        Dimension weights passed through from the parent call.
    n_subsample : int or None
        Subsampling cap per pseudo-patient in each permutation.
    n_permutations : int
        Number of random label reassignments.
    random_state : int
        Random seed.
    observed_dist : float
        The observed Wasserstein distance between the two real patients.
    _embed_fn : callable
        The embedding function (normalise + weight + subsample).
    _w2_fn : callable
        The Wasserstein-2 distance function.

    Returns
    -------
    dict
        Test result with keys: ``"test"``, ``"patient_a"``,
        ``"patient_b"``, ``"observed_distance"``, ``"null_mean"``,
        ``"null_std"``, ``"p_value"``, ``"n_permutations"``.
    """
    rng = np.random.default_rng(random_state)
    mask_a = (adata.obs[batch_key] == patients[0]).values
    mask_b = (adata.obs[batch_key] == patients[1]).values
    n_a = mask_a.sum()

    pooled_idx = np.where(mask_a | mask_b)[0]
    null_dists = np.empty(n_permutations)

    for i in range(n_permutations):
        perm = rng.permutation(pooled_idx)
        idx_a, idx_b = perm[:n_a], perm[n_a:]
        emb_a = _embed_fn(adata[idx_a].copy(), rng)
        emb_b = _embed_fn(adata[idx_b].copy(), rng)
        null_dists[i] = _w2_fn(emb_a, emb_b)

    p = (np.sum(null_dists >= observed_dist) + 1) / (n_permutations + 1)

    return {
        "test": "spot-label permutation (Wasserstein distance)",
        "patient_a": patients[0],
        "patient_b": patients[1],
        "observed_distance": observed_dist,
        "null_mean": null_dists.mean(),
        "null_std": null_dists.std(),
        "p_value": p,
        "n_permutations": n_permutations,
    }


def _permanova(dist_matrix, group_labels, n_permutations=999, random_state=0):
    """
    PERMANOVA (permutational multivariate analysis of variance) on a
    distance matrix.

    Tests whether the centroids of two (or more) groups in distance space
    differ more than expected by chance, by partitioning total sum of
    squared distances into within-group and between-group components and
    computing a pseudo-F statistic. Significance is assessed by permuting
    group labels and recomputing pseudo-F, yielding a distribution-free
    p-value. Equivalent to ``vegan::adonis2`` in R.

    This is an internal function called automatically by
    :func:`compare_wasserstein` when ``group_key`` is provided and has
    two or more levels.

    Parameters
    ----------
    dist_matrix : np.ndarray
        Square (n_samples, n_samples) distance matrix.
    group_labels : array-like
        Group label for each sample, in the same order as the rows/columns
        of ``dist_matrix``.
    n_permutations : int, default 999
        Number of label permutations for the null distribution.
    random_state : int, default 0
        Random seed.

    Returns
    -------
    dict
        Test result with keys: ``"test"``, ``"pseudo_F"``, ``"p_value"``,
        ``"n_permutations"``.
    """
    rng = np.random.default_rng(random_state)
    groups = np.unique(group_labels)
    n = len(group_labels)

    def _pseudo_f(D, labels):
        sq_D = D**2
        SS_T = sq_D.sum() / (2 * n)
        SS_W = 0
        for g in groups:
            idx = np.where(labels == g)[0]
            n_g = len(idx)
            if n_g > 1:
                SS_W += sq_D[np.ix_(idx, idx)].sum() / (2 * n_g)
        SS_A = SS_T - SS_W
        k = len(groups)
        denom = SS_W / (n - k) if (n - k) > 0 else 1e-10
        return (SS_A / (k - 1)) / denom if denom > 0 else 0.0

    obs_F = _pseudo_f(dist_matrix, group_labels)

    null_Fs = np.empty(n_permutations)
    for i in range(n_permutations):
        perm_labels = rng.permutation(group_labels)
        null_Fs[i] = _pseudo_f(dist_matrix, perm_labels)

    p = (np.sum(null_Fs >= obs_F) + 1) / (n_permutations + 1)

    return {
        "test": "PERMANOVA",
        "pseudo_F": obs_F,
        "p_value": p,
        "n_permutations": n_permutations,
    }


def _plot_wasserstein(
    dist_df: pd.DataFrame,
    adata: AnnData | None = None,
    batch_key: str | None = None,
    group_key: str | None = None,
    figsize: tuple[float, float] = (6, 6),
    cmap: str = "magma_r",
    linkage_method: str = "average",
) -> sns.matrix.ClusterGrid:
    """
    Hierarchically clustered heatmap of the pairwise Wasserstein distance
    matrix.

    Produces a ``seaborn.clustermap`` with dendrograms on both axes. When
    ``group_key`` is provided, a colour sidebar is added to annotate each
    patient's clinical group, making it easy to assess whether the
    hierarchical clustering recovers the expected grouping.

    Statistical test results stored in ``dist_df.attrs`` (spot-label
    permutation p-value and/or PERMANOVA pseudo-F and p-value) are
    displayed in the figure title when present.

    Parameters
    ----------
    dist_df : pd.DataFrame
        Square distance matrix indexed by patient name, as returned by
        :func:`compare_wasserstein`. May carry ``dist_df.attrs["pairwise_test"]``
        and/or ``dist_df.attrs["permanova"]`` dictionaries.

    adata : AnnData or None
        Original annotated data matrix. Only needed when ``group_key`` is
        set, to look up each patient's group label.

    batch_key : str or None
        Column in ``adata.obs`` identifying patients. Required when
        ``group_key`` is set.

    group_key : str or None
        Column in ``adata.obs`` for group annotation. When provided, a
        colour sidebar is drawn alongside the heatmap indicating each
        patient's group membership.

    figsize : tuple of float, default (6, 6)
        Figure size in inches.

    cmap : str, default ``"magma_r"``
        Colourmap for the heatmap.

    linkage_method : str, default ``"average"``
        Linkage method for hierarchical clustering. ``"average"`` (UPGMA)
        is the standard choice for distance matrices.

    Returns
    -------
    seaborn.matrix.ClusterGrid
        The clustermap object, which provides access to the figure
        (``g.figure``) and individual axes (``g.ax_heatmap``,
        ``g.ax_row_dendrogram``, etc.).
    """
    row_colors = None
    if group_key is not None and adata is not None and batch_key is not None:
        if group_key in adata.obs.columns:
            sample_groups = (
                adata.obs.groupby(batch_key)[group_key].first().reindex(dist_df.index)
            )
            unique_groups = sample_groups.unique()
            pal = dict(
                zip(unique_groups, sns.color_palette("Set2", len(unique_groups)))
            )
            row_colors = sample_groups.map(pal)
            row_colors.name = group_key

    g = sns.clustermap(
        dist_df,
        method=linkage_method,
        cmap=cmap,
        figsize=figsize,
        row_colors=row_colors,
        col_colors=row_colors,
        linewidths=0.3,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Wasserstein distance"},
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    title_parts = []
    if "pairwise_test" in dist_df.attrs:
        p = dist_df.attrs["pairwise_test"]["p_value"]
        title_parts.append(f"Spot-label permutation p = {p:.4f}")
    if "permanova" in dist_df.attrs:
        res = dist_df.attrs["permanova"]
        title_parts.append(
            f"PERMANOVA F = {res['pseudo_F']:.2f}, p = {res['p_value']:.4f}"
        )
    if title_parts:
        g.figure.suptitle("  |  ".join(title_parts), fontsize=9, y=1.02)

    plt.tight_layout()
    return g
