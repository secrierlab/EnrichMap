from __future__ import annotations

import logging

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram as _dendrogram
from scipy.cluster.hierarchy import leaves_list, linkage
from tqdm import tqdm

logging.getLogger("squidpy").setLevel(logging.WARNING)

plt.rcParams["axes.grid"] = False


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
    cmap: str = "magma",
    linkage_method: str = "average",
    save: str | None = None,
    save_kwargs: dict | None = None,
    return_result: bool = False,
) -> pd.DataFrame | None:
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
        provided, the heatmap is annotated with colour strips and a PERMANOVA
        test is run on the distance matrix. When ``None``, no group-level
        testing is performed.

    plot : bool, default True
        Whether to produce a hierarchically clustered heatmap of the
        distance matrix. For exactly two patients this is not very
        informative; consider setting ``plot=False``.

    figsize : tuple of float, default (7, 6)
        Controls the figure width in inches. The height is auto-computed
        from the width so that heatmap cells are exactly square; the
        height component of this tuple is therefore ignored.

    cmap : str, default ``"magma"``
        Colourmap for the distance heatmap.

    linkage_method : str, default ``"average"``
        Linkage method for hierarchical clustering of the distance matrix.
        ``"average"`` (UPGMA) is the standard choice for distance matrices.
        Other options include ``"ward"``, ``"complete"`` and ``"single"``.

    save : str or None, optional
        Filename to save the figure (e.g. ``"wasserstein.pdf"``). The file
        is written to a ``figures/`` subdirectory of the working directory.
        When ``None`` (default), the figure is not saved.

    save_kwargs : dict or None, optional
        Extra keyword arguments forwarded to ``fig.savefig``, e.g.
        ``{"dpi": 300, "bbox_inches": "tight"}``.

    return_result : bool, default False
        When ``True``, return the distance matrix DataFrame. When ``False``
        (default), return ``None``.

    Returns
    -------
    pd.DataFrame or None
        When ``return_result=True``: square distance matrix indexed and
        columned by patient name. Values are Wasserstein-2 distances in
        the joint (x_norm, y_norm, score_norm) space.

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
            save=save,
            save_kwargs=save_kwargs,
        )

    if return_result:
        return result
    return None


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
    save: str | None = None,
    save_kwargs: dict | None = None,
) -> plt.Figure:
    """
    Hierarchically clustered heatmap of the pairwise Wasserstein distance
    matrix.

    Rows and columns are reordered by hierarchical clustering (UPGMA by
    default). Dendrograms are drawn using ``no_plot=True`` so their
    coordinates can be rescaled to ``[0, n]`` before the axes are shared
    with the heatmap via ``sharex`` / ``sharey``, guaranteeing pixel-perfect
    alignment regardless of ``tight_layout``.

    Row labels are placed on the right-hand side of the heatmap since the
    left side is occupied by the dendrogram.

    Parameters
    ----------
    dist_df : pd.DataFrame
        Square distance matrix indexed by patient name.
    adata : AnnData or None
        Needed only when ``group_key`` is set.
    batch_key : str or None
        Column in ``adata.obs`` identifying patients.
    group_key : str or None
        Column in ``adata.obs`` for group annotation colour strips.
    figsize : tuple of float, default (6, 6)
    cmap : str, default ``"magma_r"``
    linkage_method : str, default ``"average"``

    Returns
    -------
    matplotlib.figure.Figure
    """
    n = len(dist_df)
    labels = list(dist_df.index)

    # --- hierarchical clustering ---
    condensed = dist_df.values[np.triu_indices(n, k=1)]
    Z = linkage(condensed, method=linkage_method)
    order = leaves_list(Z)
    ordered_labels = [labels[i] for i in order]
    mat = dist_df.iloc[order].iloc[:, order]

    # Compute dendrogram data once (no_plot=True) in default orientation='top'.
    # scipy places leaf i at x = 10*i + 5, so dividing by 10 maps leaves to
    # [0.5, 1.5, …, n-0.5] — exactly the cell centres of the heatmap.
    # Heights (dcoord) are used as-is for the top dendrogram (y growing up)
    # and negated for the left dendrogram (x growing left).
    ddata = _dendrogram(Z, no_plot=True, color_threshold=0)
    max_h = max(max(ys) for ys in ddata["dcoord"])

    # --- group colour strips ---
    group_colors: list | None = None
    legend_handles: list = []
    if group_key is not None and adata is not None and batch_key is not None:
        if group_key in adata.obs.columns:
            sample_groups = (
                adata.obs.groupby(batch_key)[group_key]
                .first()
                .reindex(dist_df.index)
            )
            unique_groups = sorted(sample_groups.dropna().unique())
            pal = dict(
                zip(unique_groups, sns.color_palette("Set2", len(unique_groups)))
            )
            group_colors = [pal[sample_groups[s]] for s in ordered_labels]
            legend_handles = [
                mpatches.Patch(color=c, label=g) for g, c in pal.items()
            ]

    has_groups = group_colors is not None

    # --- GridSpec layout (sizes in inches) ---
    # figsize[0] sets the figure width; height is AUTO-COMPUTED so that the
    # heatmap cells are exactly square (hm_size × hm_size).
    DEND = 0.7    # dendrogram thickness
    STRIP = 0.12  # colour strip thickness
    GAP = 0.45    # gap between heatmap right edge and colorbar (holds row labels)
    CBAR = 0.25   # colour bar width

    fig_w = figsize[0]
    strip = STRIP if has_groups else 0.0

    # heatmap occupies a square block; compute from available width
    hm_size = fig_w - DEND - strip - GAP - CBAR
    fig_h = hm_size + DEND + strip   # height auto-set to match

    # cols:  left_dend | [left_strip] | heatmap | gap | cbar
    col_widths = [DEND] + ([strip] if has_groups else []) + [hm_size, GAP, CBAR]

    # rows:  top_dend | [top_strip] | heatmap
    row_heights = [DEND] + ([strip] if has_groups else []) + [hm_size]

    hm_row = len(row_heights) - 1
    hm_col = len(col_widths) - 3   # before gap and cbar
    cbar_col = len(col_widths) - 1  # last column

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        len(row_heights),
        len(col_widths),
        width_ratios=col_widths,
        height_ratios=row_heights,
        hspace=0.01,
        wspace=0.01,
    )

    # --- heatmap (created first so dendrograms can share its axes) ---
    ax_hm = fig.add_subplot(gs[hm_row, hm_col])
    ax_cbar = fig.add_subplot(gs[hm_row, cbar_col])

    sns.heatmap(
        mat,
        ax=ax_hm,
        cmap=cmap,
        xticklabels=ordered_labels,
        yticklabels=ordered_labels,
        linewidths=0,
        cbar=True,
        cbar_ax=ax_cbar,
        cbar_kws={"label": "Wasserstein distance"},
    )
    ax_hm.set_xlabel("")
    ax_hm.set_ylabel("")
    ax_hm.grid(False)
    ax_hm.tick_params(axis="x", rotation=90)
    ax_hm.yaxis.tick_right()
    ax_hm.tick_params(axis="y", rotation=0)

    # --- top dendrogram — shares x with heatmap ---
    # x: leaf i → i+0.5 (via /10); y: heights growing upward
    ax_top = fig.add_subplot(gs[0, hm_col], sharex=ax_hm)
    for xs, ys in zip(ddata["icoord"], ddata["dcoord"]):
        ax_top.plot([x / 10 for x in xs], ys, color="#888888", lw=0.8)
    ax_top.set_ylim(0, max_h * 1.05)
    ax_top.axis("off")

    # --- left dendrogram — shares y with heatmap ---
    # The heatmap y-axis is inverted (row 0 at top, y=0..n with y=0 at top).
    # icoord leaf positions /10 give y = i+0.5, matching heatmap row i.
    # Heights become negative x so the dendrogram grows leftward.
    ax_left = fig.add_subplot(gs[hm_row, 0], sharey=ax_hm)
    for xs, ys in zip(ddata["icoord"], ddata["dcoord"]):
        ax_left.plot([-y for y in ys], [x / 10 for x in xs], color="#888888", lw=0.8)
    ax_left.set_xlim(-max_h * 1.05, 0)
    ax_left.axis("off")

    # --- colour strips (share axes with heatmap for automatic alignment) ---
    if has_groups:
        strip_row = 1
        strip_col_idx = 1

        ax_ts = fig.add_subplot(gs[strip_row, hm_col], sharex=ax_hm)
        for i, color in enumerate(group_colors):
            ax_ts.add_patch(plt.Rectangle((i, 0), 1, 1, color=color, linewidth=0))
        ax_ts.set_ylim(0, 1)
        ax_ts.axis("off")

        ax_ls = fig.add_subplot(gs[hm_row, strip_col_idx], sharey=ax_hm)
        for i, color in enumerate(group_colors):
            ax_ls.add_patch(plt.Rectangle((0, i), 1, 1, color=color, linewidth=0))
        ax_ls.set_xlim(0, 1)
        ax_ls.axis("off")

    # --- group legend ---
    if legend_handles:
        ax_hm.legend(
            handles=legend_handles,
            title=group_key,
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0,
            frameon=False,
        )

    # --- title ---
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
        fig.suptitle("  |  ".join(title_parts), fontsize=9, y=1.01)

    plt.tight_layout()

    if save is not None:
        if save_kwargs is None:
            save_kwargs = {}
        save_path = Path("figures") / save
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, **save_kwargs)

    return fig
