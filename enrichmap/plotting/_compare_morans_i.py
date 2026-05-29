from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq
import seaborn as sns
from anndata import AnnData
from tqdm import tqdm
from pathlib import Path

logging.getLogger("squidpy").setLevel(logging.WARNING)

plt.rcParams["axes.grid"] = False


def compare_morans_i(
    adata: AnnData,
    score_key: str,
    batch_key: str,
    n_neighbors: int = 6,
    n_permutations: int = 999,
    random_state: int = 0,
    group_key: str | None = None,
    plot: bool = True,
    figsize: tuple[float, float] = (5, 4),
    palette: str | dict | None = None,
    ax: plt.Axes | None = None,
    save: str | None = None,
    save_kwargs: dict | None = None,
) -> pd.DataFrame:
    """
    Compare spatial autocorrelation of EnrichMap scores across patients
    using permutation-standardised Moran's I.

    Moran's I measures the degree of spatial autocorrelation in a score
    field: values near +1 indicate strong positive autocorrelation
    (similar scores cluster together), values near 0 indicate spatial
    randomness, and values near −1 indicate a chequerboard-like pattern
    where neighbours tend to have dissimilar scores.

    However, raw Moran's I values depend on the spatial weights matrix
    (i.e. the tissue architecture), so they are not directly comparable
    across slides with different geometries, spot densities or tissue
    shapes. This function addresses that by computing a **permutation
    z-score** for each patient: the observed I is standardised against a
    null distribution built by randomly permuting scores on that patient's
    own spatial graph. A z-score of 25 means "the score field is 25
    standard deviations more spatially clustered than expected under
    random placement on this particular tissue", and that statement is
    valid regardless of graph topology.

    The function reports both the raw Moran's I (useful for understanding
    the absolute autocorrelation) and the z-score (useful for cross-patient
    comparison). The per-patient permutation p-value indicates whether each
    individual patient's scores are significantly spatially autocorrelated.

    Statistical testing
    ~~~~~~~~~~~~~~~~~~~
    When exactly **two patients** are present, a spot-label permutation
    test is run automatically on the difference in raw Moran's I values.
    All spots from both patients are pooled and randomly reassigned to two
    groups of the original sizes; spatial graphs are rebuilt and Moran's I
    is recomputed for each pseudo-patient. This tests whether the observed
    difference in spatial autocorrelation is larger than expected if score
    values were randomly distributed across the two tissue architectures.

    When **multiple patients** are present and ``group_key`` is provided,
    a permutation test is run on the difference in group-mean z-scores.
    Patient-level z-scores are pooled, group labels are permuted, and the
    difference in means is recomputed to build a null distribution. This
    is analogous to a two-sample t-test but without distributional
    assumptions and works with any number of patients per group, including
    n=1.

    All test results are stored in ``df.attrs`` and displayed in the plot
    title when ``plot=True``.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Must contain the EnrichMap score column in
        ``adata.obs`` and spatial coordinates in ``adata.obsm`` (used
        internally by squidpy to build the spatial neighbours graph).

    score_key : str
        Column name in ``adata.obs`` holding the EnrichMap score to
        analyse, e.g. ``"enrichmap_score"`` or ``"EMT_score"``.

    batch_key : str
        Column name in ``adata.obs`` identifying individual patients or
        slides, e.g. ``"patient_id"`` or ``"library_id"``. Each unique
        value is treated as a separate sample with its own spatial graph.

    n_neighbors : int, default 6
        Number of nearest spatial neighbours used to construct the
        connectivity graph for each patient. For Visium data, 6
        corresponds to the immediate hexagonal ring. Increasing this
        value smooths the autocorrelation estimate over a larger
        neighbourhood but may dilute local signal.

    n_permutations : int, default 999
        Number of random permutations for constructing the per-patient
        null distribution and for the between-patient/group tests. Higher
        values give more precise p-values. For exploratory analysis 499
        is sufficient; for publication-ready results use 999 or higher.
        The per-patient p-value resolution is bounded by
        1 / (n_permutations + 1).

    random_state : int, default 0
        Seed for the random number generator. Set for reproducibility.

    group_key : str or None, optional
        Column name in ``adata.obs`` for a higher-level clinical grouping,
        e.g. ``"subtype"``, ``"treatment_arm"`` or ``"response"``. When
        provided, the plot is coloured by group and a group-level
        permutation test is run on the z-scores. When ``None``, each
        patient is plotted individually and no group test is performed.

    plot : bool, default True
        Whether to produce a strip-over-boxplot of permutation z-scored
        Moran's I values, optionally grouped by ``group_key``.

    figsize : tuple of float, default (5, 4)
        Figure size in inches (width, height). Ignored if ``ax`` is
        provided.

    palette : str, dict or None, optional
        Colour palette for the plot. Can be a seaborn palette name (e.g.
        ``"Set2"``), a dictionary mapping group/patient names to colours,
        or ``None`` for the default palette.

    ax : matplotlib.axes.Axes or None, optional
        Pre-existing axes to plot on. If provided, ``figsize`` is ignored
        and the plot is drawn on the given axes, which is useful for
        embedding in multi-panel figures.

    Returns
    -------
    pd.DataFrame
        One row per patient with columns:

        - ``patient``: patient/sample identifier (from ``batch_key``).
        - ``group``: clinical group label (only if ``group_key`` is set).
        - ``morans_i``: observed Moran's I on the patient's spatial graph.
        - ``expected_i``: mean Moran's I under the permutation null.
        - ``std_i``: standard deviation of Moran's I under the null.
        - ``z_score``: ``(morans_i - expected_i) / std_i``. This is the
          value used for cross-patient comparison.
        - ``p_value``: two-sided pseudo p-value from the permutation null,
          testing whether the patient's score field is significantly
          spatially autocorrelated.

        Statistical test results, when applicable, are stored as
        dictionaries in ``df.attrs``:

        - ``df.attrs["pairwise_test"]``: spot-label permutation test
          result (two-patient case), containing keys ``"morans_i_a"``,
          ``"morans_i_b"``, ``"observed_diff"``, ``"p_value"`` and
          ``"n_permutations"``.
        - ``df.attrs["group_test"]``: group-level permutation test result
          (multi-patient with ``group_key``), containing keys
          ``"mean_a"``, ``"mean_b"``, ``"observed_diff"``, ``"p_value"``
          and ``"n_permutations"``.

    Examples
    --------
    Compare Moran's I across two patients (pairwise test):

    >>> result = compare_morans_i(
    ...     adata,
    ...     score_key="enrichmap_score",
    ...     batch_key="patient_id",
    ... )
    >>> print(result[["patient", "morans_i", "z_score", "p_value"]])
       patient  morans_i   z_score  p_value
    patient_01     0.666    26.798    0.002
    patient_02    -0.017    -0.636    0.504
    >>> print(result.attrs["pairwise_test"]["p_value"])
    0.002

    Compare Moran's I across clinical groups:

    >>> result = compare_morans_i(
    ...     adata,
    ...     score_key="enrichmap_score",
    ...     batch_key="patient_id",
    ...     group_key="subtype",
    ... )
    >>> print(result.attrs["group_test"])
    {'test': 'permutation test (difference in group means)',
     'group_a': 'luminal', 'group_b': 'basal',
     'mean_a': 26.0, 'mean_b': -0.16,
     'observed_diff': 26.16, 'p_value': 0.062, ...}

    Plot on a pre-existing axes (e.g. for a multi-panel figure):

    >>> fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    >>> compare_morans_i(adata, "EMT_score", "patient_id", ax=ax1)
    >>> compare_morans_i(adata, "hypoxia_score", "patient_id", ax=ax2)

    See Also
    --------
    compare_wasserstein : Pairwise optimal transport distance between
        spatially embedded score fields.
    compare_variograms : Semivariogram-based comparison of spatial scale
        and structure.
    """
    if score_key not in adata.obs.columns:
        raise KeyError(f"'{score_key}' not found in adata.obs")
    if batch_key not in adata.obs.columns:
        raise KeyError(f"'{batch_key}' not found in adata.obs")

    rng = np.random.default_rng(random_state)
    patients = adata.obs[batch_key].unique()
    records: list[dict] = []

    for patient in tqdm(patients, desc="Moran's I per patient"):
        mask = adata.obs[batch_key] == patient
        adata_p = adata[mask].copy()

        sq.gr.spatial_neighbors(
            adata_p,
            n_neighs=n_neighbors,
            coord_type="generic",
            key_added="spatial",
        )
        W = adata_p.obsp["spatial_connectivities"]
        scores = adata_p.obs[score_key].values.astype(np.float64)
        n = len(scores)

        # Observed Moran's I
        z = scores - scores.mean()
        zWz = z @ W.dot(z)
        z2 = z @ z
        S0 = W.sum()
        observed_I = (n / S0) * (zWz / z2) if z2 > 0 else 0.0

        # Permutation null (within this patient's own graph)
        perm_Is = np.empty(n_permutations)
        for p in range(n_permutations):
            z_perm = rng.permutation(z)
            zWz_p = z_perm @ W.dot(z_perm)
            z2_p = z_perm @ z_perm
            perm_Is[p] = (n / S0) * (zWz_p / z2_p) if z2_p > 0 else 0.0

        expected_I = perm_Is.mean()
        std_I = perm_Is.std()
        z_score = (observed_I - expected_I) / std_I if std_I > 0 else 0.0
        p_value = (
            np.sum(np.abs(perm_Is - expected_I) >= np.abs(observed_I - expected_I)) + 1
        ) / (n_permutations + 1)

        row = {
            "patient": patient,
            "morans_i": observed_I,
            "expected_i": expected_I,
            "std_i": std_I,
            "z_score": z_score,
            "p_value": p_value,
        }
        if group_key is not None and group_key in adata.obs.columns:
            row["group"] = adata_p.obs[group_key].mode().iloc[0]
        records.append(row)

    result = pd.DataFrame(records)

    # Pairwise permutation test
    if len(patients) == 2:
        result.attrs["pairwise_test"] = _pairwise_permutation_morans(
            adata,
            score_key,
            batch_key,
            patients,
            n_neighbors,
            n_permutations=n_permutations,
            random_state=random_state,
        )

    # Group-level permutation test
    if group_key is not None and group_key in adata.obs.columns:
        groups = result["group"].unique()
        if len(groups) == 2:
            result.attrs["group_test"] = _group_permutation_test(
                result,
                value_col="z_score",
                group_col="group",
                n_permutations=n_permutations,
                random_state=random_state,
            )

    if plot:
        _plot_morans_i(
            result,
            group_col="group" if "group" in result.columns else None,
            figsize=figsize,
            palette=palette,
            ax=ax,
            save=save,
            save_kwargs=save_kwargs,
        )

    return result


def _pairwise_permutation_morans(
    adata,
    score_key,
    batch_key,
    patients,
    n_neighbors,
    n_permutations=999,
    random_state=0,
):
    """
    Spot-label permutation test for the difference in Moran's I between
    two patients.

    Constructs a null distribution by pooling all spots from both patients
    and randomly reassigning them to two groups of the original sizes.
    For each permutation, spatial neighbours graphs are rebuilt from
    scratch (since the spatial layout of each pseudo-patient changes) and
    Moran's I is recomputed. The test statistic is the absolute difference
    in Moran's I between the two groups.

    The null hypothesis is that the two patients' score fields are drawn
    from the same spatial process — i.e. that any difference in spatial
    autocorrelation is attributable to the random allocation of scores
    across the two tissue architectures rather than a genuine biological
    difference.

    This is an internal function called automatically by
    :func:`compare_morans_i` when exactly two patients are present.

    Parameters
    ----------
    adata : AnnData
        Full annotated data matrix (both patients).
    score_key : str
        Column in ``adata.obs`` holding the score values.
    batch_key : str
        Column in ``adata.obs`` identifying patients.
    patients : array-like
        Two-element array of patient identifiers.
    n_neighbors : int
        Number of spatial neighbours for graph construction.
    n_permutations : int, default 999
        Number of random label reassignments.
    random_state : int, default 0
        Random seed.

    Returns
    -------
    dict
        Test result with keys: ``"test"``, ``"patient_a"``,
        ``"patient_b"``, ``"morans_i_a"``, ``"morans_i_b"``,
        ``"observed_diff"``, ``"p_value"``, ``"n_permutations"``.
    """
    rng = np.random.default_rng(random_state)
    mask_a = adata.obs[batch_key] == patients[0]
    mask_b = adata.obs[batch_key] == patients[1]
    n_a, n_b = mask_a.sum(), mask_b.sum()

    def _morans_I(ad, score_col, n_neighs):
        sq.gr.spatial_neighbors(
            ad,
            n_neighs=n_neighs,
            coord_type="generic",
            key_added="spatial",
        )
        W = ad.obsp["spatial_connectivities"]
        s = ad.obs[score_col].values.astype(np.float64)
        z = s - s.mean()
        z2 = z @ z
        if z2 == 0:
            return 0.0
        return (len(s) / W.sum()) * (z @ W.dot(z) / z2)

    obs_Ia = _morans_I(adata[mask_a].copy(), score_key, n_neighbors)
    obs_Ib = _morans_I(adata[mask_b].copy(), score_key, n_neighbors)
    obs_diff = abs(obs_Ia - obs_Ib)

    all_idx = np.where(mask_a | mask_b)[0]
    null_diffs = np.empty(n_permutations)

    for i in range(n_permutations):
        perm = rng.permutation(all_idx)
        idx_a, idx_b = perm[:n_a], perm[n_a:]
        Ia = _morans_I(adata[idx_a].copy(), score_key, n_neighbors)
        Ib = _morans_I(adata[idx_b].copy(), score_key, n_neighbors)
        null_diffs[i] = abs(Ia - Ib)

    p = (np.sum(null_diffs >= obs_diff) + 1) / (n_permutations + 1)

    return {
        "test": "spot-label permutation (Moran's I difference)",
        "patient_a": patients[0],
        "patient_b": patients[1],
        "morans_i_a": obs_Ia,
        "morans_i_b": obs_Ib,
        "observed_diff": obs_diff,
        "p_value": p,
        "n_permutations": n_permutations,
    }


def _group_permutation_test(
    df,
    value_col,
    group_col,
    n_permutations=999,
    random_state=0,
):
    """
    Permutation test on the difference in group means of a per-patient
    summary statistic.

    A distribution-free alternative to a two-sample t-test. The observed
    difference in group means is compared to a null distribution built by
    randomly reassigning the group labels across patients and recomputing
    the mean difference. Works with any number of patients per group,
    including the degenerate case of n=1 per group (though power will be
    limited).

    This is a general-purpose helper used by :func:`compare_morans_i` (on
    the z-score column) and :func:`compare_variograms` (on the effective
    range column).

    Parameters
    ----------
    df : pd.DataFrame
        Per-patient summary table with at least columns ``value_col`` and
        ``group_col``.
    value_col : str
        Column containing the numeric summary statistic to compare (e.g.
        ``"z_score"`` or ``"effective_range"``).
    group_col : str
        Column containing the two-level group label.
    n_permutations : int, default 999
        Number of label permutations.
    random_state : int, default 0
        Random seed.

    Returns
    -------
    dict
        Test result with keys: ``"test"``, ``"group_a"``, ``"group_b"``,
        ``"mean_a"``, ``"mean_b"``, ``"observed_diff"``, ``"p_value"``,
        ``"n_permutations"``.
    """
    rng = np.random.default_rng(random_state)
    groups = df[group_col].unique()
    vals_a = df.loc[df[group_col] == groups[0], value_col].values
    vals_b = df.loc[df[group_col] == groups[1], value_col].values
    obs_diff = abs(vals_a.mean() - vals_b.mean())

    all_vals = np.concatenate([vals_a, vals_b])
    n_a = len(vals_a)
    null_diffs = np.empty(n_permutations)

    for i in range(n_permutations):
        perm = rng.permutation(all_vals)
        null_diffs[i] = abs(perm[:n_a].mean() - perm[n_a:].mean())

    p = (np.sum(null_diffs >= obs_diff) + 1) / (n_permutations + 1)

    return {
        "test": "permutation test (difference in group means)",
        "group_a": groups[0],
        "group_b": groups[1],
        "mean_a": vals_a.mean(),
        "mean_b": vals_b.mean(),
        "observed_diff": obs_diff,
        "p_value": p,
        "n_permutations": n_permutations,
    }


def _plot_morans_i(
    df: pd.DataFrame,
    group_col: str | None = None,
    figsize: tuple[float, float] = (5, 4),
    palette: str | dict | None = None,
    ax: plt.Axes | None = None,
    save: str | None = None,
    save_kwargs: dict | None = None,
) -> plt.Axes:
    """
    Strip-over-boxplot of permutation z-scored Moran's I values.

    Produces a combined boxplot and stripplot showing the distribution of
    per-patient z-scores. When ``group_col`` is provided, patients are
    grouped along the x-axis and coloured by group; otherwise all patients
    are plotted on a single axis. A horizontal dashed line at zero marks
    the boundary between spatially structured (above) and spatially random
    (below) score fields.

    When statistical test results are present in ``df.attrs`` (from
    :func:`compare_morans_i`), the corresponding p-value is displayed in
    the plot title.

    Parameters
    ----------
    df : pd.DataFrame
        Output of :func:`compare_morans_i`, with one row per patient and
        at minimum a ``z_score`` column. May carry
        ``df.attrs["pairwise_test"]`` and/or ``df.attrs["group_test"]``
        dictionaries with test results.

    group_col : str or None
        Column in ``df`` used for x-axis grouping and colouring. If
        ``None``, all patients appear in a single group.

    figsize : tuple of float, default (5, 4)
        Figure size in inches. Ignored when ``ax`` is provided.

    palette : str, dict or None
        Colour palette. Accepts a seaborn palette name, a dictionary
        mapping category names to colours, or ``None`` for the default.

    ax : matplotlib.axes.Axes or None
        Pre-existing axes to draw on. If ``None``, a new figure and axes
        are created.

    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if group_col is not None and group_col in df.columns:
        x_var = group_col
    else:
        df = df.copy()
        df["_all"] = ""
        x_var = "_all"

    sns.boxplot(
        data=df,
        x=x_var,
        y="z_score",
        hue=x_var,
        palette=palette,
        showfliers=False,
        width=0.5,
        boxprops=dict(alpha=0.4),
        legend=False,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x=x_var,
        y="z_score",
        hue=x_var,
        palette=palette,
        dodge=False,
        jitter=0.15,
        size=6,
        edgecolor="k",
        linewidth=0.5,
        legend=False,
        ax=ax,
    )
    ax.set_ylabel("Moran's I  (permutation z-score)")
    ax.set_xlabel("")
    ax.axhline(0, ls="--", c="grey", lw=0.8)

    if "pairwise_test" in df.attrs:
        p = df.attrs["pairwise_test"]["p_value"]
        ax.set_title(f"Spot-label permutation p = {p:.4f}", fontsize=9)
    elif "group_test" in df.attrs:
        p = df.attrs["group_test"]["p_value"]
        ax.set_title(f"Group permutation p = {p:.4f}", fontsize=9)

    sns.despine(ax=ax)
    plt.tight_layout()
    if save is not None:
        if save_kwargs is None:
            save_kwargs = {}

        save_path = Path("figures") / save
        save_path.parent.mkdir(parents=True, exist_ok=True)

        fig.savefig(save_path, **save_kwargs)
    return ax
