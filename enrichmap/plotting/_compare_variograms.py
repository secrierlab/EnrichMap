from __future__ import annotations
from ._compare_morans_i import _group_permutation_test

import logging
import warnings
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from tqdm import tqdm

logging.getLogger("squidpy").setLevel(logging.WARNING)


def compare_variograms(
    adata: AnnData,
    score_key: str,
    batch_key: str,
    spatial_key: str = "spatial",
    n_lags: int = 20,
    maxlag: str | float = "median",
    model: Literal["spherical", "exponential", "gaussian", "matern"] = "spherical",
    n_subsample: int | None = 3000,
    n_permutations: int = 999,
    random_state: int = 0,
    group_key: str | None = None,
    plot: bool = True,
    figsize: tuple[float, float] = (10, 4),
    palette: str | dict | None = None,
) -> pd.DataFrame:
    """
    Fit empirical semivariograms to EnrichMap scores per patient and extract
    structural parameters for cross-patient comparison of spatial organisation.

    A semivariogram describes how the dissimilarity between score values
    changes as a function of the spatial distance (lag) between spots. Three
    parameters are extracted from a fitted theoretical model:

    - **Effective range**: the lag distance at which the variogram plateaus,
      i.e. the spatial scale of autocorrelation. A short range means the
      score field is composed of many small patches; a long range means
      large, coherent spatial domains. Because coordinates are min-max
      normalised per patient to [0, 1], the range is expressed as a
      fraction of tissue extent and is directly comparable across slides
      of different physical sizes.
    - **Sill**: the plateau semivariance, representing the total spatial
      variance of the score field.
    - **Nugget**: the semivariance at lag zero, capturing measurement noise
      or fine-scale variation below the resolution of the spatial graph.
      A high nugget-to-sill ratio indicates that most of the variance is
      spatially unstructured.

    Two patients can share identical score distributions yet have completely
    different variograms, making this a powerful complement to marginal
    summaries like violin plots.

    Statistical testing
    ~~~~~~~~~~~~~~~~~~~
    When exactly **two patients** are present, a spot-label permutation test
    is run automatically on the difference in effective range. All spots from
    both patients are pooled and randomly reassigned to two groups of the
    original sizes; variograms are refitted on each permutation to build a
    null distribution. This tests whether the observed difference in spatial
    scale is larger than expected if score values were randomly distributed
    across the two tissue architectures.

    When **multiple patients** are present and ``group_key`` is provided, a
    permutation test is run on the difference in group-mean effective range,
    analogous to a two-sample t-test but without distributional assumptions.

    All test results are stored in ``df.attrs`` and printed in the plot title
    when ``plot=True``.

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
        slides, e.g. ``"patient_id"`` or ``"library_id"``. Each unique value
        is treated as a separate sample with its own spatial graph and
        variogram.

    spatial_key : str, default ``"spatial"``
        Key in ``adata.obsm`` containing the 2D spatial coordinates as an
        (n_spots, 2) array.

    n_lags : int, default 20
        Number of evenly spaced lag bins for the empirical variogram. More
        bins give a smoother curve but require more spot pairs per bin.
        Values between 15 and 30 work well for Visium-scale data.

    maxlag : str or float, default ``"median"``
        Maximum lag distance considered. ``"median"`` (recommended) uses the
        median pairwise distance within each sample, which is a robust
        default that avoids noisy estimates at large lags where few spot
        pairs exist. Can also be a float in normalised coordinate units
        (e.g. ``0.5`` to consider lags up to half the tissue extent).

    model : ``{"spherical", "exponential", "gaussian", "matern"}``, default ``"spherical"``
        Theoretical variogram model fitted to the empirical semivariance.
        ``"spherical"`` is the most commonly used and has a clear transition
        from spatially structured to unstructured variance.
        ``"exponential"`` approaches the sill asymptotically (no sharp
        plateau). ``"gaussian"`` implies very smooth spatial fields.
        ``"matern"`` is the most flexible but may overfit with few lags.

    n_subsample : int or None, default 3000
        If set, subsample each patient to at most this many spots before
        fitting. Variogram estimation is O(n²) in the number of spots, so
        subsampling is recommended for slides with more than ~4000 spots.
        Set to ``None`` to use all spots.

    n_permutations : int, default 999
        Number of permutations for statistical testing. Higher values give
        more precise p-values. For exploratory analysis 499 is sufficient;
        for publication-ready results use 999 or higher.

    random_state : int, default 0
        Seed for the random number generator, used for subsampling and
        permutation tests. Set for reproducibility.

    group_key : str or None, optional
        Column name in ``adata.obs`` for a higher-level clinical grouping,
        e.g. ``"subtype"``, ``"treatment_arm"`` or ``"response"``. When
        provided, the plot is coloured by group and a group-level
        permutation test is run on the effective range. When ``None``,
        each patient is plotted individually.

    plot : bool, default True
        Whether to produce a three-panel figure showing overlaid variogram
        curves (left), effective range comparison (centre) and sill
        comparison (right).

    figsize : tuple of float, default (10, 4)
        Figure size in inches (width, height).

    palette : str, dict or None, optional
        Colour palette for the plot. Can be a seaborn palette name (e.g.
        ``"Set2"``), a dictionary mapping group/patient names to colours,
        or ``None`` for the default palette.

    Returns
    -------
    pd.DataFrame
        One row per patient with columns:

        - ``patient``: patient/sample identifier (from ``batch_key``).
        - ``group``: clinical group label (only if ``group_key`` is set).
        - ``effective_range``: spatial autocorrelation scale (normalised).
        - ``sill``: total spatial variance (plateau semivariance).
        - ``nugget``: fine-scale / measurement noise variance.
        - ``nugget_sill_ratio``: proportion of variance that is spatially
          unstructured (``nugget / sill``).

        Statistical test results, when applicable, are stored as
        dictionaries in ``df.attrs``:

        - ``df.attrs["pairwise_test"]``: spot-label permutation test result
          (two-patient case).
        - ``df.attrs["group_test"]``: group-level permutation test result
          (multi-patient with ``group_key``).

    Examples
    --------
    Compare variograms across two patients (pairwise test):

    >>> result = compare_variograms(
    ...     adata,
    ...     score_key="enrichmap_score",
    ...     batch_key="patient_id",
    ...     n_lags=15,
    ... )
    >>> print(result)
       patient  effective_range   sill  nugget  nugget_sill_ratio
    patient_01            0.532  1.170       0                0.0
    patient_02            0.069  0.972       0                0.0
    >>> print(result.attrs["pairwise_test"]["p_value"])
    0.005

    Compare variograms across clinical groups:

    >>> result = compare_variograms(
    ...     adata,
    ...     score_key="enrichmap_score",
    ...     batch_key="patient_id",
    ...     group_key="subtype",
    ... )
    >>> print(result.attrs["group_test"])
    {'test': 'permutation test (difference in group means)',
     'group_a': 'luminal', 'group_b': 'basal',
     'mean_a': 0.532, 'mean_b': 0.063,
     'observed_diff': 0.469, 'p_value': 0.088, ...}

    See Also
    --------
    compare_morans_i : Spatial autocorrelation comparison via Moran's I.
    compare_wasserstein : Pairwise optimal transport distance between
        spatially embedded score fields.
    """
    try:
        from skgstat import Variogram
    except ImportError:
        raise ImportError(
            "scikit-gstat is required. Install with: pip install scikit-gstat"
        )

    if score_key not in adata.obs.columns:
        raise KeyError(f"'{score_key}' not found in adata.obs")

    rng = np.random.default_rng(random_state)
    patients = adata.obs[batch_key].unique()
    records: list[dict] = []
    variogram_curves: list[dict] = []

    def _fit_variogram(coords, scores):
        """Normalise, subsample, fit, extract parameters."""
        coords = coords.copy().astype(np.float64)
        scores = scores.copy().astype(np.float64)

        for d in range(coords.shape[1]):
            dmin, dmax = coords[:, d].min(), coords[:, d].max()
            rng_d = dmax - dmin
            coords[:, d] = (coords[:, d] - dmin) / rng_d if rng_d > 0 else 0.0

        if n_subsample is not None and len(scores) > n_subsample:
            idx = rng.choice(len(scores), n_subsample, replace=False)
            coords, scores = coords[idx], scores[idx]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            V = Variogram(
                coords,
                scores,
                n_lags=n_lags,
                maxlag=maxlag,
                model=model,
                normalize=False,
            )
        params = V.parameters
        eff_range = params[0]
        sill = params[1]
        nugget = params[2] if len(params) > 2 else 0.0
        return V, eff_range, sill, nugget

    for patient in tqdm(patients, desc="Fitting variograms"):
        mask = adata.obs[batch_key] == patient
        coords = adata.obsm[spatial_key][mask]
        scores = adata.obs.loc[mask, score_key].values

        V, eff_range, sill, nugget = _fit_variogram(coords, scores)

        row = {
            "patient": patient,
            "effective_range": eff_range,
            "sill": sill,
            "nugget": nugget,
            "nugget_sill_ratio": nugget / sill if sill > 0 else np.nan,
        }
        if group_key is not None and group_key in adata.obs.columns:
            row["group"] = adata.obs.loc[mask, group_key].mode().iloc[0]
        records.append(row)

        variogram_curves.append(
            {
                "patient": patient,
                "group": row.get("group", None),
                "lags": V.bins,
                "experimental": V.experimental,
                "fitted": V.transform(V.bins),
            }
        )

    result = pd.DataFrame(records)

    # Pairwise permutation test (2 patients)
    if len(patients) == 2:
        obs_diff = abs(
            result["effective_range"].iloc[0] - result["effective_range"].iloc[1]
        )
        mask_a = (adata.obs[batch_key] == patients[0]).values
        mask_b = (adata.obs[batch_key] == patients[1]).values
        n_a = mask_a.sum()
        pooled_idx = np.where(mask_a | mask_b)[0]

        null_diffs = np.empty(n_permutations)
        for i in range(n_permutations):
            perm = rng.permutation(pooled_idx)
            idx_a, idx_b = perm[:n_a], perm[n_a:]

            c_a = adata.obsm[spatial_key][idx_a]
            s_a = adata.obs[score_key].values[idx_a]
            c_b = adata.obsm[spatial_key][idx_b]
            s_b = adata.obs[score_key].values[idx_b]

            _, r_a, _, _ = _fit_variogram(c_a, s_a)
            _, r_b, _, _ = _fit_variogram(c_b, s_b)
            null_diffs[i] = abs(r_a - r_b)

        p = (np.sum(null_diffs >= obs_diff) + 1) / (n_permutations + 1)
        result.attrs["pairwise_test"] = {
            "test": "spot-label permutation (effective range difference)",
            "patient_a": patients[0],
            "patient_b": patients[1],
            "range_a": result["effective_range"].iloc[0],
            "range_b": result["effective_range"].iloc[1],
            "observed_diff": obs_diff,
            "p_value": p,
            "n_permutations": n_permutations,
        }

    # Group-level permutation test
    if group_key is not None and "group" in result.columns:
        groups = result["group"].unique()
        if len(groups) == 2:
            result.attrs["group_test"] = _group_permutation_test(
                result,
                value_col="effective_range",
                group_col="group",
                n_permutations=n_permutations,
                random_state=random_state,
            )

    if plot:
        _plot_variograms(
            result,
            variogram_curves,
            group_col="group" if "group" in result.columns else None,
            figsize=figsize,
            palette=palette,
        )

    return result


def _plot_variograms(
    df: pd.DataFrame,
    curves: list[dict],
    group_col: str | None = None,
    figsize: tuple[float, float] = (10, 4),
    palette: str | dict | None = None,
) -> tuple[plt.Figure, np.ndarray]:
    """
    Three-panel figure summarising variogram comparison results.

    Generates a figure with three horizontally arranged panels:

    1. **Variogram curves** (left): overlaid empirical (dots) and fitted
       (lines) semivariograms for all patients, coloured by group or
       patient identity. A gradually rising curve indicates spatially
       structured scores; a flat curve indicates spatially random scores.
    2. **Effective range** (centre): strip-over-boxplot comparing the
       effective range parameter across groups or patients. Higher values
       correspond to larger coherent spatial domains.
    3. **Sill** (right): strip-over-boxplot comparing the sill parameter.
       The sill reflects total spatial variance and is expected to be
       similar across patients if score distributions were normalised.

    When statistical test results are present in ``df.attrs`` (from
    ``compare_variograms``), the corresponding p-value is displayed in
    the figure title.

    Parameters
    ----------
    df : pd.DataFrame
        Output of :func:`compare_variograms`, with one row per patient
        and columns ``patient``, ``effective_range``, ``sill``, and
        optionally ``group``. May carry ``df.attrs["pairwise_test"]``
        and/or ``df.attrs["group_test"]`` dictionaries with test results.

    curves : list of dict
        Variogram curve data for each patient. Each dictionary contains:

        - ``"patient"``: patient identifier.
        - ``"group"``: group label or ``None``.
        - ``"lags"``: array of lag bin centres (normalised coordinates).
        - ``"experimental"``: array of empirical semivariance values.
        - ``"fitted"``: array of fitted model semivariance values.

    group_col : str or None
        Column in ``df`` used for colouring. If ``None``, patients are
        coloured individually; if set (e.g. ``"group"``), patients sharing
        a group label receive the same colour.

    figsize : tuple of float, default (10, 4)
        Figure size in inches (width, height).

    palette : str, dict or None
        Colour palette. Accepts a seaborn palette name, a dictionary
        mapping category names to colours, or ``None`` for the default.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.

    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of the three axes objects, in order: variogram curves,
        effective range comparison, sill comparison.
    """
    if group_col is not None and group_col in df.columns:
        colour_col = group_col
    else:
        colour_col = "patient"

    categories = df[colour_col].unique()
    if palette is None:
        pal = dict(zip(categories, sns.color_palette("Set2", len(categories))))
    elif isinstance(palette, str):
        pal = dict(zip(categories, sns.color_palette(palette, len(categories))))
    else:
        pal = palette

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # Panel 1: variogram curves
    ax0 = axes[0]
    for curve in curves:
        key = curve.get("group") or curve["patient"]
        colour = pal.get(key, "grey")
        ax0.plot(
            curve["lags"],
            curve["experimental"],
            "o",
            ms=2.5,
            alpha=0.4,
            color=colour,
        )
        ax0.plot(
            curve["lags"],
            curve["fitted"],
            "-",
            alpha=0.6,
            color=colour,
            lw=1.2,
        )
        ax0.set_box_aspect(1)
        ax0.grid(False)

    for cat, col in pal.items():
        ax0.plot([], [], "-", color=col, lw=2, label=cat)
    ax0.legend(fontsize=8, frameon=False)
    ax0.set_xlabel("Lag (normalised)")
    ax0.set_ylabel("Semivariance")
    ax0.set_title("Empirical variograms", fontsize=10)
    sns.despine(ax=ax0)

    # Panels 2-3: parameter comparisons
    param_pairs = [
        ("effective_range", "Effective range\n(normalised)"),
        ("sill", "Sill\n(spatial variance)"),
    ]
    for ax, (col, ylabel) in zip(axes[1:], param_pairs):
        sns.boxplot(
            data=df,
            x=colour_col,
            y=col,
            hue=colour_col,
            palette=pal,
            showfliers=False,
            width=0.5,
            boxprops=dict(alpha=0.4),
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x=colour_col,
            y=col,
            hue=colour_col,
            palette=pal,
            dodge=False,
            jitter=0.15,
            size=5,
            edgecolor="k",
            linewidth=0.5,
            legend=False,
            ax=ax,
        )
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        ax.grid(False)
        ax.set_box_aspect(1)
        sns.despine(ax=ax)

    title_parts = []
    if "pairwise_test" in df.attrs:
        p = df.attrs["pairwise_test"]["p_value"]
        title_parts.append(f"Spot-label permutation p = {p:.4f}")
    if "group_test" in df.attrs:
        p = df.attrs["group_test"]["p_value"]
        title_parts.append(f"Group permutation p = {p:.4f}")
    if title_parts:
        fig.suptitle("  |  ".join(title_parts), fontsize=9, y=1.02)

    plt.tight_layout()
    return fig, axes
