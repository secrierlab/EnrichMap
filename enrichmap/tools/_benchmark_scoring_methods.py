from __future__ import annotations

import os
from typing import Literal, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from sklearn.metrics import (
    average_precision_score,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
)


def benchmark_scoring_methods(
    adata: AnnData,
    score_keys: Sequence[str],
    label_source: str | np.ndarray,
    positive_class: str | None = None,
    batch_key: str | None = None,
    metrics: Sequence[str] = ("auroc", "auprc", "f1_cv"),
    plot: bool = True,
    figsize: tuple[float, float] = (6, 4),
    palette: str = "muted",
    save: str | None = None,
) -> pd.DataFrame:
    """
    Benchmark multiple gene set scoring methods against a ground truth
    using threshold-free and cross-validated classification metrics.

    Given binary ground truth labels and continuous scores from different
    methods (e.g. EnrichMap, AUCell, scanpy ``score_genes``, ssGSEA), this
    function evaluates the discriminative performance of each method using
    a leave-one-patient-out cross-validation scheme that prevents any
    information leakage between threshold selection and evaluation.

    Evaluation strategy
    ~~~~~~~~~~~~~~~~~~~
    Gene set scoring methods are **unsupervised**: they take a fixed gene
    set and compute a score per spot without learning from the ground truth
    labels. There is therefore no "training" in the supervised machine
    learning sense, and no risk of label leakage in the scores themselves.
    However, converting continuous scores into binary predictions requires
    choosing a threshold, and selecting that threshold on the same data
    used for evaluation inflates performance.

    This function addresses the concern as follows:

    - **AUROC** and **AUPRC** are threshold-free ranking metrics that
      require no binarisation. They are computed per patient (when
      ``batch_key`` is provided) and macro-averaged across patients.
    - **F1 (cross-validated)** uses leave-one-patient-out CV: for each
      held-out patient, the F1-maximising threshold is learned on all
      remaining patients and applied to the held-out patient. The
      per-patient F1 values are then macro-averaged. This ensures the
      threshold is never optimised on the data it is evaluated on.

    When ``batch_key`` is ``None``, metrics are computed on the full
    dataset without cross-validation (suitable for exploratory analysis
    but not for reporting in a manuscript).

    Ground truth labels
    ~~~~~~~~~~~~~~~~~~~
    Labels can be provided in three ways:

    1. **Pre-computed array**: pass a binary ``np.ndarray`` directly.
    2. **Binary column**: pass the name of a column in ``adata.obs`` that
       already contains 0/1 values.
    3. **Categorical column + positive class**: pass a column name and set
       ``positive_class`` to the category that should be labelled 1 (e.g.
       ``label_source="cell_type", positive_class="tumor"``).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with scoring results in ``adata.obs``.

    score_keys : sequence of str
        Column names in ``adata.obs`` for each scoring method's output,
        e.g. ``["enrichmap_score", "aucell_score", "scanpy_score"]``.

    label_source : str or np.ndarray
        Ground truth labels. If a string, interpreted as a column name in
        ``adata.obs``. If the column is non-binary, ``positive_class``
        must be set. If an ``np.ndarray``, must be binary (0/1) with
        length equal to ``adata.n_obs``.

    positive_class : str or None, optional
        When ``label_source`` is a categorical column, the category to
        treat as positive (label = 1).

    batch_key : str or None, optional
        Column in ``adata.obs`` identifying patients or slides. When
        provided, leave-one-patient-out cross-validation is used for
        threshold-dependent metrics and all metrics are reported per
        patient alongside the macro average. **Required for F1 (CV).**

    metrics : sequence of str, default ``("auroc", "auprc", "f1_cv")``
        Which metrics to compute:

        - ``"auroc"``: area under the ROC curve.
        - ``"auprc"``: area under the precision-recall curve.
        - ``"f1_cv"``: cross-validated F1 (leave-one-patient-out).
          Requires ``batch_key``.

    plot : bool, default True
        Whether to produce a grouped bar chart of macro-averaged results.

    figsize : tuple of float, default (6, 4)
        Figure size in inches.

    palette : str, default ``"muted"``
        Seaborn colour palette.

    save : str or None, optional
        Path to save the figure.

    Returns
    -------
    pd.DataFrame
        Tidy dataframe with columns ``method``, ``metric``, ``value``
        and ``patient``. Macro-averaged rows have ``patient="macro_avg"``.

    Examples
    --------
    Cross-validated benchmark across patients:

    >>> results = benchmark_scoring_methods(
    ...     adata,
    ...     score_keys=["enrichmap_score", "aucell_score", "scanpy_score"],
    ...     label_source="cell_type",
    ...     positive_class="tumor",
    ...     batch_key="patient_id",
    ... )

    Quick exploratory benchmark (no CV):

    >>> results = benchmark_scoring_methods(
    ...     adata,
    ...     score_keys=["enrichmap_score", "aucell_score"],
    ...     label_source="is_tumor",
    ...     metrics=("auroc", "auprc"),
    ... )
    """
    # Resolve ground truth labels

    labels = _resolve_labels(adata, label_source, positive_class)

    # Validate

    if "f1_cv" in metrics and batch_key is None:
        raise ValueError(
            "batch_key is required for cross-validated F1 ('f1_cv'). "
            "Provide a patient/sample identifier column, or remove "
            "'f1_cv' from metrics for a non-CV evaluation."
        )

    # Route to CV or pooled computation

    if batch_key is not None:
        result = _benchmark_cv(adata, score_keys, labels, batch_key, metrics)
    else:
        result = _benchmark_pooled(adata, score_keys, labels, metrics)

    if plot and not result.empty:
        plot_df = result[result["patient"] == "macro_avg"].copy()
        if plot_df.empty:
            plot_df = result.copy()
        _plot_benchmark(plot_df, figsize=figsize, palette=palette, save=save)

    return result


# Core computation


def _benchmark_cv(adata, score_keys, labels, batch_key, metrics):
    """
    Leave-one-patient-out cross-validated benchmark.

    For threshold-free metrics (AUROC, AUPRC), per-patient values are
    computed and macro-averaged.

    For F1 (CV), the optimal threshold is learned on all patients except
    the held-out one, then applied to the held-out patient. Per-patient
    F1 values are macro-averaged.
    """
    patients = adata.obs[batch_key].unique()
    records: list[dict] = []

    for method in score_keys:
        if method not in adata.obs.columns:
            continue

        scores = adata.obs[method].values.astype(np.float64)

        patient_aurocs, patient_auprcs, patient_f1s = [], [], []

        for held_out in patients:
            test_mask = (adata.obs[batch_key] == held_out).values
            train_mask = ~test_mask

            s_test = scores[test_mask]
            y_test = labels[test_mask]
            s_train = scores[train_mask]
            y_train = labels[train_mask]

            # Skip patients with no positive or no negative labels
            if y_test.sum() == 0 or y_test.sum() == len(y_test):
                continue

            # AUROC per patient
            if "auroc" in metrics:
                try:
                    val = roc_auc_score(y_test, s_test)
                except ValueError:
                    val = np.nan
                patient_aurocs.append(val)
                records.append(
                    {
                        "method": method,
                        "metric": "AUROC",
                        "value": val,
                        "patient": held_out,
                    }
                )

            # AUPRC per patient
            if "auprc" in metrics:
                try:
                    val = average_precision_score(y_test, s_test)
                except ValueError:
                    val = np.nan
                patient_auprcs.append(val)
                records.append(
                    {
                        "method": method,
                        "metric": "AUPRC",
                        "value": val,
                        "patient": held_out,
                    }
                )

            # F1 (CV): threshold learned on training patients
            if "f1_cv" in metrics:
                _, thr = _optimal_f1(y_train, s_train)
                preds = (s_test > thr).astype(int)
                val = f1_score(y_test, preds, zero_division=0)
                patient_f1s.append(val)
                records.append(
                    {
                        "method": method,
                        "metric": "F1 (CV)",
                        "value": val,
                        "patient": held_out,
                    }
                )

        # Macro averages
        if "auroc" in metrics and patient_aurocs:
            records.append(
                {
                    "method": method,
                    "metric": "AUROC",
                    "value": np.nanmean(patient_aurocs),
                    "patient": "macro_avg",
                }
            )
        if "auprc" in metrics and patient_auprcs:
            records.append(
                {
                    "method": method,
                    "metric": "AUPRC",
                    "value": np.nanmean(patient_auprcs),
                    "patient": "macro_avg",
                }
            )
        if "f1_cv" in metrics and patient_f1s:
            records.append(
                {
                    "method": method,
                    "metric": "F1 (CV)",
                    "value": np.mean(patient_f1s),
                    "patient": "macro_avg",
                }
            )

    return pd.DataFrame(records)


def _benchmark_pooled(adata, score_keys, labels, metrics):
    """Pooled (non-CV) benchmark for exploratory use."""
    records: list[dict] = []

    for method in score_keys:
        if method not in adata.obs.columns:
            continue

        scores = adata.obs[method].values.astype(np.float64)
        valid = ~np.isnan(scores)
        s, y = scores[valid], labels[valid]

        if "auroc" in metrics:
            try:
                val = roc_auc_score(y, s)
            except ValueError:
                val = np.nan
            records.append(
                {
                    "method": method,
                    "metric": "AUROC",
                    "value": val,
                    "patient": "pooled",
                }
            )

        if "auprc" in metrics:
            try:
                val = average_precision_score(y, s)
            except ValueError:
                val = np.nan
            records.append(
                {
                    "method": method,
                    "metric": "AUPRC",
                    "value": val,
                    "patient": "pooled",
                }
            )

    return pd.DataFrame(records)


# Helpers


def _resolve_labels(adata, label_source, positive_class):
    """Resolve label_source to a binary numpy array."""
    if isinstance(label_source, np.ndarray):
        labels = label_source.astype(int)
    elif isinstance(label_source, str):
        if label_source not in adata.obs.columns:
            raise KeyError(f"'{label_source}' not found in adata.obs")
        col = adata.obs[label_source]
        if positive_class is not None:
            labels = (col == positive_class).astype(int).values
        else:
            vals = col.unique()
            if set(vals).issubset({0, 1, True, False}):
                labels = col.astype(int).values
            else:
                raise ValueError(
                    f"Column '{label_source}' is not binary. Set "
                    f"positive_class to indicate which category is positive."
                )
    else:
        raise TypeError(
            "label_source must be a column name (str) or binary np.ndarray."
        )

    if len(labels) != adata.n_obs:
        raise ValueError(
            f"Label length ({len(labels)}) != adata.n_obs ({adata.n_obs})."
        )

    n_pos = labels.sum()
    if n_pos == 0 or n_pos == len(labels):
        raise ValueError(
            "Ground truth must contain both positive and negative cases. "
            f"Found {n_pos} positives out of {len(labels)} spots."
        )

    return labels


def _optimal_f1(y_true, scores):
    """
    Find the threshold maximising F1 by scanning the precision-recall
    curve.

    Returns
    -------
    best_f1 : float
    best_threshold : float
    """
    precision, recall, thresholds = precision_recall_curve(y_true, scores)
    precision, recall = precision[:-1], recall[:-1]

    with np.errstate(divide="ignore", invalid="ignore"):
        f1 = 2 * (precision * recall) / (precision + recall)
    f1 = np.nan_to_num(f1, nan=0.0)

    best_idx = np.argmax(f1)
    return float(f1[best_idx]), float(thresholds[best_idx])


def _plot_benchmark(
    df: pd.DataFrame,
    figsize: tuple[float, float] = (6, 4),
    palette: str = "muted",
    save: str | None = None,
) -> plt.Axes:
    """Grouped bar chart of macro-averaged benchmark metrics."""
    df = df.copy()
    df["method_label"] = df["method"].str.replace("_score", "", regex=False)

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=df,
        x="method_label",
        y="value",
        hue="metric",
        palette=palette,
        edgecolor="white",
        linewidth=0.5,
        ax=ax,
    )

    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score", fontsize=9)
    ax.set_xlabel("")
    ax.set_title("Scoring method benchmark", fontsize=10)
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.legend(title="", fontsize=7, frameon=False, loc="upper right")
    ax.grid(False)
    sns.despine(ax=ax)
    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save) or "figures", exist_ok=True)
        path = os.path.join("figures", save) if not os.path.dirname(save) else save
        plt.savefig(path, dpi=300, bbox_inches="tight")

    plt.show()
    return ax
