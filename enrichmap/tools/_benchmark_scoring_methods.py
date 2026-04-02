from __future__ import annotations

import os
from typing import Sequence

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
    metrics: Sequence[str] = ("auroc", "auprc", "f1"),
    plot: bool = True,
    figsize: tuple[float, float] = (6, 4),
    palette: str = "muted",
    save: str | None = None,
) -> pd.DataFrame:
    """
    Benchmark multiple gene set scoring methods against a ground truth.

    Works seamlessly for a single sample or for multiple samples. When
    ``batch_key`` is provided, a leave-one-patient-out cross-validation
    scheme is used for threshold-dependent metrics, and all metrics are
    reported both per patient and as a macro average. When ``batch_key``
    is ``None``, metrics are computed on the full dataset.

    Evaluation strategy
    ~~~~~~~~~~~~~~~~~~~
    Gene set scoring methods are **unsupervised**: they compute a
    continuous enrichment score per spot from a fixed gene set without
    learning from ground truth labels. AUROC and AUPRC are therefore
    valid on the full dataset without a train/test split.

    For the F1 score, a binarisation threshold is needed. The function
    adapts its strategy based on whether multiple samples are available:

    - **Single sample** (``batch_key=None``): the F1-maximising threshold
      is found on the full dataset. Reported as **F1 (optimal)**. This is
      an upper bound and should be interpreted with caution.
    - **Multiple samples** (``batch_key`` provided): leave-one-patient-out
      CV is used. For each held-out patient, the F1-maximising threshold
      is learned on all other patients and applied to the held-out sample.
      Reported as **F1 (CV)**.

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
        Ground truth labels. If a string, interpreted as a column name
        in ``adata.obs``. If the column is non-binary, ``positive_class``
        must be set. If an ``np.ndarray``, must be binary (0/1) with
        length equal to ``adata.n_obs``.

    positive_class : str or None, optional
        When ``label_source`` is a categorical column, the category to
        treat as positive (label = 1).

    batch_key : str or None, optional
        Column in ``adata.obs`` identifying patients or slides. When
        provided, leave-one-patient-out cross-validation is used for F1
        and all metrics are reported per patient alongside the macro
        average. When ``None``, metrics are computed on the full dataset.

    metrics : sequence of str, default ``("auroc", "auprc", "f1")``
        Which metrics to compute:

        - ``"auroc"``: area under the ROC curve.
        - ``"auprc"``: area under the precision-recall curve.
        - ``"f1"``: F1 score. Automatically uses leave-one-patient-out
          CV when ``batch_key`` is provided, or the dataset-wide optimal
          threshold when it is not.

    plot : bool, default True
        Whether to produce a grouped bar chart.

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
        and ``patient``. When ``batch_key`` is provided, per-patient rows
        have the patient identifier and macro-averaged rows have
        ``patient="macro_avg"``. When ``batch_key`` is ``None``, all rows
        have ``patient="pooled"``.

    Examples
    --------
    Single-sample benchmark:

    >>> results = benchmark_scoring_methods(
    ...     adata,
    ...     score_keys=["enrichmap_score", "aucell_score"],
    ...     label_source="is_tumor",
    ... )

    Multi-sample benchmark with cross-validated F1:

    >>> results = benchmark_scoring_methods(
    ...     adata,
    ...     score_keys=["enrichmap_score", "aucell_score", "scanpy_score"],
    ...     label_source="cell_type",
    ...     positive_class="tumor",
    ...     batch_key="patient_id",
    ... )

    See Also
    --------
    compare_morans_i : Spatial autocorrelation comparison.
    compare_wasserstein : Optimal transport distance comparison.
    compare_variograms : Semivariogram-based spatial comparison.
    """
    labels = _resolve_labels(adata, label_source, positive_class)

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


def _benchmark_pooled(adata, score_keys, labels, metrics):
    """
    Pooled (single-sample) benchmark.

    AUROC and AUPRC are computed on the full dataset. F1 is computed at
    the optimal threshold (upper bound).
    """
    records: list[dict] = []

    for method in score_keys:
        if method not in adata.obs.columns:
            continue

        scores = adata.obs[method].values.astype(np.float64)
        valid = ~np.isnan(scores)
        s, y = scores[valid], labels[valid]

        if y.sum() == 0 or y.sum() == len(y):
            continue

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

        if "f1" in metrics:
            f1_opt, _ = _optimal_f1(y, s)
            records.append(
                {
                    "method": method,
                    "metric": "F1 (optimal)",
                    "value": f1_opt,
                    "patient": "pooled",
                }
            )

    return pd.DataFrame(records)


def _benchmark_cv(adata, score_keys, labels, batch_key, metrics):
    """
    Leave-one-patient-out cross-validated benchmark.

    AUROC and AUPRC are computed per patient and macro-averaged. F1 uses
    a threshold learned on the training patients and applied to the
    held-out patient.
    """
    patients = adata.obs[batch_key].unique()
    records: list[dict] = []

    for method in score_keys:
        if method not in adata.obs.columns:
            continue

        scores = adata.obs[method].values.astype(np.float64)
        per_patient = {"AUROC": [], "AUPRC": [], "F1 (CV)": []}

        for held_out in patients:
            test_mask = (adata.obs[batch_key] == held_out).values
            train_mask = ~test_mask

            s_test, y_test = scores[test_mask], labels[test_mask]
            s_train, y_train = scores[train_mask], labels[train_mask]

            if y_test.sum() == 0 or y_test.sum() == len(y_test):
                continue

            if "auroc" in metrics:
                try:
                    val = roc_auc_score(y_test, s_test)
                except ValueError:
                    val = np.nan
                per_patient["AUROC"].append(val)
                records.append(
                    {
                        "method": method,
                        "metric": "AUROC",
                        "value": val,
                        "patient": held_out,
                    }
                )

            if "auprc" in metrics:
                try:
                    val = average_precision_score(y_test, s_test)
                except ValueError:
                    val = np.nan
                per_patient["AUPRC"].append(val)
                records.append(
                    {
                        "method": method,
                        "metric": "AUPRC",
                        "value": val,
                        "patient": held_out,
                    }
                )

            if "f1" in metrics:
                if y_train.sum() == 0 or y_train.sum() == len(y_train):
                    continue
                _, thr = _optimal_f1(y_train, s_train)
                preds = (s_test > thr).astype(int)
                val = f1_score(y_test, preds, zero_division=0)
                per_patient["F1 (CV)"].append(val)
                records.append(
                    {
                        "method": method,
                        "metric": "F1 (CV)",
                        "value": val,
                        "patient": held_out,
                    }
                )

        for metric_name, values in per_patient.items():
            if values:
                records.append(
                    {
                        "method": method,
                        "metric": metric_name,
                        "value": np.nanmean(values),
                        "patient": "macro_avg",
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
                    f"positive_class to indicate which category is "
                    f"positive."
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
            "Ground truth must contain both positive and negative "
            f"cases. Found {n_pos} positives out of {len(labels)} "
            f"spots."
        )

    return labels


def _optimal_f1(y_true, scores):
    """
    Find the threshold maximising F1 by scanning the precision-recall
    curve.

    Returns (best_f1, best_threshold).
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
    """Grouped bar chart of benchmark metrics across scoring methods."""
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
