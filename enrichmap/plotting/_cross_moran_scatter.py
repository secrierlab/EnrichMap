from __future__ import annotations

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from anndata import AnnData
from scipy.stats import pearsonr
from libpysal.weights import KNN
from libpysal.weights.spatial_lag import lag_spatial

plt.rcParams["pdf.fonttype"] = "truetype"


# Quadrant styling based on correlation direction
def _quadrant_style(r: float) -> tuple[dict, dict]:
    """Return colours and labels appropriate for the correlation sign."""
    if r >= 0:
        colours = {
            "HH": "#c62828",  # red — concordant
            "LL": "#1565c0",  # blue — concordant
            "LH": "#b0bec5",  # grey — discordant
            "HL": "#b0bec5",
        }
        labels = {
            "HH": "Co-enriched",
            "LL": "Co-depleted",
            "LH": "Discordant",
            "HL": "Discordant",
        }
    else:
        colours = {
            "LH": "#7b1fa2",  # purple — anti-correlated
            "HL": "#7b1fa2",
            "HH": "#b0bec5",  # grey — discordant
            "LL": "#b0bec5",
        }
        labels = {
            "LH": "Anti-correlated",
            "HL": "Anti-correlated",
            "HH": "Discordant",
            "LL": "Discordant",
        }
    return colours, labels


def cross_moran_scatter(
    adata: AnnData,
    score_x: str,
    score_y: str,
    library_key: str | None = None,
    library_id: str | list[str] | None = None,
    n_neighbours: int = 6,
    save: str | None = None,
):
    """
    Plot cross-Moran scatterplots: score_x vs spatial lag of score_y, per library or globally.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in `adata.obsm["spatial"]`.
    score_x : str
        Column in `adata.obs` for X-axis variable.
    score_y : str
        Column in `adata.obs` for which the spatial lag is computed (Y-axis).
    library_key : str or None
        Key in `adata.obs` indicating batch/sample/library identifiers.
        If None, computes a single plot using the entire dataset.
    library_id : str, list of str, or None
        Specific library or libraries to plot. If None, all libraries are used.
    n_neighbours : int
        Number of nearest neighbours for spatial weights.
    save : str or None
        If provided, saves the figure to the given file name.
    """
    if library_key is None:
        fig, ax = plt.subplots(figsize=(3, 3))
        batches = [adata]
        titles = ["Cross-Moran scatter plot"]
    else:
        all_ids = adata.obs[library_key].unique().tolist()
        selected_ids = all_ids

        if library_id is not None:
            if isinstance(library_id, str):
                library_id = [library_id]
            selected_ids = [lib for lib in library_id if lib in all_ids]

        batches = [adata[adata.obs[library_key] == b].copy() for b in selected_ids]
        titles = selected_ids
        n = len(batches)
        ncols = int(np.ceil(np.sqrt(n)))
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(3 * ncols, 3 * nrows), constrained_layout=True
        )
        axes = np.ravel(axes)

    for i, (ad, title) in enumerate(zip(batches, titles)):
        coords = ad.obsm["spatial"]
        x = ad.obs[score_x].values
        y = ad.obs[score_y].values

        W = KNN.from_array(coords, k=n_neighbours)
        W.transform = "r"
        y_lag = lag_spatial(W, y)

        mask = ~np.isnan(x) & ~np.isnan(y_lag)
        if np.sum(mask) > 1:
            r, p = pearsonr(x[mask], y_lag[mask])
        else:
            r, p = np.nan, np.nan

        # Assign quadrant labels
        quadrants = np.full(len(x), "LL", dtype=object)
        quadrants[(x >= 0) & (y_lag >= 0)] = "HH"
        quadrants[(x >= 0) & (y_lag < 0)] = "HL"
        quadrants[(x < 0) & (y_lag >= 0)] = "LH"
        quadrants[(x < 0) & (y_lag < 0)] = "LL"

        # Direction-aware colours and labels
        q_colours, q_labels = _quadrant_style(r)

        current_ax = ax if library_key is None else axes[i]

        # Rasterise for high-resolution platforms (> standard Visium)
        rasterize = len(x) > 5000

        # Draw discordant quadrants first, concordant on top
        if r >= 0:
            draw_order = ["LH", "HL", "LL", "HH"]
        else:
            draw_order = ["HH", "LL", "LH", "HL"]

        for q in draw_order:
            q_mask = quadrants == q
            if np.any(q_mask):
                current_ax.scatter(
                    x[q_mask],
                    y_lag[q_mask],
                    s=10,
                    alpha=0.3,
                    color=q_colours[q],
                    edgecolors="none",
                    rasterized=rasterize,
                )

        sns.regplot(
            x=x,
            y=y_lag,
            scatter=False,
            ax=current_ax,
            color="black",
            line_kws={"lw": 1},
        )
        current_ax.axhline(0, color="grey", lw=0.8, linestyle="--")
        current_ax.axvline(0, color="grey", lw=0.8, linestyle="--")
        current_ax.set_title(f"{title}\nr = {r:.2f}, p = {p:.2g}", fontsize=10)
        current_ax.set_xlabel(score_x.replace("_score", ""), fontsize=10)
        current_ax.set_ylabel(
            f"Spatial lag of {score_y.replace('_score', '')}", fontsize=10
        )
        current_ax.grid(False)
        current_ax.set_box_aspect(1)

        # Legend — show concordant pair + discordant
        if r >= 0:
            legend_elements = [
                Patch(facecolor=q_colours["HH"], label=q_labels["HH"]),
                Patch(facecolor=q_colours["LL"], label=q_labels["LL"]),
                Patch(facecolor=q_colours["LH"], label=q_labels["LH"]),
            ]
        else:
            legend_elements = [
                Patch(facecolor=q_colours["LH"], label=q_labels["LH"]),
                Patch(facecolor=q_colours["HH"], label=q_labels["HH"]),
            ]
        current_ax.legend(
            handles=legend_elements,
            fontsize=6,
            loc="upper left",
            frameon=False,
            handletextpad=0.4,
            handlelength=1,
        )

    if library_key is not None:
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")

    if save:
        os.makedirs("figures", exist_ok=True)
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")

    plt.show()
