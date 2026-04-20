from __future__ import annotations

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from anndata import AnnData
from libpysal.weights import KNN
from esda.moran import Moran
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


def morans_correlogram(
    adata: AnnData,
    score_key: str,
    library_key: str | None = None,
    library_id: str | list[str] | None = None,
    n_neighbours: int = 6,
    save: str | None = None,
) -> None:
    """
    Plot Moran scatterplots (spatial correlograms) for one or multiple spatial libraries.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with scores in ``adata.obs[score_key]``.
    score_key : str
        Column in ``adata.obs`` containing the variable to assess spatial
        autocorrelation.
    library_key : str or None, optional
        Column in ``adata.obs`` identifying spatial libraries. If None,
        computes a single plot using the entire dataset.
    library_id : str, list of str, or None, optional
        Specific library or libraries to plot. If None, all libraries are used.
    n_neighbours : int, default 6
        Number of nearest neighbours for spatial weights.
    save : str or None, optional
        If provided, saves figure to ``figures/{save}``.
    """
    # Resolve batches
    if library_key is None:
        batches = [adata]
        titles = [None]
    else:
        all_ids = adata.obs[library_key].unique().tolist()
        if library_id is not None:
            if isinstance(library_id, str):
                library_id = [library_id]
            all_ids = [lib for lib in library_id if lib in all_ids]
        batches = [adata[adata.obs[library_key] == b].copy() for b in all_ids]
        titles = all_ids

    n = len(batches)
    ncols = min(n, 3)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(3 * ncols, 3 * nrows),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes.ravel()

    for i, (ad, title) in enumerate(zip(batches, titles)):
        coords = ad.obsm["spatial"]
        score = ad.obs[score_key].values

        mask = ~np.isnan(score)
        score = score[mask]
        coords = coords[mask]

        W = KNN.from_array(coords, k=n_neighbours)
        W.transform = "r"
        moran = Moran(score, W)
        score_lag = lag_spatial(W, score)

        # Assign quadrants
        quadrants = np.full(len(score), "LL", dtype=object)
        quadrants[(score >= 0) & (score_lag >= 0)] = "HH"
        quadrants[(score >= 0) & (score_lag < 0)] = "HL"
        quadrants[(score < 0) & (score_lag >= 0)] = "LH"
        quadrants[(score < 0) & (score_lag < 0)] = "LL"

        # Direction-aware colours and labels
        q_colours, q_labels = _quadrant_style(moran.I)

        # Rasterise for high-resolution platforms (> standard Visium)
        rasterize = len(score) > 5000

        ax = axes[i]

        # Draw discordant quadrants first, concordant on top
        if moran.I >= 0:
            draw_order = ["LH", "HL", "LL", "HH"]
        else:
            draw_order = ["HH", "LL", "LH", "HL"]

        for q in draw_order:
            q_mask = quadrants == q
            if np.any(q_mask):
                ax.scatter(
                    score[q_mask],
                    score_lag[q_mask],
                    s=10,
                    alpha=0.3,
                    color=q_colours[q],
                    edgecolors="none",
                    rasterized=rasterize,
                )

        sns.regplot(
            x=score,
            y=score_lag,
            scatter=False,
            ax=ax,
            color="black",
            line_kws={"lw": 1},
        )
        ax.axhline(0, color="grey", lw=0.8, linestyle="--")
        ax.axvline(0, color="grey", lw=0.8, linestyle="--")

        # Title
        title_str = f"Moran's I = {moran.I:.2f}, p = {moran.p_sim:.3f}"
        if title is not None:
            title_str = f"{title}\n{title_str}"
        ax.set_title(title_str, fontsize=10)
        ax.set_xlabel(score_key.replace("_score", ""), fontsize=10)
        ax.set_ylabel("Spatial lag", fontsize=10)
        ax.grid(False)
        ax.set_box_aspect(1)

        # Legend
        if moran.I >= 0:
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
        ax.legend(
            handles=legend_elements,
            fontsize=6,
            loc="upper left",
            frameon=False,
            handletextpad=0.4,
            handlelength=1,
        )

    # Remove unused axes
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    if save:
        os.makedirs("figures", exist_ok=True)
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")

    plt.show()
