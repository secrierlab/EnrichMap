from __future__ import annotations

import os
import squidpy as sq
import matplotlib.pyplot as plt
import numpy as np

from anndata import AnnData
from libpysal.weights import KNN
from esda.moran import Moran
from splot.esda import moran_scatterplot

plt.rcParams["pdf.fonttype"] = "truetype"

def morans_correlogram(
    adata: AnnData,
    score_key: str,
    library_key: str | None = None,
    library_id: str | None = None,
    save: str | None = None,
    multipanel: bool = False
) -> None:
    """
    Plot Moran scatterplots (spatial correlograms) for one or multiple spatial libraries.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with scores in `adata.obs[score_key]`.

    score_key : str
        Column in `adata.obs` containing the variable to assess spatial autocorrelation.

    library_key : str or None, optional
        Column in `adata.obs` identifying spatial libraries. Required for `multipanel=True`.

    library_id : str or None, optional
        Specific library to plot. Only used if `multipanel=False`.

    save : str or None, optional
        If provided, saves figure to `figures/{save}`.

    multipanel : bool, default False
        If True, plots Moran scatterplots for all libraries in a grid layout.
    """
    scatter_kwds = {
        "s": 10,
        "alpha": 0.3,
        "edgecolor": "k",
        "linewidth": 0,
    }
    fitline_kwds = {
        "color": "darkred",
        "linestyle": "-",
        "linewidth": 1,
        "label": "Fit line"
    }

    if multipanel:
        if library_key is None:
            raise ValueError("`library_key` must be provided for multipanel plotting.")

        libraries = adata.obs[library_key].unique().tolist()
        n_panels = len(libraries)
        ncols = min(n_panels, 3)
        nrows = int(np.ceil(n_panels / ncols))

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(3 * ncols, 3 * nrows),
            constrained_layout=True,
            sharex=True,
            sharey=True
        )
        axes = axes.flatten()

        for i, lib in enumerate(libraries):
            adata_subset = adata[adata.obs[library_key] == lib].copy()

            if "spatial_neighbors" not in adata_subset.uns:
                sq.gr.spatial_neighbors(adata_subset, n_neighs=6, coord_type="generic", key_added="spatial")

            score = adata_subset.obs[score_key].values
            spatial = adata_subset.obsm["spatial"]
            mask = ~np.isnan(score)
            score = score[mask]
            spatial = spatial[mask]

            W = KNN.from_array(spatial)
            W.transform = "r"
            moran = Moran(score, W)

            moran_scatterplot(moran, p=0.05, ax=axes[i], aspect_equal=True,
                              scatter_kwds=scatter_kwds, fitline_kwds=fitline_kwds)
            axes[i].set_title(f"Slide {lib}\nMoran’s I = {moran.I:.2f}, p = {moran.p_sim:.3f}", fontsize=10)
            axes[i].set_xlabel("Score")
            axes[i].set_ylabel("Spatial lag")
            axes[i].grid(False)

        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        if save:
            os.makedirs("figures", exist_ok=True)
            if not os.path.dirname(save):
                save = os.path.join("figures", save)
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()

    else:
        if library_key is not None:
            if library_id is None:
                raise ValueError("If library_key is provided, library_id must also be specified.")
            adata = adata[adata.obs[library_key] == library_id].copy()

        if "spatial_neighbors" not in adata.uns:
            sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="generic", key_added="spatial")

        score = adata.obs[score_key].values
        spatial = adata.obsm["spatial"]
        mask = ~np.isnan(score)
        score = score[mask]
        spatial = spatial[mask]

        W = KNN.from_array(spatial)
        W.transform = "r"
        moran = Moran(score, W)

        fig, ax = moran_scatterplot(moran, p=0.05, aspect_equal=True,
                                    scatter_kwds=scatter_kwds, fitline_kwds=fitline_kwds)
        fig.set_size_inches(3, 3)
        ax.set_title(f"Moran’s I = {moran.I:.2f}, p = {moran.p_sim:.3f}", fontsize=10)
        ax.set_ylabel("Spatial lag")
        ax.set_xlabel("Score")
        ax.grid(False)

        if save:
            os.makedirs("figures", exist_ok=True)
            if not os.path.dirname(save):
                save = os.path.join("figures", save)
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
