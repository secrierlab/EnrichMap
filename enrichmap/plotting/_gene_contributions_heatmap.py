from __future__ import annotations

import os
from typing import Literal

from anndata import AnnData
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
import squidpy as sq

plt.rcParams["pdf.fonttype"] = "truetype"


def gene_contributions_heatmap(
    adata: AnnData,
    score_key: str,
    top_n_genes: int = 5,
    bottom_n_genes: int = 5,
    cluster_genes: bool = True,
    order_spots: Literal["spatial", "cluster", "coherence"] | None = "coherence",
    groupby: str | None = None,
    spatial_key: str = "spatial",
    cmap: str = "RdBu_r",
    fontsize: int = 8,
    center: float = 0.0,
    batch_key: str | None = None,
    library_id: str | list | None = None,
    ncols: int = 2,
    figsize: tuple = (12, 6),
    save: str | None = None,
    n_neighbors: int = 6,
    group_palette: dict | None = None,
) -> None:
    """
    Spatial gene contribution heatmap with expression scaling and spatial ordering.
    """

    if "gene_contributions" not in adata.uns:
        raise ValueError(
            "Gene contributions not found in adata.uns. "
            "Run enrichmap.tools.score first."
        )

    contribution_matrix = adata.uns["gene_contributions"][score_key]

    def _select_extreme_genes(contributions, n_top, n_bottom):

        mean_abs = {
            gene: np.mean(np.abs(scores)) for gene, scores in contributions.items()
        }

        sorted_genes = sorted(
            mean_abs,
            key=mean_abs.get,
            reverse=True,
        )

        top = sorted_genes[:n_top]

        bottom = sorted_genes[-n_bottom:] if n_bottom > 0 else []

        bottom = [g for g in bottom if g not in top]

        return top + bottom, mean_abs

    def _cluster_order(matrix):

        if matrix.shape[0] < 2:
            return np.arange(matrix.shape[0])

        dist = pdist(
            matrix,
            metric="correlation",
        )

        dist = np.nan_to_num(
            dist,
            nan=0.0,
        )

        Z = linkage(
            dist,
            method="average",
        )

        return leaves_list(Z)

    def _spatial_order(coords):

        if coords is None or coords.shape[0] < 2:
            return np.arange(coords.shape[0])

        pc1 = PCA(n_components=1).fit_transform(coords).ravel()

        return np.argsort(pc1)

    def _scale_rows(X):

        mean = X.mean(
            axis=1,
            keepdims=True,
        )

        std = X.std(
            axis=1,
            keepdims=True,
        )

        std[std == 0] = 1.0

        return (X - mean) / std

    def _compute_spot_coherence(data, coords):

        from scipy.sparse import csr_matrix
        from anndata import AnnData as _AnnData

        signal = np.mean(
            np.abs(data),
            axis=0,
        )

        if coords is None or len(signal) < 3:
            return np.zeros_like(signal)

        ad = _AnnData(X=np.zeros((coords.shape[0], 1)))

        ad.obsm[spatial_key] = coords

        sq.gr.spatial_neighbors(
            ad,
            coord_type="generic",
            n_neighs=n_neighbors,
        )

        W = ad.obsp["spatial_connectivities"]

        if not isinstance(W, csr_matrix):
            W = csr_matrix(W)

        row_sum = np.array(W.sum(axis=1)).flatten()

        row_sum[row_sum == 0] = 1.0

        W = W.multiply(1 / row_sum[:, None])

        smoothed = W.dot(signal)

        return smoothed

    def _compute_gene_morans(
        contributions,
        genes,
        coords,
    ):

        from anndata import AnnData as _AnnData

        X = np.array([contributions[g] for g in genes]).T

        ad = _AnnData(X=X)

        ad.obsm[spatial_key] = coords
        ad.var_names = genes

        sq.gr.spatial_neighbors(
            ad,
            coord_type="generic",
            n_neighs=n_neighbors,
        )

        sq.gr.spatial_autocorr(
            ad,
            mode="moran",
        )

        return ad.uns["moranI"]["I"].to_dict()

    def _order_genes(
        contributions,
        genes,
        mean_abs,
        coords,
    ):

        if not cluster_genes:
            return genes

        moran = _compute_gene_morans(
            contributions,
            genes,
            coords,
        )

        ordered = sorted(
            genes,
            key=lambda g: (
                mean_abs[g],
                moran.get(g, 0.0),
            ),
            reverse=True,
        )

        X = np.array([contributions[g] for g in ordered])

        gene_order = _cluster_order(X)

        return [ordered[i] for i in gene_order]

    def _order_spots_within(data, coords):

        if order_spots == "spatial":
            return _spatial_order(coords)

        elif order_spots == "cluster":
            return _cluster_order(data.T)

        elif order_spots == "coherence":
            score = _compute_spot_coherence(
                data,
                coords,
            )

            return np.argsort(score)[::-1]

        return np.arange(data.shape[1])

    def _order_spots_grouped(
        data,
        coords,
        group_labels,
    ):

        unique_groups = list(dict.fromkeys(group_labels))

        final_order = []
        group_boundaries = []

        running = 0

        for grp in unique_groups:
            grp_mask = np.array(group_labels) == grp

            grp_indices = np.where(grp_mask)[0]

            if len(grp_indices) == 0:
                continue

            grp_data = data[:, grp_indices]

            grp_coords = coords[grp_indices] if coords is not None else None

            within_order = _order_spots_within(
                grp_data,
                grp_coords,
            )

            final_order.extend(grp_indices[within_order].tolist())

            running += len(grp_indices)

            group_boundaries.append((grp, running))

        return (
            np.array(final_order),
            unique_groups,
            group_boundaries,
        )

    def _build_heatmap(
        contributions,
        genes,
        mean_abs,
        coords,
        group_labels=None,
    ):

        ordered_genes = _order_genes(
            contributions,
            genes,
            mean_abs,
            coords,
        )

        X = np.array([contributions[g] for g in ordered_genes])

        X = _scale_rows(X)

        if group_labels is not None:
            (
                spot_order,
                groups,
                boundaries,
            ) = _order_spots_grouped(
                X,
                coords,
                group_labels,
            )

            X = X[:, spot_order]

            return (
                X,
                ordered_genes,
                groups,
                boundaries,
            )

        else:
            spot_order = _order_spots_within(
                X,
                coords,
            )

            X = X[:, spot_order]

            return (
                X,
                ordered_genes,
                None,
                None,
            )

    def _get_group_colours(groups):

        if group_palette is not None:
            return [group_palette.get(g, "#cccccc") for g in groups]

        palette = sns.color_palette(
            "tab20",
            n_colors=len(groups),
        )

        return [palette[i] for i in range(len(groups))]

    def _draw_heatmap(
        data,
        genes,
        ax,
        title,
        groups=None,
        boundaries=None,
    ):

        fig = ax.figure

        parent_spec = ax.get_subplotspec()

        ax.remove()

        inner = gridspec.GridSpecFromSubplotSpec(
            3,
            1,
            subplot_spec=parent_spec,
            height_ratios=[0.08, 0.06, 1],
            hspace=0.02,
        )

        ax_title = fig.add_subplot(inner[0])
        ax_group = fig.add_subplot(inner[1])
        ax_heat = fig.add_subplot(inner[2])

        # -------------------------------------------------
        # Title
        # -------------------------------------------------

        ax_title.text(
            0.5,
            0.5,
            title,
            ha="center",
            va="center",
            fontsize=fontsize + 2,
            fontweight="bold",
            transform=ax_title.transAxes,
        )

        ax_title.axis("off")

        # -------------------------------------------------
        # Group annotation bar
        # -------------------------------------------------

        if groups is not None and boundaries is not None:
            colours = _get_group_colours(groups)

            group_vector = np.zeros(
                data.shape[1],
                dtype=int,
            )

            start = 0

            for idx, (_, end) in enumerate(boundaries):
                group_vector[start:end] = idx

                start = end

            cmap_groups = ListedColormap(colours)

            ax_group.imshow(
                group_vector[np.newaxis, :],
                aspect="auto",
                cmap=cmap_groups,
                interpolation="nearest",
            )

            ax_group.set_xlim(0, data.shape[1])
            ax_heat.set_xlim(0, data.shape[1])
            ax_group.set_yticks([])
            ax_group.set_xticks([])

            ax_group.grid(False)

            for spine in ax_group.spines.values():
                spine.set_visible(False)

        else:
            ax_group.axis("off")

        # -------------------------------------------------
        # Heatmap
        # -------------------------------------------------

        vmax = np.nanpercentile(
            np.abs(data),
            99,
        )

        vmin = -vmax

        sns.heatmap(
            data,
            yticklabels=genes,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            xticklabels=False,
            ax=ax_heat,
        )

        xlabel = f"Spots ({order_spots or 'original order'})"

        if groups is not None:
            xlabel = (
                f"Spots (grouped by {groupby}, ordered by {order_spots or 'original'})"
            )

        ax_heat.set_xlabel(
            xlabel,
            fontsize=fontsize,
        )

        ax_heat.set_ylabel(
            "Genes (z-scored)",
            fontsize=fontsize,
        )

        ax_heat.tick_params(labelsize=fontsize)

        ax_heat.grid(False)

        # -------------------------------------------------
        # Group labels on bottom axis only
        # -------------------------------------------------

        if groups is not None and boundaries is not None:
            ticks = []
            labels = []

            start = 0

            for grp, end in boundaries:
                midpoint = (start + end) / 2

                ticks.append(midpoint)
                labels.append(str(grp))

                start = end

            ax_heat.set_xticks(ticks)

            ax_heat.set_xticklabels(
                labels,
                rotation=45,
                ha="right",
                fontsize=max(fontsize - 1, 5),
            )

    def _get_group_labels(mask=None):

        if groupby is None:
            return None

        if groupby not in adata.obs.columns:
            raise ValueError(f"'{groupby}' not found in adata.obs")

        labels = adata.obs[groupby].values

        if mask is not None:
            labels = labels[mask]

        return list(labels)

    # -------------------------------------------------
    # Single panel
    # -------------------------------------------------

    if batch_key is None:
        coords = adata.obsm.get(spatial_key)

        group_labels = _get_group_labels()

        genes, mean_abs = _select_extreme_genes(
            contribution_matrix,
            top_n_genes,
            bottom_n_genes,
        )

        (
            data,
            ordered_genes,
            groups,
            boundaries,
        ) = _build_heatmap(
            contribution_matrix,
            genes,
            mean_abs,
            coords,
            group_labels,
        )

        fig, ax = plt.subplots(
            figsize=figsize,
            constrained_layout=True,
        )

        _draw_heatmap(
            data,
            ordered_genes,
            ax,
            score_key,
            groups,
            boundaries,
        )

    # -------------------------------------------------
    # Multi-panel
    # -------------------------------------------------

    else:
        if batch_key not in adata.obs.columns:
            raise ValueError(f"{batch_key} not found in adata.obs")

        if library_id is not None:
            if isinstance(library_id, str):
                library_id = [library_id]

            unique_batches = library_id

        else:
            unique_batches = list(adata.obs[batch_key].unique())

        n_batches = len(unique_batches)

        n_rows = (n_batches + ncols - 1) // ncols

        fig, axes = plt.subplots(
            n_rows,
            ncols,
            figsize=(
                ncols * figsize[0] / 2,
                n_rows * figsize[1] / 2,
            ),
            constrained_layout=True,
        )

        axes = np.atleast_2d(axes).flatten()

        for i, batch in enumerate(unique_batches):
            mask = (adata.obs[batch_key] == batch).values

            coords = adata.obsm[spatial_key][mask]

            group_labels = _get_group_labels(mask)

            batch_contributions = {g: v[mask] for g, v in contribution_matrix.items()}

            (
                genes,
                mean_abs,
            ) = _select_extreme_genes(
                batch_contributions,
                top_n_genes,
                bottom_n_genes,
            )

            (
                data,
                ordered_genes,
                groups,
                boundaries,
            ) = _build_heatmap(
                batch_contributions,
                genes,
                mean_abs,
                coords,
                group_labels,
            )

            _draw_heatmap(
                data,
                ordered_genes,
                axes[i],
                f"{score_key} ({batch})",
                groups,
                boundaries,
            )

        for j in range(
            n_batches,
            len(axes),
        ):
            axes[j].axis("off")

    if save:
        os.makedirs(
            "figures",
            exist_ok=True,
        )

        if not os.path.dirname(save):
            save = os.path.join(
                "figures",
                save,
            )

        plt.savefig(
            save,
            dpi=300,
            bbox_inches="tight",
        )

    plt.show()
