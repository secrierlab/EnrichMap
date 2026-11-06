from __future__ import annotations

import logging
import warnings
from itertools import product

logging.getLogger("fontTools").setLevel(logging.WARNING)

from anndata import AnnData
import numpy as np
import squidpy as sq
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, PathCollection
from matplotlib.colors import CenteredNorm

from squidpy.pl._spatial_utils import _SeqStr
from squidpy.pl._utils import save_fig

plt.rcParams["pdf.fonttype"] = "truetype"


def _fdr_strokes(coll, significant: np.ndarray) -> PathCollection:
    """Black outlines, without fill, of the significant spots of a score scatter (points or shapes)."""
    style = {
        "facecolors": "none",
        "edgecolors": "black",
        "linewidths": 0.6,
        "transform": coll.get_transform(),
        "zorder": coll.get_zorder() + 0.1,
    }
    if isinstance(coll, PatchCollection):
        return PathCollection([p for p, s in zip(coll.get_paths(), significant) if s], **style)
    # points: the stroke follows the visible outline of the point, which includes its edge only when that is drawn
    sizes = np.broadcast_to(coll.get_sizes(), significant.shape)[significant]
    edges = coll.get_edgecolors()
    edge = np.max(coll.get_linewidths()) if len(edges) and np.max(edges[:, 3]) > 0 else 0.0
    grow = edge / coll.get_paths()[0].get_extents().width
    return PathCollection(
        coll.get_paths(),
        sizes=(np.sqrt(sizes) + grow) ** 2,
        offsets=coll.get_offsets()[significant],
        offset_transform=coll.get_offset_transform(),
        **style,
    )


def spatial_enrichmap(
    adata: AnnData,
    score_key: str | list | None = None,
    cmap: str = "RdBu_r",
    library_key: str | None = "library_id",
    library_id: _SeqStr | None = None,
    img_alpha: float = 0.5,
    ncols: int | None = None,
    size: int = 1.5,
    vcenter: float | None = 0,
    plot_fdr: bool = False,
    fdr_threshold: float = 0.05,
    save: str | None = None,
    **kwargs: dict,
) -> None:
    """
    Visualise spatial enrichment maps for given signatures using spatial scatter plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing spatial and gene expression information.
    score_key : str or list of str, optional
        A list of signature names for which enrichment scores are plotted.
        Defaults to ["enrichmap_score"].
    cmap : str, optional
        Colormap for visualisation. Defaults to "RdBu_r".
    library_key : str or None, optional
        Key in `adata.obs` for library identifiers. Defaults to "library_id".
    library_id : str or list of str, optional
        Specific library ID(s) to visualise. Defaults to all libraries.
    img_alpha : float, optional
        Alpha value for image overlay. Defaults to 0.5.
    ncols : int, optional
        Number of columns in the grid layout. Defaults to length of score_key.
    size : int, optional
        Scatter plot marker size. Defaults to 1.5.
    vcenter : float or None, optional
        Central value for colour scaling. Defaults to 0. Each panel is scaled symmetrically around it
        (vcenter ± the largest distance from vcenter), so equal shades mean equal distances on both sides.
        If `vmin` or `vmax` is given, the range is mirrored from the one further from vcenter
        (`vmin=-4, vmax=10` draws -10 to 10). None scales each panel from its minimum to its maximum,
        or uses `vmin` and `vmax` exactly as given.
    plot_fdr : bool, optional
        Whether to draw a black stroke around spots called significant, i.e. with q-value below
        `fdr_threshold` in `adata.obs[f"{signature}_q"]` (written by
        ``em.tl.score(..., compute_q_vals=True)``). Defaults to False.
    fdr_threshold : float, optional
        FDR level for `plot_fdr`. Defaults to 0.05.
    save : str or None, optional
        Filename to save the plot. If None, plot is not saved.
    **kwargs : dict, optional
        Additional arguments passed to `sq.pl.spatial_scatter`.
    """

    # --- Validate and normalise score keys ---
    if score_key is None:
        score_key = ["enrichmap_score"]
    elif isinstance(score_key, str):
        score_key = [score_key]

    if ncols is None:
        ncols = len(score_key)

    # --- Determine libraries to plot ---
    if library_key is None or library_key not in adata.obs.columns:
        # Single-library mode
        libs_to_plot = [None]
        library_key = None
    else:
        all_libs = adata.obs[library_key].unique().tolist()
        if library_id is None:
            libs_to_plot = all_libs
        elif isinstance(library_id, str):
            libs_to_plot = [library_id]
        else:
            libs_to_plot = library_id

    # --- Construct titles for plots ---
    titles = []
    for lib in libs_to_plot:
        for score in score_key:
            clean_score = score.replace("_score", "")
            if lib is not None:
                titles.append(f"{lib}: {clean_score}")
            else:
                titles.append(clean_score)

    # --- Symmetric colour scale around vcenter ---
    # squidpy turns vcenter into a TwoSlopeNorm with vmin and vmax or the data min and max, which stretches
    # the smaller side over half of the colour map (e.g. background z of -2 drawn as dark as a peak of +30)
    if vcenter is not None and "norm" not in kwargs:
        limits = [v for v in (kwargs.pop("vmin", None), kwargs.pop("vmax", None)) if v is not None]
        halfrange = max(abs(v - vcenter) for v in limits) if limits else None
        kwargs["norm"] = CenteredNorm(vcenter=vcenter, halfrange=halfrange)
        vcenter = None

    # --- q-value columns for plot_fdr ---
    q_keys = {}
    if plot_fdr:
        q_keys = {s: (s[: -len("_score")] if s.endswith("_score") else s) + "_q" for s in score_key}
        missing = [s for s, q in q_keys.items() if q not in adata.obs.columns]
        if len(missing) == len(score_key):
            raise ValueError(
                "plot_fdr=True needs q-values in adata.obs (e.g. 'enrichmap_q'); "
                "score with em.tl.score(..., compute_q_vals=True)."
            )
        if missing:
            warnings.warn(f"No q-values for {missing}; these panels are drawn without FDR strokes.", stacklevel=2)

    # --- Generate spatial plots ---
    axes = sq.pl.spatial_scatter(
        adata,
        color=score_key,
        size=size,
        cmap=cmap,
        ncols=ncols,
        vcenter=vcenter,
        img_alpha=img_alpha,
        library_id=None if libs_to_plot == [None] else libs_to_plot,
        library_key=library_key,
        frameon=False,
        save=None if plot_fdr else save,
        title=titles,
        return_ax=plot_fdr,
        **kwargs,
    )

    # --- Black stroke around spots with q < fdr_threshold ---
    if plot_fdr:
        axes = list(axes) if isinstance(axes, (list, tuple, np.ndarray)) else [axes]
        # panels are drawn library by library, then signature by signature, with spots in adata order
        for ax, (lib, score) in zip(axes, product(libs_to_plot, score_key)):
            if q_keys[score] not in adata.obs.columns:
                continue
            in_lib = np.ones(adata.n_obs, bool) if lib is None else (adata.obs[library_key] == lib).values
            significant = adata.obs[q_keys[score]].values[in_lib] < fdr_threshold
            # the score scatter is the last collection with one element per spot; the strokes go in their own
            # collection, as in a PDF even a zero-width edge on the other spots would be drawn (as a hairline)
            for coll in reversed(ax.collections):
                n = len(coll.get_paths()) if isinstance(coll, PatchCollection) else len(coll.get_offsets())
                if n == in_lib.sum():
                    ax.add_collection(_fdr_strokes(coll, significant), autolim=False)
                    break
        if save:
            save_fig(axes[0].figure, path=save)
