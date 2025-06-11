from __future__ import annotations

from anndata import AnnData
import squidpy as sq
import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = "truetype"


def spatial_enrichmap(
    adata: AnnData,
    score_key: str | list | None = None,
    cmap: str = "seismic",
    library_key: str = "library_id",
    library_id: str | list | None = None,
    img_alpha: float = 0.5,
    ncols: int | None = None,
    size: int = 2,
    vcenter: int = 0,
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
        Colormap for visualisation. Defaults to "seismic".

    library_key : str, optional
        Key in `adata.obs` for library identifiers. Defaults to "library_id".

    library_id : str or list of str, optional
        Specific library ID(s) to visualise. Defaults to all libraries.

    img_alpha : float, optional
        Alpha value for image overlay. Defaults to 0.5.

    ncols : int, optional
        Number of columns in the grid layout. Defaults to length of score_key.

    size : int, optional
        Scatter plot marker size. Defaults to 2.

    vcenter : int, optional
        Central value for colour scaling. Defaults to 0.

    save : str or None, optional
        Filename to save the plot. If None, plot is not saved.

    **kwargs : dict, optional
        Additional arguments passed to `sq.pl.spatial_scatter`.
    """
    if score_key is None:
        score_key = ["enrichmap_score"]
    elif isinstance(score_key, str):
        score_key = [score_key]

    if ncols is None:
        ncols = len(score_key)

    # Determine libraries to plot
    all_libs = adata.obs[library_key].unique()
    if library_id is None:
        libs_to_plot = all_libs
    elif isinstance(library_id, str):
        libs_to_plot = [library_id]
    else:
        libs_to_plot = library_id

    # Construct titles: one per combination of library and score key
    titles = []
    for lib in libs_to_plot:
        for score in score_key:
            clean_score = score.replace("_score", "")
            titles.append(f"{lib}: {clean_score}")

    sq.pl.spatial_scatter(
        adata,
        color=score_key,
        size=size,
        cmap=cmap,
        ncols=ncols,
        vcenter=vcenter,
        img_alpha=img_alpha,
        library_id=libs_to_plot,
        library_key=library_key,
        frameon=False,
        save=save,
        title=titles,
        **kwargs,
    )
