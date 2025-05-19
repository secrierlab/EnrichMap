from __future__ import annotations

from anndata import AnnData
import squidpy as sq
import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = "truetype"

def spatial_enrichmap(
    adata: AnnData,
    score_key: str | list | None = None,
    cmap: str = "seismic",
    library_key: str = None,
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

    This function generates scatter plots for spatial enrichment scores of specified gene sets.
    It uses the `sq.pl.spatial_scatter` method to create the plots and allows for customisation of plot
    properties such as colour maps, marker size, and plot layout.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing spatial and gene expression information.

    score_key : str or list of str, optional
        A list of signature names (e.g. "enrichmap") for which enrichment scores are plotted. 
        If None, defaults to plotting the "enrichmap" signature. If a single string is provided, it will be 
        converted into a list with one element.

    cmap : str, optional
        The colormap for visualising the enrichment scores. Defaults to "seismic".

    library_key : str, optional
        The key in `adata.obs` that contains the library identifiers for different spatial libraries. 
        Defaults to "library_id" if not provided.

    library_id : str or list of str, optional
        A specific library ID or a list of library IDs to visualise. If None, visualises all libraries.

    img_alpha : float, optional
        The alpha blending value for image overlay. Defaults to 0.5.

    ncols : int, optional
        The number of columns for arranging the plots in a grid layout. Defaults to 2.

    size : int, optional
        The size of the scatter plot markers. Defaults to 2.

    vcenter : int, optional
        The central value for colour scaling. Defaults to 0.

    save : str or None, optional
        If specified, the plot will be saved to the file with the given filename (including extension). 
        If None, the plot will not be saved.

    **kwargs : dict, optional
        Additional keyword arguments passed to `sq.pl.spatial_scatter`.
        This can include parameters like `color_map`, `size`, etc.
        Refer to the `squidpy` documentation for more details on available parameters.
    """
    if score_key is None:
        score_key = ["enrichmap_score"]
    elif isinstance(score_key, str):
        score_key = [score_key]

    if ncols is None:
        ncols = len(score_key)

    sq.pl.spatial_scatter(
        adata,
        color=score_key,
        size=size,
        cmap=cmap,
        ncols=ncols,
        vcenter=vcenter,
        img_alpha=img_alpha,
        library_id=library_id,
        library_key=library_key,
        frameon=False,
        save=save,
        **kwargs
    )