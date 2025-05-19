from __future__ import annotations

import os
from anndata import AnnData
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = "truetype"

def gene_contributions_heatmap(
    adata: AnnData,
    score_key: str,
    top_n_genes: int = 10,
    cmap: str = "seismic",
    fontsize: int = 8,
    center: int = 0,
    batch_key: str | None = None,
    library_id: str | list | None = None,
    ncols: int = 2,
    figsize: tuple = (12, 6),
    save: str | None = None
) -> None:
    """
    Visualise gene contributions from spatial gene set scoring as heatmaps.

    This function displays a heatmap of the top contributing genes for a given gene set signature.
    If `batch_key` is provided, one subplot is generated per batch (or specified `library_id`), enabling
    batch-wise comparison. Otherwise, a single heatmap is shown with mean contributions.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Gene contributions must be stored in `adata.uns["gene_contributions"]`.

    score_key : str
        Name of the gene signature as used in `enrichmap.tools.score`.

    top_n_genes : int, default=10
        Number of top contributing genes to display in the heatmap.

    cmap : str, default="seismic"
        Colour map for the heatmap (passed to seaborn.heatmap).

    fontsize : int, default=8
        Font size used for labels and titles.

    center : int, default=0
        Value at which to centre the colormap. Typically 0 for contributions.

    batch_key : str, optional
        Column in `adata.obs` identifying batches or libraries. Required to produce batch-wise subplots.

    library_id : str or list of str, optional
        Specific library IDs (within `batch_key`) to include in the visualisation. If not specified, all batches are included.

    ncols : int, default=2
        Number of columns when displaying multiple heatmaps (used when `batch_key` is set).

    figsize : tuple, default=(12, 6)
        Size of the overall figure for plotting.

    save : str, optional
        If provided, path to save the resulting figure (e.g. "figures/heatmap.pdf").

    Returns
    -------
    None
        Displays heatmaps and optionally saves the plot.

    Notes
    -----
    - The function assumes `spatial_gene_set_scoring` has been run and results stored in `adata.uns["gene_contributions"]`.
    - Mean contributions are used to identify top genes per batch or globally.
    - The resulting heatmaps show spatial spot-level contributions per gene, not aggregated statistics.
    """
    if "gene_contributions" not in adata.uns:
        raise ValueError("Gene contributions not found in adata.uns. Run spatial_gene_set_scoring first.")
    
    contribution_matrix = adata.uns["gene_contributions"][score_key]
    
    if batch_key is None:
        batch_contributions = {gene: scores for gene, scores in contribution_matrix.items()}
    
        mean_contributions = {gene: np.mean(scores) for gene, scores in batch_contributions.items()}
        sorted_genes = sorted(mean_contributions, key=mean_contributions.get, reverse=True)[:top_n_genes]
        
        heatmap_data = np.array([batch_contributions[gene] for gene in sorted_genes])
        
        plt.figure(figsize=figsize)
        sns.heatmap(heatmap_data, yticklabels=sorted_genes, cmap=cmap, center=center, annot=False, xticklabels=False)
        plt.xlabel("Spots", fontsize=fontsize)
        plt.ylabel("Top contributing genes", fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.title(f"{score_key}", fontsize=fontsize)
        plt.grid(False)
        if save:
        # Ensure 'figures/' directory exists
            os.makedirs("figures", exist_ok=True)

            # If 'save' has no directory path, prepend 'figures/'
            if not os.path.dirname(save):
                save = os.path.join("figures", save)
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
        return
    
    if batch_key not in adata.obs.columns:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs.")
    
    # Handle the case where library_id is a string or a list of strings
    if library_id is not None:
        if isinstance(library_id, str):
            library_id = [library_id]  # Convert to list if it"s a single string
        # Check that library_id values exist
        if not all(item in adata.obs[batch_key].values for item in library_id):
            raise ValueError(f"Some library ids in '{library_id}' not found in adata.obs[{batch_key}].")
        
        # Filter based on values in library_id
        unique_batches = [batch for batch in library_id if batch in adata.obs[batch_key].values]
    else:
        unique_batches = adata.obs[batch_key].unique()
    
    # Determine the number of rows for subplots
    n_batches = len(unique_batches)
    n_rows = (n_batches + ncols - 1) // ncols  # Calculate number of rows to fit all subplots
    
    fig, axes = plt.subplots(n_rows, ncols, figsize=(ncols * 3, n_rows * 3), constrained_layout=True)
    axes = axes.flatten()  # Flatten to easily iterate through the axes
    
    for i, batch in enumerate(unique_batches):
        batch_mask = adata.obs[batch_key] == batch
        
        batch_contributions = {gene: scores[batch_mask] for gene, scores in contribution_matrix.items()}
        
        mean_contributions = {gene: np.mean(scores) for gene, scores in batch_contributions.items()}
        sorted_genes = sorted(mean_contributions, key=mean_contributions.get, reverse=True)[:top_n_genes]
        
        heatmap_data = np.array([batch_contributions[gene] for gene in sorted_genes])
        
        sns.heatmap(heatmap_data, yticklabels=sorted_genes, cmap=cmap, center=center, annot=False, xticklabels=False, ax=axes[i])
        axes[i].set_xlabel("Spots", fontsize=fontsize)
        axes[i].set_ylabel("Top contributing genes", fontsize=fontsize)
        axes[i].tick_params(axis="both", labelsize=fontsize)
        axes[i].set_title(f"Top {top_n_genes} contributing genes for {score_key} ({batch})", fontsize=fontsize)
    
    # Hide any unused subplots (if there are fewer batches than subplots)
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    plt.grid(False)
    if save:
        # Ensure 'figures/' directory exists
        os.makedirs("figures", exist_ok=True)

        # If 'save' has no directory path, prepend 'figures/'
        if not os.path.dirname(save):
            save = os.path.join("figures", save)
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()