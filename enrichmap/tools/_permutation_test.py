from __future__ import annotations

from anndata import AnnData
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["pdf.fonttype"] = "truetype"

def permutation_test(
    adata: AnnData,
    signature_name: str = "enrichmap",
    fig_size: tuple = (8, 6),
    n_permutations: int = 100
) -> AnnData | None:
    """
    Perform a permutation test to assess the significance of a gene signature score.

    This test compares observed signature scores with scores obtained from randomly selected gene sets,
    providing an empirical p-value for each cell.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing signature scores in `adata.obs`.

    signature_name : str, default "enrichmap"
        Name of the gene signature, corresponding to the key `{signature_name}_score` in `adata.obs`.

    fig_size : tuple, default (8, 6)
        Size of the histogram plot showing the distribution of permutation scores.

    n_permutations : int, default 100
        Number of random permutations to perform for generating the null distribution.

    Returns
    -------
    AnnData
        The original `AnnData` object with an additional column `{signature_name}_pval` in `adata.obs`,
        containing empirical p-values for each cell.
    """

    observed_scores = adata.obs[f"{signature_name}_score"].values
    permuted_scores = np.zeros((len(observed_scores), n_permutations))
    
    for i in range(n_permutations):
        shuffled_genes = np.random.permutation(adata.var_names)
        random_scores = adata[:, shuffled_genes[:len(observed_scores)]].X.mean(axis=1)
        permuted_scores[:, i] = random_scores.flatten()

    # Compute p-values
    p_values = np.mean(permuted_scores >= observed_scores[:, None], axis=1)

    # Ensure correct shape for assignment
    if len(p_values) != len(adata.obs):
        raise ValueError(f"Mismatch between computed p-values ({len(p_values)}) and number of cells ({len(adata.obs)}).")

    adata.obs[f"{signature_name}_pval"] = p_values

    # Visualisation
    plt.figure(figsize=fig_size)
    sns.histplot(permuted_scores.flatten(), bins=50, kde=True, color='grey', alpha=0.5)
    plt.axvline(observed_scores.mean(), color='red', linestyle='dashed', linewidth=2, label='Observed mean')
    plt.xlabel("Permutation scores")
    plt.ylabel("Frequency")
    plt.legend()
    plt.title(f"Permutation test for {signature_name}")
    plt.show()