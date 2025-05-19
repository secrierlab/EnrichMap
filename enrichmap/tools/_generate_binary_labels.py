from __future__ import annotations

import numpy as np

from anndata import AnnData
from typing import Literal

def generate_binary_labels(
    adata: AnnData,
    signature_name: str,
    threshold: Literal["median", "quantile"] | float = "median",
    quantile: float | None = None
) -> np.ndarray:
    """
    Generate binary labels for cells based on a gene signature score and a given threshold.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing signature scores in `adata.obs`.

    signature_name : str
        Name of the gene signature, corresponding to the key `{signature_name}_score` in `adata.obs`.

    threshold : {"median", "quantile"} or float, default "median"
        Method or value to determine the threshold for binary labelling.
        - "median": uses the median of the score distribution.
        - "quantile": uses a quantile specified by the `quantile` argument.
        - float: uses the given numeric threshold directly.

    quantile : float, optional
        Quantile to use when `threshold="quantile"`. Must be a float between 0 and 1.

    Returns
    -------
    np.ndarray
        A NumPy array of binary labels (1 for high score, 0 for low score), based on the specified threshold.
    """
    score = adata.obs[f"{signature_name}_score"]
    if threshold == "median":
        threshold_value = np.median(score)
    elif threshold == "quantile" and quantile is not None:
        threshold_value = np.quantile(score, quantile)
    else:
        threshold_value = threshold
    return (score > threshold_value).astype(int)
