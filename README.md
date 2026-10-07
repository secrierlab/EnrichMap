# EnrichMap: Spatially-aware gene set enrichment

<p align="center">
  <img src="https://github.com/secrierlab/enrichmap/raw/main/img/enrichmap_logo.jpg" alt="EnrichMap" width="400" />
</p>

`EnrichMap` is a lightweight tool designed to compute and visualise enrichment scores of a given gene set or signature in spatial transcriptomics datasets across different platforms. It offers flexible scoring, batch correction, spatial smoothing and visual outputs for intuitive exploration of biological signatures.

<img src="https://github.com/secrierlab/enrichmap/raw/main/img/enrichmap_workflow.jpg" alt="EnrichMap workflow" style="width: auto; height: auto;">

## Features

- Fast computation of enrichment scores
- Support for batch correction and spatial covariates
- Built-in spatial smoothing
- Visualisation tools for intuitive mapping
- Easy integration with AnnData (.h5ad) objects

## System requirements

- **Operating system:** developed and tested on macOS 26 (Apple M2 Pro, 16 GB RAM) and Ubuntu 22.04 Server. `EnrichMap` is a pure-Python package with pip-installable dependencies, so it is expected to run on Linux and Windows as well, but these have not been tested.
- **Python:** ≥ 3.10. Tested with Python 3.11.9 and 3.12.4.
- **Hardware:** no non-standard hardware is required (no GPU). The demo below peaks at about 1.4 GB of RAM; memory use grows with the size of the dataset.
- **Dependencies:** installed automatically by `pip`. `EnrichMap` 0.2.4 was tested with the following resolved versions (dependencies are not pinned; versions resolved by `pip` on 7 October 2026):

| Package | Python 3.11.9 | Python 3.12.4 |
|---|---|---|
| numpy | 2.4.6 | 2.4.6 |
| pandas | 2.3.3 | 3.0.6 |
| scipy | 1.16.3 | 1.16.3 |
| scikit-learn | 1.9.1 | 1.9.1 |
| scanpy | 1.11.5 | 1.12.4 |
| anndata | 0.12.19 | 0.13.4 |
| squidpy | 1.8.2 | 1.8.3 |
| matplotlib | 3.11.2 | 3.11.2 |
| seaborn | 0.13.2 | 0.13.2 |
| dask | 2026.1.1 | 2026.8.0 |
| pygam | 0.12.0 | 0.12.0 |
| POT | 0.9.7.post1 | 0.9.7.post1 |
| esda | 2.9.0 | 2.10.0 |
| libpysal | 4.14.1 | 4.15.0 |
| scikit-gstat | 1.0.24 | 1.0.24 |
| statannotations | 0.7.2 | 0.7.2 |
| adjustText | 1.4.0 | 1.4.0 |

## Installation

A `conda` environment is strongly recommended with `python` ≥ 3.10.

```bash
conda create -n enrichmap_env python=3.11
conda activate enrichmap_env
```

Then, install `enrichmap` via `pip`.

```bash
pip install enrichmap
```

or directly from GitHub:

```bash
pip install git+https://github.com/secrierlab/enrichmap.git
```

Typical install time on a normal desktop computer: about 5 minutes (measured 4.5–5 minutes into a fresh environment with an empty `pip` cache, mostly spent downloading the dependencies).

## Basic usage

```python
import scanpy as sc
import enrichmap as em

# Load your AnnData object
adata = sc.read_h5ad("PATH/TO/YOUR/DATA.h5ad")

# Define a gene set
gene_set = ["CD3D", "CD3E", "CD8A"]

# Run scoring
em.tl.score(
    adata=adata,
    gene_set=gene_set,
    score_key="T_cell_signature",
    smoothing=True,  # by default,
    correct_spatial_covariates=True,  # by default
    batch_key=None,  # Set batch_key if working with multiple slides
)

# Visualise
em.pl.spatial_enrichmap(adata=adata, score_key="T_cell_signature_score")
```

> Important note: EnrichMap currently does not support reading in `SpatialData` format. However, users can simply convert `SpatialData`  to legacy `AnnData` to use EnrichMap.
```python
import spatialdata_io as sd

# read in SpatialData
sdata = sd.visium_hd("PATH_TO_DATA_FOLDER/")
# convert to AnnData
adata = to_legacy_anndata(
    sdata,
    include_images=True,
    table_name="square_008um",
    coordinate_system="downscaled_hires",
)
```

## Demo

The demo scores the hybrid epithelial-to-mesenchymal transition (EMT) signature ([Malagoli Tagliazucchi et al., 2023](https://doi.org/10.1038/s41467-023-36439-7)) on a small spatial transcriptomics dataset of two breast cancer Visium slides (5,870 spots × 14,664 genes, raw counts, 72 MB). The dataset is included in this repository at `tests/dataset/adata_breast.h5ad` and is downloaded automatically by the code below. It corresponds to `adata_slides_0_3.h5ad` in the Zenodo record [10.5281/zenodo.15438169](https://doi.org/10.5281/zenodo.15438169) (CC BY 4.0), with the metadata reduced to keep the file small.

**Instructions.** Install `EnrichMap` (see above), then run the following code as a Python script or in a notebook:

```python
import scanpy as sc
import enrichmap as em

# Download the demo dataset (72 MB; cached in the working directory)
adata = sc.read(
    "adata_breast.h5ad",
    backup_url="https://github.com/secrierlab/EnrichMap/raw/main/tests/dataset/adata_breast.h5ad",
)
print(
    f"{adata.n_obs} spots x {adata.n_vars} genes; slides: {sorted(adata.obs['batch'].unique())}"
)

# The data are raw counts: normalise and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Hybrid EMT signature
hybrid = [
    "PDPN",
    "ITGA5",
    "ITGA6",
    "TGFBI",
    "LAMC2",
    "MMP10",
    "LAMA3",
    "CDH13",
    "SERPINE1",
    "P4HA2",
    "TNC",
    "MMP1",
]

# Score; batch_key makes smoothing and spatial correction run per slide
em.tl.score(adata, gene_set=hybrid, score_key="Hybrid", batch_key="batch")
print(adata.obs["Hybrid_score"].describe().round(3).to_string())

# Visualise both slides and save the figure to figures/demo_hybrid_score.png
em.pl.spatial_enrichmap(
    adata,
    score_key=["Hybrid_score"],
    size=2,
    library_key="batch",
    library_id=["0", "3"],
    cmap="RdBu_r",
    shape=None,
    save="demo_hybrid_score.png",
)
```

**Expected output.** The console prints:

```
5870 spots x 14664 genes; slides: ['0', '3']
Scoring Hybrid: 12/12 genes found: 100%|██████████| 1/1 [...]
count    5870.000
mean       -0.000
std         0.460
min        -1.498
25%        -0.295
50%        -0.049
75%         0.236
max         3.486
```

and a two-panel figure is saved (one spatial map per slide, spots coloured by the `Hybrid` enrichment score on a diverging red-blue scale, centred at zero). The scores are stored in `adata.obs["Hybrid_score"]`. Results were identical on Python 3.11.9 and 3.12.4, and with both the PyPI release (0.2.3) and the current GitHub version (0.2.4). Depending on your `squidpy` and `anndata` versions, harmless `FutureWarning` deprecation messages may also be printed.

**Expected run time.** About 10 seconds on a normal desktop computer once the dataset is cached (the scoring itself takes under 1 second). The very first run after installation takes longer (about 1 minute in our tests) because of one-off library start-up and compilation costs, plus the 72 MB dataset download.

## Reproducing the manuscript results

The notebooks used to generate the figures and supplementary figures are available in a separate repository: [secrierlab/enrichmap-reproducibility](https://github.com/secrierlab/enrichmap-reproducibility). Its README lists, for each figure, the notebook that produces it, covering the simulated benchmarks, Visium (mouse brain and human breast cancer), Visium HD, Xenium, MERFISH and imaging mass cytometry analyses. Figure numbering in that repository follows the bioRxiv preprint. Its README describes the environment (`requirements.txt`), where to obtain each dataset (a Zenodo record, built-in `squidpy` datasets and public downloads) and the order in which to run the notebooks.

## Documentation

Comprehensive documentation is available at:
https://enrichmap.readthedocs.io/en/latest

## Contributing

If you have ideas for new features or spot a bug, please open an issue or submit a pull request.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0-only); see the [LICENSE](LICENSE) file.

## Citation

Celik C & Secrier M (2025). EnrichMap: Spatially-informed enrichment analysis for functional interpretation of spatial transcriptomics. [biorxiv.com](https://www.biorxiv.org/content/10.1101/2025.05.30.656960v1)

### Copyright

This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.