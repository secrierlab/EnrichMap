Welcome to EnrichMap's documentation!
=====================================

**EnrichMap** performs gene set enrichment analysis on spatial transcriptomics data, with gene weight inference and geostatistical modelling of spatial structure. It enables pathway-level mapping and comparison of opposing programmes across platforms, with outputs designed for direct interpretation.

.. image:: https://github.com/secrierlab/enrichmap/raw/main/img/enrichmap_workflow.jpg
   :alt: EnrichMap workflow
   :align: center

Features
--------

- Efficient computation of gene set enrichment scores in spatial data
- Gene weight inference within pathways and across opposing programmes
- Spatial smoothing and correction for spatial covariates
- Optional batch correction (recommended mode) for multi-sample integration
- Geostatistical diagnostics to quantify spatial structure
- Compatibility with major platforms (e.g. Visium, Visium HD, Xenium, MERFISH, etc.)
- Native support for AnnData (.h5ad) objects
- Direct visualisation of spatial enrichment maps

Documentation
-------------

Comprehensive documentation is available at:
https://enrichmap.readthedocs.io/en/latest

Contributing
------------

If you have ideas for new features or spot a bug, please open an issue or submit a pull request.

License
-------

This project is licensed under the GNU GENERAL PUBLIC LICENSE.

Citation
--------

Celik C & Secrier M (2025). *EnrichMap: Spatially-informed enrichment analysis for functional interpretation of spatial transcriptomics*. `biorxiv.com <https://www.biorxiv.org/content/10.1101/2025.05.30.656960v1>`_

Copyright
---------

This program is free software: you can redistribute it and/or modify it under the terms of the `GNU General Public License <https://github.com/secrierlab/EnrichMap/blob/main/LICENSE>`_ as published by the Free Software Foundation, either version 3 of the License or any later version.
This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Contents
   
   install
   usage
   tutorials/index
   api