Welcome to EnrichMap's documentation!
=====================================

**EnrichMap** is a lightweight tool designed to compute and visualise enrichment scores of a given gene set or signature in spatial transcriptomics datasets across different platforms. It offers flexible scoring, batch correction, spatial smoothing and visual outputs for intuitive exploration of biological signatures.

.. image:: https://github.com/secrierlab/enrichmap/raw/main/img/enrichmap_workflow.jpg
   :alt: EnrichMap workflow
   :align: center

Features
--------

- Fast computation of enrichment scores
- Support for batch correction and spatial covariates
- Built-in spatial smoothing
- Visualisation tools for intuitive mapping
- Easy integration with AnnData (.h5ad) objects

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

This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Contents
   
   install
   usage
   tutorials/index
   api