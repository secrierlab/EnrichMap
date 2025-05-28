API Reference
=============

.. currentmodule:: enrichmap

Tools: ``tl``
=============

.. module:: enrichmap.tl

Main scoring function
---------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.score

Helper functions
----------------

These functions are used by the ``score`` function.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   tl.generate_binary_labels
   tl.cluster_gene_correlation
   tl.infer_gene_weights


Plotting: ``pl``
================

.. module:: enrichmap.pl

Main plotting function
----------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.spatial_enrichmap

Geostatistical evaluations
--------------------------

Functionalities to assess the smoothness of scores or input gene signature.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.spatial_metrics
   pl.variogram
   pl.variogram_all
   pl.morans_correlogram
   pl.cross_moran_scatter

Signature associations
----------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.signature_correlation_heatmap

Gene-level visualisations
-------------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   pl.gene_contributions_heatmap
   pl.gene_contributions_pca
