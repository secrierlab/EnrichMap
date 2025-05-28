API Reference
=============

Main scoring function
---------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   enrichmap.tl.score

Helper functions
----------------

These functions are used by the ``score`` function.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   enrichmap.tl.generate_binary_labels
   enrichmap.tl.cluster_gene_correlation
   enrichmap.tl.infer_gene_weights


Main plotting function
----------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   enrichmap.pl.spatial_enrichmap

Geostatistical evaluations
--------------------------

Functionalities to assess the smoothness of scores or input gene signature.

.. autosummary::
   :nosignatures:
   :toctree: generated/

   enrichmap.pl.spatial_metrics
   enrichmap.pl.variogram
   enrichmap.pl.variogram_all
   enrichmap.pl.morans_correlogram
   enrichmap.pl.cross_moran_scatter

Signature associations
----------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   enrichmap.pl.signature_correlation_heatmap

Gene-level visualisations
-------------------------

.. autosummary::
   :nosignatures:
   :toctree: generated/

   enrichmap.pl.gene_contributions_heatmap
   enrichmap.pl.gene_contributions_pca
