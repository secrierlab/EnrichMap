Basic tutorial
==============================================

.. code-block:: python

   import scanpy as sc
   import enrichmap as em

   # Load your AnnData object
   adata = sc.read_h5ad("PATH/TO/YOUR/DATA.h5ad")

   # Define a gene set
   gene_set = ["CD3D", "CD3E", "CD8A"]

   # Run scoring
   em.tl.score(
       adata=adata,
       gene_list=gene_set,
       score_key="T_cell_signature",
       smoothing=True,  # by default
       correct_spatial_covariates=True,  # by default
       batch_key=None  # Set batch_key if working with multiple slides
   )

   # Visualise
   em.pl.spatial_enrichmap(
       adata=adata,
       score_key="T_cell_signature"
   )

.. note::

   EnrichMap currently does not support reading in ``SpatialData`` format. However, users can simply convert ``SpatialData`` to legacy ``AnnData`` to use EnrichMap.

.. code-block:: python

   import spatialdata_io as sd
   # Read in SpatialData
   sdata = sd.visium_hd("PATH_TO_DATA_FOLDER/")
   # Convert to AnnData
   adata = to_legacy_anndata(
       sdata,
       include_images=True,
       table_name="square_008um",
       coordinate_system="downscaled_hires"
   )