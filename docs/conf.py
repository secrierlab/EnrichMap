# Configuration file for the Sphinx documentation builder.

import os
import sys

# -- Path setup --------------------------------------------------------------

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------

project = "EnrichMap"
author = "Cenk Celik"
release = "0.1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "myst_nb",
]

autosummary_generate = True

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

# Prevent Sphinx from importing heavy scientific dependencies
autodoc_mock_imports = [
    "scanpy",
    "anndata",
    "numpy",
    "scipy",
    "sklearn",
    "matplotlib",
    "pandas",
    "seaborn",
    "pygam",
    "scikit_gstat",
    "squidpy",
    "spatialdata",
    "xarray",
    "xarray_schema",
    "dask",
]

# -- MyST / Notebook configuration ------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
]

# Disable execution of notebooks during build
myst_nb_execute = "off"

# -- Autodoc configuration ---------------------------------------------------

autodoc_member_order = "bysource"
autodoc_typehints = "description"

# -- HTML output -------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
