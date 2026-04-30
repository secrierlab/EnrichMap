"""Configuration for EnrichMap's Sphinx documentation."""

from __future__ import annotations
from pathlib import Path
from functools import partial
from docutils import nodes
from typing import TYPE_CHECKING

HERE = Path(__file__).parent

import os
import sys
from unittest.mock import MagicMock

sys.path.insert(0, os.path.abspath(".."))

# Mock heavy/broken dependencies before importing enrichmap so that:
# 1. enrichmap can be imported successfully (registering enrichmap.pl / enrichmap.tl)
# 2. xarray_schema's broken `from pkg_resources import ...` never executes
_MOCK_MODULES = [
    "squidpy",
    "squidpy.pl",
    "squidpy.pl._spatial_utils",
    "spatialdata",
    "spatialdata._logging",
    "xarray_schema",
    "anndata",
    "scanpy",
    "numpy",
    "scipy",
    "scipy.sparse",
    "scipy.spatial",
    "scipy.spatial.distance",
    "scipy.cluster",
    "scipy.cluster.hierarchy",
    "scipy.stats",
    "pandas",
    "matplotlib",
    "matplotlib.pyplot",
    "matplotlib.gridspec",
    "matplotlib.colors",
    "matplotlib.patches",
    "matplotlib.cm",
    "mpl_toolkits",
    "mpl_toolkits.axes_grid1",
    "seaborn",
    "sklearn",
    "sklearn.decomposition",
    "sklearn.metrics",
    "tqdm",
    "pygam",
    "libpysal",
    "libpysal.weights",
    "libpysal.weights.spatial_lag",
    "esda",
    "esda.moran",
    "esda.geary",
    "esda.getisord",
    "splot",
    "splot.esda",
    "skgstat",
    "adjustText",
]
for _mod in _MOCK_MODULES:
    sys.modules[_mod] = MagicMock()

import enrichmap  # noqa: E402  — must come after mocks

if TYPE_CHECKING:
    from sphinx.application import Sphinx

# -- Project information
project = "EnrichMap"
copyright = "2026, Cenk Celik"
author = "Cenk Celik"

release = enrichmap.__version__
version = ".".join(release.split(".")[:2])

master_doc = "index"

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "myst_nb",
]

nb_execution_mode = "off"

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "deflist",
    "html_admonition",
    "html_image",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

# -- Options for EPUB output

epub_show_urls = "footnote"

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_mock_imports = [
    "squidpy",
    "spatialdata",
    "xarray_schema",
    "anndata",
    "scanpy",
    "numpy",
    "scipy",
    "pandas",
    "matplotlib",
    "seaborn",
    "sklearn",
    "tqdm",
    "pygam",
    "libpysal",
    "esda",
    "splot",
    "skgstat",
    "adjustText",
]

# -- Options for HTML output ----------------------------------------------
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_logo = "_static/logo-light.svg"

# Use light and dark logos for sphinx-book-theme
html_theme_options = {
    "repository_url": "https://github.com/secrierlab/enrichmap",
    "use_repository_button": True,
    "use_issues_button": True,
    "logo": {"image_dark": "_static/logo-dark.svg"},
}

html_show_sphinx = False
html_title = "EnrichMap"
html_favicon = "_static/favicon.png"


def setup(app: Sphinx):
    """App setup hook."""
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))


# -- Options for other output formats ------------------------------------------

htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        "One line description of project.",
        "Miscellaneous",
    )
]
