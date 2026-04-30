"""Configuration for EnrichMap's Sphinx documentation."""

from __future__ import annotations

import os
import sys
from functools import partial
from typing import TYPE_CHECKING
from unittest.mock import MagicMock

from docutils import nodes

if TYPE_CHECKING:
    from sphinx.application import Sphinx

sys.path.insert(0, os.path.abspath(".."))


_MOCK_MODULES = [
    # scverse ecosystem
    "anndata",
    "scanpy",
    "squidpy",
    "squidpy.pl",
    "squidpy.pl._spatial_utils",
    "spatialdata",
    "spatialdata._logging",
    # numerics
    "numpy",
    "pandas",
    "scipy",
    "scipy.cluster",
    "scipy.cluster.hierarchy",
    "scipy.sparse",
    "scipy.spatial",
    "scipy.spatial.distance",
    "scipy.stats",
    # plotting
    "matplotlib",
    "matplotlib.cm",
    "matplotlib.colors",
    "matplotlib.gridspec",
    "matplotlib.patches",
    "matplotlib.pyplot",
    "mpl_toolkits",
    "mpl_toolkits.axes_grid1",
    "seaborn",
    "adjustText",
    # ML / stats
    "sklearn",
    "sklearn.decomposition",
    "sklearn.metrics",
    "pygam",
    # spatial statistics
    "libpysal",
    "libpysal.weights",
    "libpysal.weights.spatial_lag",
    "esda",
    "esda.geary",
    "esda.getisord",
    "esda.moran",
    "skgstat",
    # misc
    "tqdm",
]


class _SmartMock(MagicMock):
    """Mock that preserves class names so Sphinx renders clean type hints."""

    def __getattr__(self, name):
        parent = self._mock_name or ""
        child = _SmartMock(name=f"{parent}.{name}" if parent else name)
        child.__name__ = name
        child.__qualname__ = name
        child.__module__ = parent
        return child

    def __repr__(self):
        return getattr(self, "__name__", "Mock")


for _mod in _MOCK_MODULES:
    sys.modules[_mod] = _SmartMock(name=_mod)

import enrichmap

project = "EnrichMap"
copyright = "2026, Cenk Celik"  # noqa: A001
author = "Cenk Celik"
release = enrichmap.__version__
version = ".".join(release.split(".")[:2])
master_doc = "index"

# General configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.duration",
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

# Autodoc

autosummary_generate = True
autodoc_member_order = "bysource"

# Use the same list for autodoc_mock_imports (top-level packages only)
autodoc_mock_imports = sorted({mod.split(".")[0] for mod in _MOCK_MODULES})

# HTML output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_logo = "_static/logo-light.svg"
html_favicon = "_static/favicon.png"
html_title = "EnrichMap"
html_show_sphinx = False

html_theme_options = {
    "repository_url": "https://github.com/secrierlab/enrichmap",
    "use_repository_button": True,
    "use_issues_button": True,
    "logo": {"image_dark": "_static/logo-dark.svg"},
}

# Other output formats

epub_show_urls = "footnote"
htmlhelp_basename = f"{project}doc"

_doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", _doc_title, author, "manual")]
man_pages = [(master_doc, project, _doc_title, [author], 1)]

# App setup


def setup(app: Sphinx):
    """Register custom roles."""
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))
