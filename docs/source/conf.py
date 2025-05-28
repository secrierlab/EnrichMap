"""Configuration for EnrichMap's Sphinx documentation."""

from __future__ import annotations
from pathlib import Path
from functools import partial
from docutils import nodes
from typing import TYPE_CHECKING

HERE = Path(__file__).parent

import enrichmap

if TYPE_CHECKING:
    from sphinx.application import Sphinx

# -- Project information
project = "EnrichMap"
copyright = "Cenk Celik"
author = "Cenk Celik"

release = "0.1"
version = "0.1.5"

master_doc = "index"

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx_book_theme",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

# -- Options for EPUB output

epub_show_urls = "footnote"

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
api_dir = HERE / "api"
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto", "ftp")
myst_heading_anchors = 3
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True


ogp_site_url = "https://enrichmap.readthedocs.io/en/stable/"
ogp_image = "https://github.com/secrierlab/enrichmap/raw/main/img/enrichmap_logo.svg"

typehints_defaults = "braces"

pygments_style = "default"
pygments_dark_style = "native"


# -- Options for HTML output ----------------------------------------------
html_theme = "sphinx_book_theme"
# The theme is sphinx-book-theme, with patches for readthedocs-sphinx-search
html_theme_options = {
    "repository_url": "https://github.com/secrierlab/enrichmap",
    "use_repository_button": True,
    "use_issues_button": True,
}
html_show_sphinx = False
html_logo = "https://github.com/secrierlab/enrichmap/raw/main/img/enrichmap_logo.svg"
html_title = "EnrichMap"


def setup(app: Sphinx):
    """App setup hook."""
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))
    app.add_config_value(
        "recommonmark_config",
        {
            "auto_toc_tree_section": "Contents",
            "enable_auto_toc_tree": True,
            "enable_math": True,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,  # noqa: FBT003
    )


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
