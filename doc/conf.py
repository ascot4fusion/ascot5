# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Set path to where the Python source can be found---------------------------
import os
import sys
sys.path.insert(0, os.path.abspath("../"))
sys.path.append(os.path.abspath("_ext"))

from multiproject.utils import get_project

import a5py

extensions = [
    "breathe",                   # Sphinx can access Doxygen output
    "nbsphinx",                  # Embed Jupyter notebooks
    "multiproject",              # For making separate user and dev docs
    "sphinx_design",             # Tabs and other nice widgets
    "options_table",             # Turns ASCOT5 options to RST table (home-made).
    "sphinx.ext.autodoc",        # For generating doc from Python source
    "numpydoc",                  # Source docs are done in numpy style (must be loaded after autodoc)
    "sphinxcontrib.bibtex",      # Can use bibtex
    "sphinxcontrib.mermaid",     # Graphs and diagrams
    "sphinx_design_elements",    # Additional nice widgets such as tables
    "sphinx.ext.autosummary",    # Creating summary tables
    "sphinx.ext.intersphinx",    # Link to external libraries
    "sphinx_gallery.load_style", # Examples thumbnails
]

multiproject_projects = {
   "user": {
       "use_config_file": False,
       },
   "dev": {
       "use_config_file": False,
   },
}

current_project = get_project(multiproject_projects)

locale_dirs = [f"{current_project}"]
if current_project == "user":
    project = "User Documentation"
elif current_project == "dev":
    project = "Development documentation"

# -- Project information -------------------------------------------------------
author = "Ascot Group"
project = "ASCOT5"
release = "6.0.0"#a5py.data.access.leaf.VERSION
copyright = "2023, Ascot Group"

# -- General configuration -----------------------------------------------------
exclude_patterns = []
autoclass_content = "both"
autodoc_typehints = "none"
autosummary_generate = True
add_module_names = False            # Removes "a5py..." path from class names
autodoc_member_order = "bysource"   # Autodoc lists things in same order as in source
numpydoc_xref_param_type = True     # Automatically link str, array_like, etc.
numpydoc_show_class_members = False # Removes table summarizing class methods

# -- Where Doxygen generated xml files are located -----------------------------
breathe_default_project = "ascot5"
breathe_projects = {"ascot5": "_static/doxygen/xml"}
intersphinx_mapping = {
    "h5py": ("https://docs.h5py.org/en/stable/", None),
    "unyt": ("https://unyt.readthedocs.io/en/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "python": ("https://docs.python.org/3", None),
    "typing": ("https://typing-extensions.readthedocs.io/en/latest/", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    }

# -- Options for HTML output ---------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": True,
    "navigation_depth": 6,
    }
html_logo = "_static/logo.png"

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Bibtex
bibtex_bibfiles = ["ascotwork.bib"]
nbsphinx_execute = "never"
nbsphinx_requirejs_path = "''"

mermaid_version = "10.9.0"

# Thumbnails (also link images somewhere as otherwise they are not copied to
# _images) (deprecated, now possible to display from the notebook)
nbsphinx_thumbnails = {
    "tutorials/introduction"  : "_images/iconlarge.png",
    "tutorials/poincare"      : "_images/poincare.png",
    "tutorials/slowingdown"   : "_images/slowingdown.png",
    "tutorials/wallload"      : "_images/wallload.png",
    "tutorials/mhd"           : "_images/mhd.png",
    "tutorials/orbits"        : "_images/orbits.png",
    "tutorials/biosaw"        : "_images/biosaw.png",
    "tutorials/distributions" : "_images/distributions.png",
    "tutorials/reversetime"   : "_images/reversetime.png",
}
