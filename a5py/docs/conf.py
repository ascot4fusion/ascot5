# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Set path to where the Python source can be found---------------------------
import os
import sys
sys.path.insert(0, os.path.abspath("../../a5py"))

# -- Project information -------------------------------------------------------
project   = "ASCOT5"
copyright = "2023, Ascot Group"
author    = "Ascot Group"
release   = "5.4"

# -- General configuration -----------------------------------------------------
extensions = ["sphinx.ext.autodoc",     # For generating docs from Python source
              "numpydoc",               # Source docs are done in numpy style
              "nbsphinx",               # Embed Jupyter notebooks
              "breathe",]               # Sphinx can access Doxygen output

exclude_patterns = []

numpydoc_show_class_members = False # Removes table summarizing class methods

# -- Where Doxygen generated xml files are located -----------------------------
breathe_projects = {"ascot5": "doc/xml"}
breathe_default_project = "ascot5"

# -- Options for HTML output ---------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": True,
    }
html_logo = "logo.png"

html_static_path = ['_static']
html_css_files = [
    'custom.css',
]
