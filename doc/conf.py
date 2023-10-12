# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Set path to where the Python source can be found---------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

import a5py

# -- Project information -------------------------------------------------------
project   = 'ASCOT5'
copyright = '2023, Ascot Group'
author    = 'Ascot Group'
release   = a5py.ascot5io.coreio.fileapi.VERSION

# -- General configuration -----------------------------------------------------
extensions = ['sphinx.ext.autodoc',      # For generating doc from Python source
              'numpydoc',                # Source docs are done in numpy style
              'nbsphinx',                # Embed Jupyter notebooks
              'breathe',                 # Sphinx can access Doxygen output
              'sphinxcontrib.bibtex',    # Can use bibtex
              'sphinx.ext.autosummary',  # Creating summary tables
              'sphinx.ext.intersphinx',  # Link to external libraries
              "sphinx_gallery.load_style",] # Examples thumbnails

exclude_patterns = []
numpydoc_xref_param_type = True # Automatically link str, array_like, etc.
numpydoc_show_class_members = False # Removes table summarizing class methods

# -- Where Doxygen generated xml files are located -----------------------------
breathe_projects = {'ascot5': '_static/doxygen/xml'}
breathe_default_project = 'ascot5'

intersphinx_mapping = {
    'python'    : ('https://docs.python.org/3', None),
    'numpy'     : ('https://docs.scipy.org/doc/numpy', None),
    'scipy'     : ('https://docs.scipy.org/doc/scipy/reference', None),
    'matplotlib': ('https://matplotlib.org/stable', None)}

# -- Options for HTML output ---------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'navigation_depth':6,
    }
html_logo = '_static/logo.png'

html_static_path = ['_static']
html_css_files = [
    'custom.css',
]

# Bibtex
bibtex_bibfiles = ['ascotwork.bib']
nbsphinx_execute = 'never'

# Thumbnails (also link images somewhere as otherwise they are not copied to
# _images)
nbsphinx_thumbnails = {
    "tutorials/introduction" : "_images/iconlarge.png",
    "tutorials/poincare"     : "_images/poincare.png",
    "tutorials/slowingdown"  : "_images/slowingdown.png",
    "tutorials/wallload"     : "_images/wallload.png",
    "tutorials/mhd"          : "_images/mhd.png",
    "tutorials/orbits"       : "_images/orbits.png",
    "tutorials/biosaw"       : "_images/biosaw.png",
}
