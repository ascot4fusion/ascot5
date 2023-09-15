# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Set path to where the Python source can be found---------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../../a5py'))

import a5py

from inspect import getsourcefile

# Get path to directory containing this file, conf.py.
DOCS_DIRECTORY = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))

def ensure_pandoc_installed(_):
    import pypandoc

    # Download pandoc if necessary. If pandoc is already installed and on
    # the PATH, the installed version will be used. Otherwise, we will
    # download a copy of pandoc into docs/bin/ and add that to our PATH.
    pandoc_dir = os.path.join(DOCS_DIRECTORY, "bin")
    # Add dir containing pandoc binary to the PATH environment variable
    if pandoc_dir not in os.environ["PATH"].split(os.pathsep):
        os.environ["PATH"] += os.pathsep + pandoc_dir
    pypandoc.ensure_pandoc_installed(
        quiet=True,
        targetfolder=pandoc_dir,
        delete_installer=True,
    )


def setup(app):
    app.connect("builder-inited", ensure_pandoc_installed)

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
              'sphinx.ext.intersphinx',] # Link to external libraries

exclude_patterns = []
numpydoc_xref_param_type = True # Automatically link str, array_like, etc.
numpydoc_show_class_members = False # Removes table summarizing class methods

# -- Where Doxygen generated xml files are located -----------------------------
breathe_projects = {'ascot5': '_static/xml'}
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
html_logo = 'figs/logo.png'

html_static_path = ['_static']
html_css_files = [
    'custom.css',
]

# Bibtex
bibtex_bibfiles = ['ascotwork.bib']
nbsphinx_execute = 'never'
