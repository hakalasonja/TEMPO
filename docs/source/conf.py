# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TEMPO'
copyright = '2024, Jner Tzern Oon, Sonja Hakala, George Witt, Ronald Walsworth'
author = 'Jner Tzern Oon, Sonja Hakala, George Witt, Ronald Walsworth'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'myst_parser',  
]

# Allow Sphinx to parse Markdown files
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
# Set the theme to Read the Docs theme (installed via pip)
html_theme = 'sphinx_rtd_theme'

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
  "external_links": [
      {"name": "Home", "url": "https://tempo-documentation.readthedocs.io/en/latest/index.html"},
  ]
}