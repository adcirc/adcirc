# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ADCIRC'
# copyright = '2025, ADCIRC Developers and Users Community'
# author = 'ADCIRC Developers and Users Community'
copyright = ''
author = ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']

# Theme options
html_theme_options = {
    'navigation_depth': 2,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'titles_only': False
}

# -- Custom roles -----------------------------------------------------------
from docutils import nodes

def version_info_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    """Role for displaying version compatibility information with styling."""
    node = nodes.inline(rawtext, text, classes=['version-info'])
    return [node], []

def setup(app):
    app.add_role('version_info', version_info_role)
