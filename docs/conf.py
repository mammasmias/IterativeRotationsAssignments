# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
import sys; sys.setrecursionlimit(1500)
# sys.path.insert(0, os.path.abspath('.'))
# sys.path.insert(0, os.path.abspath('./..'))
from sphinx.builders.html import StandaloneHTMLBuilder
import subprocess, os
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('..'))

# Doxygen
subprocess.call('doxygen Doxyfile.in', shell=True)


project = 'lib_IRA'
copyright = '2024, MAMMASMIAS Consortium'
author = 'MG'
release = 'v1.6'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
#    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
#    'sphinx.ext.autosectionlabel',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx_sitemap',
    'sphinx.ext.inheritance_diagram',
    'breathe',
    'sphinxfortran.fortran_domain',
    'sphinxfortran.fortran_autodoc',
    'sphinx_rtd_size'
]

source_suffix = {
    '.rst': 'restructuredtext',
#    '.md': 'markdown',
}

fortran_src=[os.path.abspath('../src/*.f90'), ]
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


sys.path.insert(0, os.path.abspath('../interface'))


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

sphinx_rtd_size_width = "90%"
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


# -- Breathe configuration -------------------------------------------------

breathe_projects = {
	"lib_IRA": "_build/xml/"
}
breathe_default_project = "lib_IRA"
breathe_default_members = ('members', 'undoc-members')
