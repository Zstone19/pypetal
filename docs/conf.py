
import os
import sys
sys.path.insert( 0, os.path.abspath('../..') )
sys.path.insert( 1, os.path.abspath("../../src/"))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'petl'
copyright = '2022, Zachary Stone'
author = 'Zachary Stone'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'nbsphinx',
    'sphinx_copybutton'
]
suppress_warnings = ["nbsphinx"]
master_doc = 'index'

templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]





# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None
napoleon_attr_annotations = True

#GitHub config
html_context = dict(
    display_github=True,
    github_user="Zstone19",
    github_repo="petl",
    conf_py_path="/docs/",
)

#nbsphinx config
nbsphinx_kernel_name = 'python3'
nbsphinx_execute = "never"
