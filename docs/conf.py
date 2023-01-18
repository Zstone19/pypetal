
import os
import sys
sys.path.insert( 0, os.path.abspath('../..') )
sys.path.insert( 1, os.path.abspath("../src/"))
sys.path.insert( 2, os.path.abspath("../src/pypetal/"))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pypetal'
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
    'sphinx_copybutton',
    'sphinxcontrib.bibtex'
]
suppress_warnings = ["nbsphinx"]
master_doc = 'index'

templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

bibtex_bibfiles = ['refs.bib']
bibtex_default_style = 'plain'
bibtex_reference_style = 'author_year'



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['css/custom.css']


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
    github_repo="pypetal",
    conf_py_path="/docs/",
)

#nbsphinx config
nbsphinx_kernel_name = 'python3'
nbsphinx_execute = "never"


# -- Options for mock output ---------------------------------------------

import sys
from unittest.mock import MagicMock

class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()

MOCK_MODULES = ['javelin', 'javelin.zylc', 'javelin.lcio', 'javelin.lcmodel']
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)
