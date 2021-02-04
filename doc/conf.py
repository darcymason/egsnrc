# -*- coding: utf-8 -*-
#

import sys

import sphinx_rtd_theme

import egsnrc


extensions = [
    'sphinx.ext.napoleon',  # Numpy style docstrings
    'sphinx.ext.autodoc',
    'sphinx_autodoc_typehints',
    'sphinx.ext.extlinks',
]

autosummary_generate = True

autodoc_default_options = {
    'members': None,
    'no-inherited-members': None,
}

add_module_names = False

# Shortcuts for sphinx.ext.extlinks
extlinks = {
    # 'alias' : (url_prefix, caption)
    'gh': (
        'https://github.com/egsnrc/%s',
        None
    ),
}

# intersphinx configuration
intersphinx_mapping = {
    'python': ('https://docs.python.org/{.major}'.format(
        sys.version_info), None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
}


napoleon_google_docstring = False
napoleon_numpy_docstring = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
# source_encoding = 'utf-8'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'egsnrc'
copyright = '2020-2021, Darcy Mason and egsnrc contributors'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = egsnrc.__version__  # noQA
# The full version, including alpha/beta/rc tags.
release = egsnrc.__version__

# List of directories, relative to source directory, that shouldn't be searched
# for source files.
exclude_trees = ['_build']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# Custom style
html_style = 'css/egsnrc.css'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# -- Options for HTML output -----------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
# html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Output file base name for HTML help builder.
htmlhelp_basename = 'egsnrcdoc'


# Config for sphinx_issues
issues_github_path = 'darcymason/egsnrc'


def setup(app):
    app.add_css_file('css/egsnrc.css')

# Example configuration for intersphinx: refer to
# the Python standard library.
# intersphinx_mapping = {'http://docs.python.org/': None}
