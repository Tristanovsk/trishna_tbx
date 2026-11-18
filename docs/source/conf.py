# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import shutil
import sys
from pathlib import Path

# Make the package importable without installation (local builds);
# on Read the Docs the package is also pip-installed (see .readthedocs.yaml).
DOCS_SOURCE = Path(__file__).resolve().parent
REPO_ROOT = DOCS_SOURCE.parents[1]
sys.path.insert(0, str(REPO_ROOT))

import trishna_tbx

# -- Project information -----------------------------------------------------

project = 'trishna_tbx'
copyright = '2025-2026, CNES/Magellium'
author = 'Tristan Harmel'
version = trishna_tbx.__version__
release = trishna_tbx.__version__
today_fmt = '%Y-%m-%d'

# -- Tutorials ---------------------------------------------------------------

# The tutorials are the notebooks of the VRE gallery (notebooks/<section>/*.ipynb, written by
# notebooks/build_notebooks.py). They are copied here at each build so that there is only one
# version of them in the repository.
TUTORIALS = DOCS_SOURCE / 'tutorials'
for nb in sorted((REPO_ROOT / 'notebooks').glob('[0-9]_*/*.ipynb')):
    dest = TUTORIALS / nb.parent.name / nb.name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or dest.stat().st_mtime < nb.stat().st_mtime:
        shutil.copy2(nb, dest)

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',
    'sphinxcontrib.mermaid',
    'myst_nb',
]

templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', '**.ipynb_checkpoints', 'Thumbs.db', '.DS_Store']

# -- Autodoc / autosummary ---------------------------------------------------

autosummary_generate = True
autoclass_content = 'both'          # the parameters are documented in __init__
autodoc_typehints = 'none'          # types are given in the docstrings
autodoc_member_order = 'bysource'
# 'members' is set in the autosummary templates (_templates/autosummary/) to
# avoid documenting objects twice
add_module_names = False

# The sibling packages are imported inside the functions that need them: they are
# mocked so that the API reference builds without them (Read the Docs).
autodoc_mock_imports = ['GRSdriver', 'grs', 'GRSl2bgen', 'radcalnet_oc', 'invRrs', 'hgrs']

# The docstrings use reST fields (":param x:"); napoleon also converts the
# NumPy sections ("Parameters", "Notes") if any.
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False
napoleon_use_ivar = True
napoleon_preprocess_types = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'xarray': ('https://docs.xarray.dev/en/stable', None),
    'pandas': ('https://pandas.pydata.org/docs', None),
    'radcalnet_oc': ('https://radcalnet-oc.readthedocs.io/en/latest', None),
    'grs': ('https://grs.readthedocs.io/en/latest', None),
}

# labelled equations are numbered by page and cited with :eq:
math_eqref_format = '({number})'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_book_theme'
pygments_style = 'sphinx'

html_theme_options = {
    'repository_url': 'https://github.com/Tristanovsk/trishna_tbx',
    'repository_branch': 'main',
    'path_to_docs': 'docs/source',
    'use_repository_button': True,
    'use_issues_button': True,
    'use_edit_page_button': True,
    'use_download_button': True,
    'navigation_with_keys': True,
    'show_toc_level': 2,
    'logo': {
        'image_light': '_static/trishna_tbx_logo.svg',
        'image_dark': '_static/trishna_tbx_logo_dark.svg',
    },
}

html_title = ''
html_logo = '_static/trishna_tbx_logo.svg'
html_favicon = '_static/trishna_tbx_favicon.svg'

html_static_path = ['_static']
html_css_files = ['custom.css']
html_show_sourcelink = False
html_last_updated_fmt = today_fmt

htmlhelp_basename = 'trishna_tbx_doc'

# -- Options for LaTeX output (PDF of Read the Docs) -------------------------

latex_elements = {
    # long file names and code-like words do not break: stretch the lines rather
    # than overflow the margin, and do not report the loose lines
    'preamble': r'''
\setlength{\emergencystretch}{3em}
\hbadness=10000
\vbadness=10000
\DeclareUnicodeCharacter{2212}{\ensuremath{-}}
''',
}


def setup(app):
    from docutils import nodes
    from sphinx.transforms.post_transforms import SphinxPostTransform

    class DropRemoteImages(SphinxPostTransform):
        """Remote images (the version badge of the notebooks) cannot be embedded in the PDF and EPUB:
        they are removed for the non-HTML builders, before Sphinx downloads them (priority 100)."""
        default_priority = 50

        def run(self, **kwargs):
            if self.app.builder.name in ('html', 'dirhtml', 'singlehtml'):
                return
            for node in list(self.document.findall(nodes.image)):
                if '://' in node.get('uri', ''):
                    node.parent.remove(node)

    app.add_post_transform(DropRemoteImages)

# -- MyST / notebook rendering -----------------------------------------------

myst_enable_extensions = [
    'amsmath',
    'colon_fence',
    'deflist',
    'dollarmath',
    'html_admonition',
    'html_image',
    'linkify',
    'smartquotes',
]
myst_heading_anchors = 3

# The tutorial notebooks need Sentinel-2 images, CAMS files, HYPERNETS data and the
# look-up tables of the sibling packages, which are not available on Read the Docs:
# they are rendered with the outputs stored in the notebooks, not executed.
nb_execution_mode = 'off'
nb_merge_streams = True
suppress_warnings = ['mystnb.unknown_mime_type', 'myst.header']
