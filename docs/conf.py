# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from sphinx.ext.autosummary import Autosummary
from sphinx.ext.autosummary import get_documenter
from docutils.parsers.rst import directives
from sphinx.util.inspect import safe_getattr
import os
import sys
import sphinx
import sphinx_rtd_theme
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = u'HPS'
copyright = u'2021, Quyen Vu'
author = u'Quyen Vu'

# The short X.Y version
version = '2021.0'
# The full version, including alpha/beta/rc tags
release = 'v1.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'numpydoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage'
]

autosummary_generate = True
#numpydoc_show_class_members = False
autodoc_default_flags = ['members', 'inherited-members']


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
# If true, the current module name will be prepended to all description unit titles.
add_module_names = True
master_doc = 'index'
html_sidebars = {
        '**': ['localtoc.html', 'sourcelink.html', 'searchbox.html'],
}

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': 'white',
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': True
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = "_static/logo.svg"


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/': None}


class AutoAutoSummary(Autosummary):

    option_spec = {
        'methods': directives.unchanged,
        'attributes': directives.unchanged
    }

    required_arguments = 1

    @staticmethod
    def get_members(obj, typ, include_public=None):
        if not include_public:
            include_public = []
        items = []
        for name in dir(obj):
            try:
                documenter = get_documenter(safe_getattr(obj, name), obj)
            except AttributeError:
                continue
            if documenter.objtype == typ:
                items.append(name)
        public = [x for x in items if x in include_public or not x.startswith('_')]
        return public, items

    def run(self):
        clazz = str(self.arguments[0])
        try:
            (module_name, class_name) = clazz.rsplit('.', 1)
            m = __import__(module_name, globals(), locals(), [class_name])
            c = getattr(m, class_name)
            if 'methods' in self.options:
                _, methods = self.get_members(c, 'method', ['__init__'])

                self.content = ["~%s.%s" % (clazz, method) for method in methods if not method.startswith('_')]
            if 'attributes' in self.options:
                _, attribs = self.get_members(c, 'attribute')
                self.content = ["~%s.%s" % (clazz, attrib) for attrib in attribs if not attrib.startswith('_')]
        finally:
            return super(AutoAutoSummary, self).run()


def setup(app):
    app.add_directive('autoautosummary', AutoAutoSummary)


############################
# SETUP THE RTD LOWER-LEFT #
############################
try:
    html_context
except NameError:
    html_context = dict()
html_context['display_lower_left'] = True

templates_path = ['_templates']

if 'REPO_NAME' in os.environ:
    REPO_NAME = os.environ['REPO_NAME']
else:
    REPO_NAME = 'hpsOpenMM'

# SET CURRENT_LANGUAGE
if 'current_language' in os.environ:
    # get the current_language env var set by buildDocs.sh
    current_language = os.environ['current_language']
else:
    # the user is probably doing `make html`
    # set this build's current language to english
    current_language = 'en'

# tell the theme which language to we're currently building
html_context['current_language'] = current_language

# SET CURRENT_VERSION
from git import Repo

repo = Repo(search_parent_directories=True)

if 'current_version' in os.environ:
    # get the current_version env var set by buildDocs.sh
    current_version = os.environ['current_version']
else:
    # the user is probably doing `make html`
    # set this build's current version by looking at the branch
    current_version = repo.active_branch.name

# tell the theme which version we're currently on ('current_version' affects
# the lower-left rtd menu and 'version' affects the logo-area version)
html_context['current_version'] = current_version
html_context['version'] = current_version

# POPULATE LINKS TO OTHER LANGUAGES
html_context['languages'] = [('en', '/' + REPO_NAME + '/en/' + current_version + '/')]


# POPULATE LINKS TO OTHER VERSIONS
html_context['versions'] = list()

versions = [branch.name for branch in repo.branches]
for version in versions:
    html_context['versions'].append((version, '/' + REPO_NAME + '/' + current_language + '/' + version + '/'))

# POPULATE LINKS TO OTHER FORMATS/DOWNLOADS

# settings for creating PDF with rinoh
rinoh_documents = [(
    master_doc,
    'target',
    project + ' Documentation',
    '© ' + copyright,
)]
today_fmt = "%B %d, %Y"

##########################
# "EDIT ON GITHUB" LINKS #
##########################

html_context['display_github'] = True
html_context['github_user'] = 'qvv5013'
html_context['github_repo'] = 'rtd-github-pages'
html_context['github_version'] = 'main/docs/'
