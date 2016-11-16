#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
import sys
import os
import shlex
import sphinx_rtd_theme


extensions = [
    'sphinx.ext.mathjax',
]


templates_path = ['_templates']
source_suffix = '.rst'
source_encoding = 'utf-8-sig'
master_doc = 'index'

project = 'gams-matlab'
copyright = '2016, ojdo'
author = 'ojdo'
version = '0.1'
release = '0.1'

language = None
today_fmt = '%B %d, %Y'

exclude_patterns = ['_build']
#add_function_parentheses = True
#add_module_names = True
pygments_style = 'sphinx'
todo_include_todos = False


# html

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']

htmlhelp_basename = 'gams-matlabdoc'


# latex

latex_elements = { }
latex_documents = [
  (master_doc, 'gams-matlab.tex', 'gams-matlab Documentation', 
   'ojdo', 'manual'),
]

# manpages

man_pages = [
    (master_doc, 'gams-matlab', 'gams-matlab Documentation',
     [author], 1)
]


# texinfo
texinfo_documents = [
  (master_doc, 'gams-matlab', 'gams-matlab Documentation',
   author, 'gams-matlab', 'One line description of project.',
   'Miscellaneous'),
]
