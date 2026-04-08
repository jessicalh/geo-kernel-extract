project = "nmr_extract"
author = "Jessica"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
]

autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstrings = True
napoleon_numpy_docstrings = False
add_module_names = False

html_theme = "sphinx_rtd_theme"
html_static_path = []
