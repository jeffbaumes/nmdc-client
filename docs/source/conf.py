project = "nmdc-client"
copyright = "2024, NMDC Collaborators"
author = "NMDC Collaborators"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_mdinclude",
]

templates_path = ["_templates"]
exclude_patterns = []

# html_theme = 'sphinx_rtd_theme'
html_theme = "furo"
html_static_path = ["_static"]

# Add the favicon
html_favicon = "_static/favicon.ico"
