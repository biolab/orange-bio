"""
Bioinformatics
==============

Bioinformatics widgets for Orange3

"""

NAME = "Bioinformatics"

ICON = "../widgets/icons/Category-Bioinformatics.svg"

BACKGROUND = "light-grass"

intersphinx = (
    ("{DEVELOP_ROOT}/docs/build/html/", None),
    ("http://pythonhosted.org/Orange-Bioinformatics/", None)
)

WIDGET_HELP_PATH = (
    # Used for development.
    # You still need to build help pages using
    # make htmlhelp
    # inside doc folder
    ("{DEVELOP_ROOT}/doc/build/htmlhelp/index.html", None),

    ("http://pythonhosted.org/Orange-Bioinformatics/", "")
)
