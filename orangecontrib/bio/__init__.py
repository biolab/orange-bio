# namespace stub
__import__("pkg_resources").declare_namespace(__name__)

import warnings

alreadyWarned = False
disabledMsg = "Some features will be disabled due to failing modules\n"
def _import(name):
    global alreadyWarned
    try:
        __import__(name, globals(), locals(), [], -1)
    except ImportError, err:
        warnings.warn("%sImporting '%s' failed: %s" %
            (disabledMsg if not alreadyWarned else "", name, err),
            UserWarning, 2)
        alreadyWarned = True

_import("arrayexpress")
_import("biomart")
_import("gene")
_import("geo")
_import("go")
_import("mesh")
_import("taxonomy")
_import("omim")
_import("ontology")
_import("utils")

del _import
del alreadyWarned
del disabledMsg
