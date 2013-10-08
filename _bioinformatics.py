"""
Backwards compatibility stub. Should be removed by the 2.7 release.
"""
import sys
import warnings

from orangecontrib import bio

# Change the PendingDeprecation to Deprecated in 2.6 release.
warnings.warn(
    "'_bioinformatics' import name is deprecated.\n" +
    "Please use 'orangecontrib.bio'.",
    PendingDeprecationWarning,
    stacklevel=2
)

sys.modules["_bioinformatics"] = bio
