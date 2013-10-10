"""

"""
from __future__ import absolute_import

KEGG_WDSL = "http://soap.genome.jp/KEGG.wsdl"

REST_API = "http://rest.kegg.jp/"


def slumber_service():
    """
    Return a rest based service using `slumber` package
    """
    import slumber
    if not hasattr(slumber_service, "_cached"):
        slumber_service._cached = slumber.API(REST_API)
    return slumber_service._cached


from . import conf

default_service = slumber_service

web_service = slumber_service
