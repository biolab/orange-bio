"""
SOAP service client. Needs suds library.
 
"""

KEGG_WDSL = "http://soap.genome.jp/KEGG.wsdl"

def suds_service():
    """ Return an suds service object with kegg api service methods.
    
    >>> service = web_service()
    >>> service.list_databases()
    [(Definition){...

    """
    from suds.client import Client
    # TODO: extend suds.transport HttpTransport
    # to support keep-alive 
    client = Client(KEGG_WDSL)
    return client.service

def SOAPy_service():
    import SOAPy

web_service = suds_service
