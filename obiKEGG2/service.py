"""
SOAP service client. Needs suds library.
 
"""
from __future__ import absolute_import

KEGG_WDSL = "http://soap.genome.jp/KEGG.wsdl"

import urllib2
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

def suds_service_with_requests():
    from suds.client import Client
    from suds.transport import Transport, TransportError
    from suds.properties import Unskin
    
    import requests

    class RequestsResponse(object):
        pass

    class RequestsTransport(Transport):
        def __init__(self, **kwargs):
            Transport.__init__(self)
            Unskin(self.options).update(kwargs)
            self.cookies = {}
            
        def send(self, request):
            result = None
            url = request.url
#            print "URL", url
            message = request.message
            headers = request.headers
            headers["Connection"] = "Keep-Alive"
            try:
                response = requests.post(url, data=message, 
                                         headers=headers,
                                         cookies=self.cookies)
                    
                self.proxy = self.options.proxy
                response.raise_for_status()
                
                self.cookies.update(response.cookies)
                result = RequestsResponse()
                result.message = response.raw.read()
                result.headers = response.headers
                result.code = response.status_code
                return result
                
            except urllib2.HTTPError, e:
                if e.code in [202, 204]:
                    return None
                else:
                    raise TransportError(str(e), e.code, e.fp)
                
        def open(self, request):
            url = request.url
            message = request.message
            try:
                respose = requests.get(url)
                self.proxy = self.options.proxy
                
                respose.raise_for_status()
                return response.raw
            except urllib2.HTTPError, e:
                raise TransportError(str(e), e.code, e.fp)
            
    transport = RequestsTransport()
    client = Client(KEGG_WDSL, transport=transport)
    return client.service
    
def SOAPy_service():
    import SOAPy

from . import conf

default_service = suds_service
web_service = suds_service

if conf.params["service.transport"] == "requests":
    try:
        import requests    
        web_service = suds_service_with_requests
    except ImportError:
        import warnings
        warnings.warn("requests package not installed.")
elif conf.params["service.transport"] == "urllib2":
    pass


