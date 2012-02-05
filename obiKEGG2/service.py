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
    client = Client(KEGG_WDSL)
    return client.service

def suds_service_with_requests():
    import StringIO
    
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
            message = request.message
            headers = request.headers
            # This does not seem to work, every request seems
            # to open a new connection
            headers["Connection"] = "Keep-Alive"
            try:
                response = requests.post(url, data=message, 
                                         headers=headers,
                                         cookies=self.cookies)
                    
                self.proxy = self.options.proxy
                response.raise_for_status()
                
                self.cookies.update(response.cookies)
                result = RequestsResponse()
                result.code = response.status_code
                result.headers = response.headers
                result.message = response.raw.read()
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
                response = requests.get(url)
                self.proxy = self.options.proxy
                
                response.raise_for_status()
                return StringIO.StringIO(response.raw.read())
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


