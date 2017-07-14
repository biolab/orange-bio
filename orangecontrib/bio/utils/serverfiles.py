"""ServerFiles"""
from __future__ import absolute_import

import serverfiles

try:
    from Orange.utils import environ
except ImportError:
    from . import environ

from orangecontrib.bio.utils import serverfile_path

server_url = "http://orange.biolab.si/serverfiles-bio/"


class ServerFiles(serverfiles.ServerFiles):

    def __init__(self, server=server_url):
        serverfiles.ServerFiles.__init__(self, server)

PATH = serverfile_path()
LOCALFILES = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())


def localpath(*args, **kwargs):
    return LOCALFILES.localpath(*args, **kwargs)

    
def listfiles(*args, **kwargs):
    return [fname for domain, fname in LOCALFILES.listfiles(*args, **kwargs)]


def localpath_download(*path, **kwargs):
    return LOCALFILES.localpath_download(*path, **kwargs)


def download(*args, **kwargs):
    return LOCALFILES.download(*args, **kwargs)


def info(*args, **kwargs):
    return LOCALFILES.info(*args, **kwargs)


def update(*path, **kwargs):
    return LOCALFILES.update(*path, **kwargs)


def sizeformat(size):
    return serverfiles.sizeformat(size)
