"""ServerFiles"""
from __future__ import absolute_import

import serverfiles
import os

try:
    from Orange.utils import environ
except ImportError:
    from . import environ


_server_url = "http://orange.biolab.si/serverfiles-bio/"


class ServerFiles(serverfiles.ServerFiles):

    def __init__(self, server=_server_url):
        serverfiles.ServerFiles.__init__(self, server)


PATH = os.path.join(environ.buffer_dir, "serverfiles-bio")
LOCALFILES = serverfiles.LocalFiles(PATH, serverfiles=ServerFiles())


def localpath(*args):
    return LOCALFILES.localpath(*args)

    
def listfiles(*args):
    return [fname for domain, fname in LOCALFILES.listfiles(*args)]


def localpath_download(*path, **kwargs):
    return LOCALFILES.localpath_download(*path, **kwargs)


def download(*path):
    return LOCALFILES.download(*path)


def info(*path):
    return LOCALFILES.info(*path)


def update(*path, **kwargs):
    return LOCALFILES.update(*path, **kwargs)


def sizeformat(size):
    return serverfiles.sizeformat(size)

