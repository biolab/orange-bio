"""ServerFiles"""
from __future__ import absolute_import

import serverfiles
import os

try:
    from Orange.utils import environ
except ImportError:
    from . import environ

_server_url = "http://193.2.72.57/newsf/"


class ServerFiles(serverfiles.ServerFiles):

    def __init__(self, server=_server_url):
        serverfiles.ServerFiles.__init__(self, server)

_path = os.path.join(environ.buffer_dir, "testServerFiles")
_localFiles = serverfiles.LocalFiles(_path, serverfiles=ServerFiles())


def localpath(*args):
    return _localFiles.localpath(*args)

    
def listfiles(*args):
    return [fname for domain, fname in _localFiles.listfiles(*args)]


def localpath_download(*path, **kwargs):
    return _localFiles.localpath_download(*path, **kwargs)


def download(*path):
    return _localFiles.download(*path)


def info(*path):
    return _localFiles.info(*path)


def update(*path, **kwargs):
    return _localFiles.update(*path, **kwargs)


def sizeformat(size):
    return serverfiles.sizeformat(size)
