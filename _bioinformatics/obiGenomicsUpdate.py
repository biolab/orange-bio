from __future__ import absolute_import, with_statement

import os
import shelve
import datetime
import tarfile
import pickle
from UserDict import DictMixin

from Orange.orng import orngServerFiles

class UpdateShelve(DictMixin):
    def __init__(self, path):
        self.__shelve = shelve.open(path)

    def __expand(self, (func, args)):
        return func.im_class, func.__name__, args

    def __collapse(self, (class_, funcname, args)):
        return getattr(class_, funcname), args

    def __dumps(self, arg):
        return pickle.dumps(self.__expand(arg))

    def __loads(self, arg):
        return self.__collapse(pickle.loads(arg))
        
    def __getitem__(self, arg):
        return self.__shelve.__getitem__(self.__dumps(arg))

    def __setitem__(self, arg, val):
        self.__shelve.__setitem__(self.__dumps(arg), val)

    def __delitem__(self, arg):
        self.__shelve.__delitem__(self.__dumps(arg))

    def keys(self):
        return [self.__loads(key) for key in self.__shelve.keys()]

    def __contains__(self, arg):
        return self.__shelve.__contains__(self.__dumps(arg))

    def __iter__(self):
        for key in self.__shelve:
            yield self.__loads(key)

    def iteritems(self):
        for key, val in self.__shelve.iteritems():
            yield self.__loads(key), val

    def sync(self):
        self.__shelve.sync()
            
from threading import Lock
from functools import wraps

def synchronized(lock):
    def syncfunc(func):
        @wraps(func)
        def f(*args, **kwargs):
            with lock:
                return func(*args, **kwargs)
        return f
    return syncfunc

from datetime import datetime

class Update(object):
    def __init__(self, local_database_path, progressCallback=None):
        """Base update object.
        Each obi module that uses updatable resources should subclass this class and at minimum reimplement
        the IsUpdatable and GetDownloadable methods.
        """
        self.local_database_path = local_database_path
        self.progressCallback = progressCallback
        import os
        try:
            os.mkdir(self.local_database_path)
        except:
            pass
        self.shelve = UpdateShelve(os.path.join(self.local_database_path, ".updates.shelve"))

    def IsUpdatable(self, *args):
        """Return True if the local can be updated else return False
        """
        raise NotImplementedError
    
    def GetLocal(self):
        """Return a list [(callable, args), ...] that have been called in the past.
        """
        return self.shelve.keys()
    
    def GetUpdatable(self):
        """Return a list [(callable, args), ...] that can be updated.
        """
        return [item for item in self.GetLocal() if self.IsUpdatable(*item)]
        
    def GetDownloadable(self):
        """Return a list [(callable, args), ...] that can be downloaded. Must no contain any
        (callable, args) that are returnd by GetLocal.
        """
        raise NotImplementedError
    
    def GetLastUpdateTime(self, callable, args):
        """Returns the last update time of callable with args.
        """
        return self.shelve.get((callable, args), datetime(1,1,1))

    def _update(self, callable, args, date=None):
        """Insert the record of update with date. if date is None use the current date.
        """
        self.shelve[callable, args] = date if date else datetime.now()
        self.shelve.sync()

from Orange.orng import orngServerFiles

class PKGUpdate(Update):
    def __init__(self, domain, wrappedUpdater, *args, **kwargs):
        """Wrap the subclass of Update to download data as packages from our own server.
        Arguments:
            - domain : should be a string that identifies a resource e.g. "kegg", "go" ...
            - wrappedUpdater : instance of a subclass of Update that is wrapped
        Example:
        >>> pkg = PKGUpdate("go", obiGO.Update())
        >>> for item in pkg.GetUpdatable():
        ...     pkg.Apply(*item)
        ...
        >>>
        """
        Update.__init__(self, wrappedUpdater.local_database_path, wrappedUpdater.progressCallback)
        self.domain = domain
        self.serverFiles = orngServerFiles.ServerFiles()
        self.wrappedUpdater = wrappedUpdater
        self.wrappedUpdaterClass = wrappedUpdater.__class__

    def _GetServerFileName(self, func, args):
        return func.__name__+ ("_" + str(args) if args else "") + ".tar.gz"

    def _GetServerDateTime(self, *args):
        time = self.serverFiles.info(self.domain, self._GetServerFileName(*args))["datetime"].split(".")[0]
        return datetime.strptime(time, "%Y-%m-%d %H:%M:%S")
    
    def IsUpdatable(self, *args):
        try:
            return args in self.shelve and self.GetLastUpdateTime(*args) < self._GetServerDateTime(*args)
        except Exception:
            return False
        
    def GetDownloadable(self):
        files = self.serverFiles.listfiles(self.domain)
        downloadable = []
        for filename in files:
            if "_" in filename:
                funcname, args = filename.split("_",1)
                args = args.split(".")[0]
                args = eval(args) 
            else:
                funcname, args = filename.split(".")[0], ()
            downloadable.append((getattr(self.wrappedUpdaterClass, funcname), args))
        return [item for item in downloadable if item not in self.shelve]

    def UpdateWrapper(self, func):
        @wraps(func)
        def update(*args):
            filename = self._GetServerFileName(func, args)
            time = self._GetServerDateTime(func, args)
            orngServerFiles.download(self.domain, filename)
            filepath = orngServerFiles.localpath(self.domain, filename)
            tar = tarfile.open(filepath)
            tar.extractall(self.local_database_path)
            self._update(func, args, time)
        return update

    def Apply(self, func, args):
        """Use this method to apply a return values of GetUpdatable ... which are unwrapped
        """
        if func.im_class == type(self):
            return func(self, *args)
        else:
            return getattr(self, func.__name__)(*args)
            
    def __getattr__(self, name):
        try:
            return self.UpdateWrapper(getattr(self.wrappedUpdaterClass, name))
        except AttributeError:
            raise AttributeError(name)

def firstUpdateConsole():
    from . import obiKEGG, obiGO, obiGenomicsUpdate

    pkgUpdate = obiGenomicsUpdate.PKGUpdate("go", obiGO.Update())

    print "Updating GO ontology"
    pkgUpdate.UpdateOntology()

    for org, name in [("goa_human", "Homo sapiens"), ("sgd", "Yeast")]:
        print "Updating GO anotations for", name
        pkgUpdate.UpdateAnnotation(org)

def firstUpdateQt():
    pass
    
