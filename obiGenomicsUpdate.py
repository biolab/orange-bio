from __future__ import with_statement
import os
import shelve
import datetime
import tarfile
import pickle
import orngServerFiles
from UserDict import DictMixin

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
        
class Singleton(object):
    __single_instances = {}            
    def __init__(self, *args):
        if self.__class__ in getattr(self, "_Singleton__single_instances"):
            raise Exception("Creating a singleton class, use getinstance() instead")
        else:
            getattr(self, "_Singleton__single_instances")[self.__class__] = self
    @classmethod
    def getinstance(cls, *args, **kw):
        instances = getattr(cls, "_Singleton__single_instances")
        if cls in instances:
            return instances[cls]
        else:
            return cls(*args, **kw)
            
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
        self.local_database_path = local_database_path
        self.progressCallback = progressCallback
        import os
        try:
            os.mkdir(self.local_database_path)
        except:
            pass
        self.shelve = UpdateShelve(os.path.join(self.local_database_path, ".updates.shelve"))

    def IsUpdatable(self, *args):
        raise NotImplementedError
    
    def GetLocal(self):
        return self.shelve.keys()
    
    def GetUpdatable(self):
        return [item for item in self.GetLocal() if self.IsUpdatable(*item)]
        
    def GetDownloadable(self):
        raise NotImplementedError
    
    def GetLastUpdateTime(self, callable, args):
        return self.shelve.get((callable, args), datetime(1,1,1))

    def _update(self, callable, args, date=None):
        self.shelve[callable, args] = date if date else datetime.now()
        self.shelve.sync()

import orngServerFiles

class PKGUpdate(Update):
    def __init__(self, domain, wrappedUpdater, *args, **kwargs):
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
        if func.im_class == type(self):
            return func(self, *args)
        else:
            return getattr(self, func.__name__)(*args)
            
    def __getattr__(self, name):
        try:
            return self.UpdateWrapper(getattr(self.wrappedUpdaterClass, name))
        except AttributeError:
            raise AttributeError(name)
        
class PKGManager(object):
    def __init__(self, updater, tarballs=[], compression="gz", serverFiles=None, domain=None):
        self.local_database_path = os.path.normpath(updater.local_database_path)
        self.updater = updater
        self.tarballs = [os.path.normpath(tar) for tar in tarballs]
        self.compression = compression or ""
        self.serverFiles = serverFiles
        self.domain = domain
        self.initialState = {}
        self.endState = {}

    def Update(self):
        for func, args in self.updater.GetUpdatable() + self.updater.GetDownloadable():
            try:
                self.MakePKG(func, args)
            except Exception, ex:
                print "Update package failed due to:", ex

    def MakePKG(self, func, args):
        realpath = os.path.realpath(os.curdir)
        os.chdir(self.local_database_path)
        self.InitWatch()
        func(self.updater, *args)
        self.EndWatch()
        name = self.Create(func, args)
        if name and self.serverFiles:
            try:
                self.serverFiles.create_domain(self.domain)
            except:
                pass
            print "Uploading %s to server" % name
            self.serverFiles.upload(self.domain, name, open(name, "rb"), title=name, tags=[func.im_class.__module__, func.__name__] + [str(arg) for arg in args] )
            self.serverFiles.unprotect(self.domain, name)
        os.chdir(realpath)
        
    def Collect(self, state):
        for dirpath, dirnames, filenames in os.walk("."):
            for filename in filenames:
                if not filename.startswith("."):
                    time = os.stat(os.path.join(dirpath, filename)).st_mtime
                    state[os.path.normpath(os.path.join(dirpath, filename))] = time
                
    def InitWatch(self):
        self.initialState = {}
        self.endState = {}
        self.Collect(self.initialState)

    def EndWatch(self):
        self.endState = {}
        self.Collect(self.endState)

    def Diff(self):
        return list(set(self.endState.items()) - set(self.initialState.items()))

    def Create(self, func, args):
        """Create the package and return its name otherwise return None"""
        name = func.__name__ + ("_" + str(args) if args else "")
        tarDirs = set()
        files = set()
        for filename, time in self.Diff():
            if self.InTarball(filename):
                tarDirs.add(self.InTarball(filename))
            else:
                files.add(filename)
        if len(files)+len(tarDirs) > 0:
            name = name+".tar" + (("." + self.compression) if self.compression else "")
            tarFile = tarfile.open(name, "w:"+self.compression)
            print "Creating:",  name
            for file in files:
                tarFile.add(file)
            for tarDir in tarDirs:
                tarFile.add(tarDir)
            tarFile.close()
            return name
        
    def InTarball(self, filename):
        for tarball in self.tarballs:
            if filename.startswith(tarball):
                return tarball
