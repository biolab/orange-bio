import os
import shelve
import datetime
import tarfile

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
        def f(*args, **kw):
            lock.acquire()
            ret = func(*args, **kw)
            lock.release()
            return ret
        return f
    return syncfunc

class Update(Singleton):
    def __init__(self, local_database_path, progressCallback=None):
        Singleton.__init__(self)
        self.local_database_path = local_database_path
        self.progressCallback = progressCallback
        self.shelve = shelve.open(os.path.join(self.local_database_path, ".updates.shelve"))

    def GetUpdatable(self):
        raise NotImplementedError

    def GetDownloadable(self):
        raise NotImplementedError
    
    def GetLastUpdateTime(self, callable, args):
        return self.shelve.get(str((callable, args)), datetime.datetime(1900,1,1))

    def _update(self, callable, args, date=None):
        self.shelve[str((callable, args))] = date if date is None else datetime.datetime.now()
        self.shelve.sync()


class PKGUpdate(Update):
    def __init__(self, server_address, *args, **kwargs):
        Update.__init__(self, *args, **kwargs)
        self.server_address = server_address
    def GetUpdatable(self):
        pass
    def GetDownloadable(self):
        pass

class PKGManager(object):
    def __init__(self, updater, tarballs=[], compression="gz"):
        self.local_database_path = os.path.normpath(updater.local_database_path)
        self.updater = updater
        self.tarballs = [os.path.normpath(tar) for tar in tarballs]
        self.compression = compression or ""
        self.initialState = {}
        self.endState = {}

    def Update(self):
        os.chdir(self.local_database_path)
        for func, desc, argList in self.updater.GetUpdatable() + self.updater.GetDownloadable():
            for args in argList:
                self.InitWatch()
                print desc, "Calling:", func, "with:", args
                func(self.updater,args)
                self.EndWatch()
                self.Create(func.__name__+str(args))

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

    def Create(self, name):
        initial = set(self.initialState.items())
        end = set(self.endState.items())
        update = end - initial
        tarDirs = set()
        files = set()
        for filename, time in list(update):
            if self.InTarball(filename):
                tarDirs.add(self.InTarball(filename))
            else:
                files.add(filename)
        if len(files)+len(tarDirs) > 0:
            tarFile = tarfile.open(name+".tar" + (("." + self.compression) if self.compression else ""), "w:"+self.compression)
            print "Creating: " + os.path.join(self.local_database_path, name+".tar" + (("." + self.compression) if self.compression else ""))
            for file in files:
                tarFile.add(file)
            for tarDir in tarDirs:
                tarFile.add(tarDir)
        
    def InTarball(self, filename):
        for tarball in self.tarballs:
            if filename.startswith(tarball):
                return tarball
