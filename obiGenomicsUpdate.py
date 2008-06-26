import os
import shelve
import datetime

class Singleton(object):
    __single_instances = {}
    def __init__(self):
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
        self.shelve = shelve.open(os.path.join(self.local_database_path, "updates.shelve"))

    def GetUpdatable(self):
        raise Exception("Not implemented")

    def GetDownloadable(self):
        raise Exception("Not implemented")
    
    def GetLastUpdateTime(self, callable, args):
        return self.shelve.get(str((callable, args)), datetime.datetime(1900,1,1))

    def _update(self, callable, args):
        self.shelve[str((callable, args))] = datetime.datetime.now()
        self.shelve.sync()
        
    
    