import shelve
import datetime

from functools import wraps

class Update(object):
    def __init__(self, local_database_path, progressCallback=None):
        self.local_database_path = local_database_path
        self.progressCallback = progressCallback
        self.shelve = shelve.open(self.local_database_path+"//updates.shelve")

    def GetUpdatable(self):
        raise Exception("Not implemented")

    def GetDownloadable(self):
        raise Exception("Not implemented")
    
    def GetLastUpdateTime(self, callable, args):
        return self.shelve.get(str((callable.func_name, args)), None)

    def _update(self, callable, args):
        self.shelve[str((callable.func_name, args))] = datetime.datetime.now()

    @classmethod
    def _auto_updater(cls, f):
        @wraps(f)
        def mywrapper(self, *args):
            r = f(self, *args)
            self._update(f, args)
            return r
        return mywrapper
        
    
    