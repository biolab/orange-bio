import shelve
import datetime
import os

class Update(object):
    def __init__(self, local_database_path, progressCallback=None):
        self.local_database_path = local_database_path
        self.progressCallback = progressCallback
        self.shelve = shelve.open(os.path.join(self.local_database_path, "updates.shelve"))

    def GetUpdatable(self):
        raise Exception("Not implemented")

    def GetDownloadable(self):
        raise Exception("Not implemented")
    
    def GetLastUpdateTime(self, callable, args):
        return self.shelve.get((str(callable), args), None)

    def _update(self, callable, args):
        self.shelve[(str(callable), args)] = datetime.datetime.now()
        
    
    