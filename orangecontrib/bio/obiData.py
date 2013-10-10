from __future__ import absolute_import, with_statement

import ftplib
import urllib, urllib2
import threading
import os, sys
import time
import socket
from Queue import Queue
from StringIO import StringIO
from datetime import datetime
from threading import RLock

from .obiGenomicsUpdate import synchronized

class FileNotFoundError(IOError):
    pass

class SharedCache(object):
    def __init__(self, *args, **kwargs):
        self._dict = dict(*args, **kwargs)
        self._lock = RLock()

    def __getitem__(self, key):
        with self._lock:
            return self._dict.__getitem__(key)

    def __setitem__(self, key, item):
        with self._lock:
            return self._dict.__setitem__(key, item)

    def __delitem__(self, key):
        with self._lock:
            return self._dict.__delitem__(key)

    def __contains__(self, key):
        with self._lock:
            return self._dict.__contains__(key)

    def __getattr__(self, name):
        try:
            attr = getattr(self._dict, name)
            if callable(attr):
                return synchronized(self._lock)(attr)
            else:
                return attr
        except AttributeError:
            raise AttributeError(name)
        
class _ftpCallbackWrapper(object):
    def __init__(self, size, callback, progressCallback):
        self.size = max(size, 1)
        self.callback = callback
        self.progressCallback = progressCallback
        self.transferCount = 0

    def __call__(self, block):
        self.callback(block)
        self.progressCallback(min(int(100.0*self.transferCount/self.size), 100))
        self.transferCount+=len(block)

_monthDict = {"Jan":1, "Feb":2, "Mar":3, "Apr":4, "May":5, "Jun":6, "Jul":7, "Aug":8, "Sep":9, "Oct":10, "Nov":11, "Dec":12}

class FtpWorker(object):
    def __init__(self,  ftpAddr, statCache = None):
        self.ftpAddr = ftpAddr
        self.ftp = None
        self.now = datetime.now()
        self.statCache = statCache if statCache != None else SharedCache()

    def connect(self):
        if not self.ftp:
            self.ftp = ftplib.FTP()
        self.ftp.connect(self.ftpAddr)
        self.ftp.sock.settimeout(30)
        self.ftp.login()
        
    def retrieve(self, filename, local, update, progressCallback=None):
        local = os.path.normpath(local)
        isLocal = self.isLocal(local)
        if not update and isLocal:
            return
        retryCount = 0
        while retryCount<3:
            if not self.ftp:
                self.ftp = ftplib.FTP()
                self.connect()
            try:
                retryCount+=1
                size, date = self.statFtp(filename)
                #if update and update!="force":
                if isLocal:
                    sizeLocal, dateLocal = self.statLocal(local)
                    if sizeLocal!=size:
                        update = "force"
                    #print dateLocal, date, dateLocal < date
                    if dateLocal < date: # TODO fix date comparison
                        update = "force"
                if update=="force" or not isLocal:
                    s = StringIO()
##                    f = open(local + ".tmp", "wb")
                    if progressCallback:
                        self.ftp.retrbinary("RETR "+filename, _ftpCallbackWrapper(size, s.write, progressCallback), )
                    else:
                        self.ftp.retrbinary("RETR "+filename, s.write)
                    s.getvalue()
                    if s.len>size:
                        raise Exception("Wrong size of file "+filename)
                    f = open(local, "wb")
                    f.write(s.buf)
                    f.flush()
                    f.close()
##                    try:
##                        if os.path.exists(local):
##                            os.remove(local)
##                        os.rename(local + ".tmp", local)
##                    except Exception, ex:
##                        print ex, local
##                        raise
                    break
            except ftplib.error_perm, ex:
                if ex.args[0].startswith("550"):
                    self.connect()
                else:
                    raise
            except ftplib.error_temp, ex:
                if retryCount >= 3:
                    raise
                else:
                    time.sleep(3)
            except socket.error:
                if retryCount >= 3:
                    raise
                else:
                    self.connect()
            except FileNotFoundError:
                raise
    
    def isLocal(self, filename):
        try:
            open(filename)
            return True
        except:
            return False
        
    def statFtp(self, filename):
        if not self.ftp:
            self.connect()
        dir, file = os.path.split(filename)
        dir = dir + "/"
        with self.statCache._lock:
            if dir in self.statCache:
                if file:
                    try:
                        s = self.statCache[dir][file].split()
                    except KeyError:
                        raise FileNotFoundError(filename)
                else:
                    return
            else:
                lines = []
##                print "Ftp Stat:", dir, file
                self.ftp.dir(dir, lines.append)
                self.statCache[dir] = dict([(line.split()[-1].strip(), line.strip()) for line in lines if line.strip()])
                if file:
                    try:
                        s = self.statCache[dir][file].split()
                    except KeyError:
                        raise FileNotFoundError(filename)
                else:
                    return
##                print dir ,file, s
    ##            s = s.getvalue().split()
        size, date = int(s[-5]), s[-4:-1]
        if ":" in date[-1]:
            date = datetime(self.now.year, _monthDict.get(date[0], 1), int(date[1]))
            if date>self.now:
                date = datetime(date.year-1, date.month, date.day)
        else:
            date = datetime(int(date[-1]), _monthDict.get(date[0], 1), int(date[1]))
        return size, date

    def statLocal(self, filename):
        stat = os.stat(filename)
        return stat.st_size, datetime.fromtimestamp(stat.st_mtime)
    
    def listdir(self, ftp_dir):
        """ List the contents of a remote ftp directory (similar to os.listdir)
        """
        if not self.ftp:
            self.connect()
        lines = []
        self.ftp.dir(ftp_dir, lines.append)
        self.statCache[dir] = dict([(line.split()[-1].strip(), line.strip()) for line in lines if line.strip()])
        contents = [line.split()[-1] for line in lines]
        return [name for name in contents if name not in [".", ".."]] 
        

class FtpThreadWorker(threading.Thread, FtpWorker):
    def __init__(self, ftpAddr, queue, statCache=None, group=None, target=None, name=None, args=(), kwargs={}):
        threading.Thread.__init__(self, group, target, name, args, kwargs)
        FtpWorker.__init__(self, ftpAddr, statCache)
        self.queue = queue

    def run(self):
        while True:
            filename, local, update, retryCount, progressCallback = self.queue.get()
            try:
                self.retrieve(filename, local, update, progressCallback)
            except Exception, ex:
                sys.excepthook(*sys.exc_info())
            self.queue.task_done()

class FtpDownloader(object):
    def __init__(self, ftpAddr, localDir, ftpDir="", numOfThreads=5):
        self.ftpAddr = ftpAddr
        self.localDir = localDir
        self.ftpDir = ftpDir
        self.numOfThreads = numOfThreads
        self.queue = Queue(0)
        self.workers = []
        self.statCache = SharedCache()
        self.ftpWorker = FtpWorker(self.ftpAddr, statCache=self.statCache)

    def initWorkers(self):
        if self.workers:
            return
        for i in range(self.numOfThreads):
            try:
                t = FtpThreadWorker(self.ftpAddr, self.queue, self.statCache, name="FtpThread(%s):#%i"%(self.ftpAddr, i))
                t.setDaemon(True)
                t.start()
                self.workers.append(t)
            except Exception, ex:
                print ex
                print "Cannot create thread num: "+str(i)
                
    def massRetrieve(self, filenames, update=False, blocking=True, progressCallback=None):
        self.initWorkers()
        for filename in filenames:
            self.retrieve(filename, update, blocking=False)
##            localDir = os.path.split(self.localDir+filename)[0]
##            try:
##                os.makedirs(localDir)
##            except:
##                pass
##            self.queue.put((self.ftpDir+filename, self.localDir+filename, update, 0, None))
        if blocking:
            self.wait(progressCallback)
            
    def retrieve(self, filename, update=False, blocking=True, progressCallback=None):
        if type(filename) == str:
            filename = (filename, filename)
        localDir = os.path.split(os.path.join(self.localDir, filename[1]))[0]
        try:
            os.makedirs(localDir)
        except:
            pass
        if blocking:
            self.ftpWorker.retrieve(self.ftpDir+filename[0], os.path.join(self.localDir, filename[1]), update, progressCallback)
        else:
            self.queue.put((self.ftpDir+filename[0], os.path.join(self.localDir, filename[1]), update, 0, progressCallback))

    def wait(self, progressCallback=None):
        count = self.queue.qsize()
        while not self.queue.empty():
            if progressCallback:
                progressCallback(min(100.0, 100.0*(float(count)-self.queue.qsize())/count))
            time.sleep(0.1)
        self.queue.join()
        
    def listdir(self, ftp_dir):
        """ List the contents of the remote ftp dir (similar to os.listdir)
        """
        return self.ftpWorker.listdir(self.ftpDir + ftp_dir)
    
