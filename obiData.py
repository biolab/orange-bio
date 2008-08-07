import ftplib
import urllib, urllib2
import threading
import os
import time
import socket
from Queue import Queue
from StringIO import StringIO
from datetime import datetime

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
    def __init__(self,  ftpAddr):
        self.ftpAddr = ftpAddr
        self.ftp = None
        self.now = datetime.now()
        self.statCache = {}

    def connect(self):
        self.ftp.connect(self.ftpAddr)
        self.ftp.sock.settimeout(30)
        self.ftp.login()
        
    def retrieve(self, filename, local, update, progressCallback=None):
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
                    if progressCallback:
                        self.ftp.retrbinary("RETR "+filename, _ftpCallbackWrapper(size, s.write, progressCallback))
                    else:
                        self.ftp.retrbinary("RETR "+filename, s.write)
                    s.getvalue()
                    if s.len>size:
                        raise Exception("Wrong size of file "+filename)
                    f = open(local, "wb")
                    f.write(s.buf)
                    f.close()
                    break
            except ftplib.error_perm, ex:
                #print ex
                if ex.args[0].startswith("550"):
                    self.connect()
                else:
                    break
            except ftplib.error_temp, ex:
                #print ex
                break
            except socket.error:
                self.connect()
    
    def isLocal(self, filename):
        try:
            open(filename)
            return True
        except:
            return False
        
    def statFtp(self, filename):
        dir, file = os.path.split(filename)
##        print dir, file
        if (dir, file) in self.statCache:
            s = self.statCache[(dir, file)].split()
        else:
##            s = StringIO()
##            self.ftp.dir(filename, s.write)
            lines = []
            self.ftp.dir(dir, lines.append)
##            print "Lines : #", len(lines)
            self.statCache.update(dict([((dir, line.split()[-1].strip()), line.strip()) for line in lines if line.strip()]))
            s = self.statCache.get((dir, file), "-- 0 1 1 1 1").split()
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

class FtpThreadWorker(threading.Thread, FtpWorker):
    def __init__(self, ftpAddr, queue, group=None, target=None, name=None, args=(), kwargs={}):
        threading.Thread.__init__(self, group, target, name, args, kwargs)
        FtpWorker.__init__(self, ftpAddr)
        self.queue = queue

    def run(self):
        while True:
            filename, local, update, retryCount, progressCallback = self.queue.get()
            self.retrieve(filename, local, update, progressCallback)
            self.queue.task_done()

class DownloaderBase(object):
    def __init__(self, remoteAddr, localDir, remoteDir="", numOfThreads=5):
        self.remoteAddr = remoteAddr
        self.localDir = localDir
        self.remoteDir = remoteDir
        self.numOfThreads = numOfThreads

    def retrieve(self, *args):
        raise NotImplementedError

    def massRetrieve(self, *args):
        raise NotImplementedError

    def GetIfModified(self, filename, date):
        raise NotImplementedError
    
class FtpDownloader(object):
    def __init__(self, ftpAddr, localDir, ftpDir="", numOfThreads=5):
        self.ftpAddr = ftpAddr
        self.localDir = localDir
        self.ftpDir = ftpDir
        self.numOfThreads = numOfThreads
        self.queue = Queue(0)
        self.workers = []
        self.ftpWorker = FtpWorker(self.ftpAddr)

    def initWorkers(self):
        if self.workers:
            return
        for i in range(self.numOfThreads):
            try:
                t = FtpThreadWorker(self.ftpAddr, self.queue, name="FtpThread(%s):#%i"%(self.ftpAddr, i))
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
            while not self.queue.empty():
                if progressCallback:
                    progressCallback(min(100.0, 100.0*(float(len(filenames))-self.queue.qsize())/len(filenames)))
                time.sleep(0.1)

    def retrieve(self, filename, update=False, blocking=True, progressCallback=None):
        localDir = os.path.split(self.localDir+filename)[0]
        try:
            os.makedirs(localDir)
        except:
            pass
        if blocking:
            self.ftpWorker.retrieve(self.ftpDir+filename, self.localDir+filename, update, progressCallback)
        else:
            self.queue.put((self.ftpDir+filename, self.localDir+filename, update, 0, progressCallback))

##class HTTPDownloader(DownloaderBase):
##    def __init__(self, *args, **kwargs):
##        DownloaderBase.__init__(self, *args, **kwargs)
##
##    def retrieve(self, filename):
##        
