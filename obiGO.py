from go import *

from obiGenomicsUpdate import Update as UpdateBase
from obiGenomicsUpdate import PKGManager as PKGManagerBase
##from obiGenomicsUpdate import synchronized
##
##from threading import Lock
import urllib2
##updateLock = Lock()

##class Update(UpdateBase):
##    def __init__(self, local_database_path=None, progressCallback=None):
##        UpdateBase.__init__(self, local_database_path or getDataDir(), progressCallback)
##    @synchronized(updateLock)
##    def GetUpdatable(self):
##        orgs = listDownloadedOrganisms()
##        return [(Update.UpdateAnnotation, "Update organism annotation", orgs), (Update.UpdateOntology, "Update ontology", [])]
##
##    @synchronized(updateLock)
##    def GetDownloadable(self):
##        orgs = set(listOrganisms())-set(listDownloadedOrganisms())
##        
##        return [(Update.UpdateAnnotation, "Download organism annotation", list(orgs)), (Update.UpdateOntology, "Downlaod ontology", [])]
##
##    @synchronized(updateLock)
##    def UpdateAnnotation(self, org):
##        downloadAnnotationTo(org, self.local_database_path+"//gene_association."+org, self.progressCallback)
##        self._update(Update.UpdateAnnotation, (org,))
##
##    @synchronized(updateLock)      
##    def UpdateOntology(self):
##        downloadGOTo(self.local_database_path+"//gene_ontology.obo", self.progressCallback)
##        self._update(Update.UpdateOntology, ())

class Update(UpdateBase):
    def __init__(self, local_database_path=None, progressCallback=None):
        UpdateBase.__init__(self, local_database_path or getDataDir(), progressCallback)
    def CheckModified(self, addr, date=None):
        req = urllib2.Request(addr, headers=date and {"If-Modified-Since":date} or {})
        try:
            urllib2.urlopen(req)
            return False
        except urllib2.HTTPError, er:
            if er.errno != 304:
                print er
                return False
            return True
        
    def CheckModifiedOrg(self, org):
        return self.CheckModified("http://www.geneontology.org/gene-associations/gene_association." + org + ".gz", self.LastModifiedOrg(org))
    
    def LastModifiedOrg(self, org):
        return self.shelve.get((Update.UpdateAnnotation, (org,)), None)

    def GetLastModified(self, addr):
        stream = urllib2.urlopen(addr)
        return stream.headers.get("Last-Modified")

    def IsUpdatable(self, func, args):
        if func == Update.UpdateOntology:
            return self.CheckModified("http://www.geneontology.org/ontology/gene_ontology.obo", self.shelve.get((Update.UpdateOntology, ()), None))
        elif func == Update.UpdateAnnotation:
            return self.CheckModifiedOrg(args[0])
            
##    @synchronized(updateLock)
##    def GetUpdatable(self):
##        orgs = [org for org in listDownloadedOrganisms() if self.CheckModifiedOrg(org)]
##        ret = []
##        if (Update.UpdateOntology, ()) in self.shelve and self.CheckModified("http://www.geneontology.org/ontology/gene_ontology.obo", self.shelve.get((Update.UpdateOntology, ()), None)):
##            ret.append((Update.UpdateOntology, ()))
##        if orgs:
##            ret.extend([(Update.UpdateAnnotation, (org,)) for org in orgs])
##        return ret

##    @synchronized(updateLock)
    def GetDownloadable(self):
        orgs = set(listOrganisms())-set(listDownloadedOrganisms())
        ret = []
        if (Update.UpdateOntology, ()) not in self.shelve:
            ret.append((Update.UpdateOntology, ()))
        if orgs:
            ret.extend([(Update.UpdateAnnotation, (org,)) for org in orgs])
        return ret

##    @synchronized(updateLock)
    def UpdateAnnotation(self, org):
        downloadAnnotationTo(org, self.local_database_path+"//gene_association."+org, self.progressCallback)
        self._update(Update.UpdateAnnotation, (org,), self.GetLastModified("http://www.geneontology.org/gene-associations/gene_association." + org + ".gz"))

##    @synchronized(updateLock)
    def UpdateOntology(self):
        downloadGOTo(self.local_database_path+"//gene_ontology.obo", self.progressCallback)
        self._update(Update.UpdateOntology, (), self.GetLastModified("http://www.geneontology.org/ontology/gene_ontology.obo"))
    

class PKGManager(PKGManagerBase):
    def __init__(self, updater, compression="gz"):
        PKGManagerBase.__init__(self, updater, [], compression)
