from go import *

from obiGenomicsUpdate import Update as UpdateBase
from obiGenomicsUpdate import PKGManager as PKGManagerBase

import urllib2

from datetime import datetime
class Update(UpdateBase):
    def __init__(self, local_database_path=None, progressCallback=None):
        UpdateBase.__init__(self, local_database_path or getDataDir(), progressCallback)
    def CheckModified(self, addr, date=None):
        return date > self.GetLastModified(addr)
##        req = urllib2.Request(addr, headers=date and {"If-Modified-Since":date} or {})
##        try:
##            urllib2.urlopen(req)
##            return True
##        except urllib2.HTTPError, er:
####            if er.errno != 304:
####                print er
####                return False
##            return False
        
    def CheckModifiedOrg(self, org):
        return self.CheckModified("http://www.geneontology.org/gene-associations/gene_association." + org + ".gz", self.LastModifiedOrg(org))
    
    def LastModifiedOrg(self, org):
        return self.shelve.get((Update.UpdateAnnotation, (org,)), None)

    def GetLastModified(self, addr):
        stream = urllib2.urlopen(addr)
        return datetime.strptime(stream.headers.get("Last-Modified"), "%a, %d %b %Y %H:%M:%S %Z")
##        return stream.headers.get("Last-Modified")

    def IsUpdatable(self, func, args):
        if func == Update.UpdateOntology:
            return self.CheckModified("http://www.geneontology.org/ontology/gene_ontology.obo", self.shelve.get((Update.UpdateOntology, ()), None))
        elif func == Update.UpdateAnnotation:
            return self.CheckModifiedOrg(args[0])
            
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
        downloadAnnotationTo(org, os.path.join(self.local_database_path, "gene_association." + org), self.progressCallback)
        self._update(Update.UpdateAnnotation, (org,), self.GetLastModified("http://www.geneontology.org/gene-associations/gene_association." + org + ".gz"))

##    @synchronized(updateLock)
    def UpdateOntology(self):
        downloadGOTo(os.path.join(self.local_database_path, "gene_ontology.obo"), self.progressCallback)
        self._update(Update.UpdateOntology, (), self.GetLastModified("http://www.geneontology.org/ontology/gene_ontology.obo"))
    

class PKGManager(PKGManagerBase):
    def __init__(self, updater, *args, **kwargs):
        PKGManagerBase.__init__(self, updater, [], *args, **kwargs)
