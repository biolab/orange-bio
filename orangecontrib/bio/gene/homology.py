import sys, os
import urllib2

from collections import defaultdict

from Orange.orng import orngServerFiles

class _homolog(object):
    __slots__ = ["group_id", "taxonomy_id", "gene_id", "gene_symbol"]
    def __init__(self, homolog_line):
        for attr, val in zip(self.__slots__, homolog_line.split("\t")):
            setattr(self, attr, val)
    
    
class _Homologs(object):
    """ A base class for homolog mappers
    """
    def all_genes(self):
        """ Return all genes in this class instance
        """
        raise NotImplemented
    
    def homologs(self, gene, taxid):
        """ Return all (homolotaxid, homolog) pairs of gene from organism with taxid
        """
        raise NotImplemented
    
    def homolog(self, gene, taxid, homolotaxid):
        """ Return homolog of gene from organism with *taxid* in organism with "homolotaxid*
        """
        homologs = dict(self.homologs(gene, taxid))
        return homologs.get(homolotaxid, None)
    
class HomoloGene(_Homologs):
    DEFAULT_DATABASE_PATH = orngServerFiles.localpath("HomoloGene")
    VERSION = 1
    DOMAIN = "HomoloGene"
    FILENAME = "homologene.data"
    def __init__(self, local_database_path=None):
        self.local_database_path = local_database_path if local_database_path else self.DEFAULT_DATABASE_PATH
        self.load()

    @classmethod
    def download_from_NCBI(cls, file=None):
        data = urllib2.urlopen("ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/current/homologene.data").read()
        if file is None:
            try:
                os.mkdir(orngServerFiles.localpath("HomoloGene"))
            except OSError:
                pass
            file = open(orngServerFiles.localpath("HomoloGene", "homologene.data"), "wb")
        elif type(file) in [str, unicode]:
            file = open(file, "wb")
        file.write(data)
        file.flush()
        
    @classmethod    
    def get_instance(cls):
        if not hasattr(cls, "_shared_dict"):
            h = cls()
            cls._shared_dict = h.__dict__
        h = cls.__new__(cls)
        h.__dict__ = cls._shared_dict
        return h
    
    def load(self):
        path = orngServerFiles.localpath_download(self.DOMAIN, self.FILENAME)
        lines = open(path, "rb").read().splitlines()[:-1]
        self._homologs = {} 
        self._homologs = dict([((h.taxonomy_id, h.gene_symbol), h) for h in [_homolog(line) for line in lines]])
        self._homologs_by_group = reduce(lambda dict, h: dict[h.group_id].append(h) or dict, self._homologs.values(), defaultdict(list))
#        for line in lines:
#            h = _homolog(line)
#            self._homologs[h.taxonomy_id, h.gene_symbol] = h
#            self._homologs_by_group[h.group_id].append(h)
        
    def all_genes(self, taxid=None):
        return [homolog.gene_symbol for (tid, id), homolog in self._homologs.iteritems() if tid == taxid]
    
    def homologs(self, gene, taxid):
        group = self._homologs.get((taxid, gene), _homolog("")).group_id
        homologs = self._homologs_by_group[group]
        return [(h.taxonomy_id, h.gene_symbol) for h in homologs]
        
    def homolog(self, gene, taxid, homolotaxid):
        homologs = dict(self.homologs(gene, taxid))
        return homologs.get(homolotaxid, None)
        
def _parseOrthoXML(file):
    """ Return (cluster_id, taxid, gene_id) tuples from orthoXML file 
    """
    from xml.dom.minidom import parse
    doc = parse(file)
    species = doc.getElementsByTagName("species")
    geneIds = {}
    geneIdToTaxid = {}
    for sp in species:
        taxid = sp.attributes.get("NCBITaxId").value
        genes = sp.getElementsByTagName("gene")
        geneIds.update(dict([(gene.attributes.get("id").value, (gene.attributes.get("geneId").value,
                  gene.attributes.get("protId").value)) for gene in genes]))
        geneIdToTaxid.update(dict.fromkeys([gene.attributes.get("geneId").value for gene in genes], taxid))
        
    orthologs = []
    clusters = doc.getElementsByTagName("cluster")
    for cl in clusters:
        clId = cl.attributes.get("id").value
        geneRefs = cl.getElementsByTagName("geneRef")
        ids = [ref.attributes.get("id").value for ref in geneRefs]
        orthologs.extend([(clId, geneIdToTaxid[geneIds[id][0]], geneIds[id][0]) for id in ids])
    return orthologs
        
class InParanoid(object):
    """ InParanoid: Eukaryotic Ortholog Groups
    """
    VERSION = 1
    def __init__(self):
        import sqlite3
        self.con = sqlite3.connect(orngServerFiles.localpath_download("HomoloGene", "InParanoid.sqlite"))
        
    def all_genes(self, taxid):
        """ Return all genes in the database for the given taxid
        """
        return [t[0] for t in self.con.execute("select distinct geneid from homologs where homologs.taxid=?", (taxid,)).fetchall()]
    
    def all_taxids(self):
        """ Return all taxids in the database
        """
        return [t[0] for t in self.con.execute("select distinct taxid from homologs").fetchall()]
    
    def _groups(self, gene, taxid):
        """ Return all group identifiers for gene, taxid pair
        """
        return self.con.execute("select distinct groupid from homologs where homologs.taxid=? and homologs.geneid=?", (taxid, gene)).fetchall()
    
    def orthologs(self, gene, taxid, ortholog_taxid=None):
        """ Return all orthologs of genename from organism with taxid. 
        If ortholog_taxid is given limit to orthologs from that organism only
        """
        groups = self._groups(gene, taxid)
        res = []
        for group in groups:
            if ortholog_taxid:
                res.extend(self.con.execute("select distinct taxid, geneid from homologs where homologs.groupid=? and homologs.taxid=?", (group[0], ortholog_taxid)).fetchall())
            else:
                res.extend(self.con.execute("select distinct taxid, geneid from homologs where homologs.groupid=?", group).fetchall())
        res = sorted(set(res))
        if ortholog_taxid:
            res = [r[1] for r in res]
        return res
        
def all_genes(taxid):
    """ Return a set of all genes for organism taxid.
    """
    return HomoloGene.get_instance().all_genes(taxid)

def homologs(genename, taxid):
    """ Return a list of homologs (taxid, genename) for a homolog group of gene (organism taxid).
    """ 
    return HomoloGene.get_instance().homologs(genename, taxid)

def homolog(genename, taxid, homolotaxid):
    """ Return a homolog of genename (for taxid) in organism holomotaxid. If the homolog does not exist, return None.
    """
    return HomoloGene.get_instance().homolog(genename, taxid, homolotaxid)

def all_genes_inParanoid(taxid):
    """ Return a set of all genes for organism with taxid in the InParanoid database.
    """
    return InParanoid().all_genes(taxid)

def orthologs(genename, taxid, ortholog_taxid=None):
    """ Return all InParanoid orthologs of genename from organism with taxid. 
    If ortholog_taxid is given limit to orthologs from that organism only.
    """
    return InParanoid().orthologs(genename, taxid, ortholog_taxid)
