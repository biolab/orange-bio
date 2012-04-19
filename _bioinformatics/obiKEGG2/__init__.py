"""\
==============================================
KEGG - Kyoto Encyclopedia of Genes and Genomes
==============================================

This is a python module for access to `KEGG`_ using its web services. To use this module you need to have
`SUDS`_ python library installed (other backends are planed). 

.. _`KEGG`: http://www.genome.jp/kegg/

.. _`SUDS`: http://pypi.python.org/pypi/suds/

"""
from __future__ import absolute_import


import os, sys
from collections import defaultdict

from datetime import datetime

from Orange.utils import lru_cache, serverfiles

from . import databases
from . import entry

from .brite import BriteEntry, Brite

from . import api
from . import conf
from . import pathway

KEGGGenome = databases.Genome
KEGGGenes = databases.Genes
KEGGEnzymes = databases.Enzymes
KEGGReaction = databases.Reactions
KEGGPathways = databases.Pathways

KEGGBrite = Brite
KEGGBriteEntry = BriteEntry

KEGGPathway = pathway.Pathway

DEFAULT_CACHE_DIR = conf.params["cache.path"]


from .. import obiProb
from Orange.utils import deprecated_keywords, deprecated_attribute

class Organism(object):
    def __init__(self, org, genematcher=None):
        self.org_code = self.organism_name_search(org)
        self.genematcher = genematcher
        self.api = api.CachedKeggApi()
        
    @property
    def org(self):
        return self.org_code
    
    @property
    def genes(self):
        if not hasattr(self, "_genes"):
            genes = KEGGGenes(self.org_code)
            self._genes = genes
        return self._genes
    
    def gene_aliases(self):
        return self.genes().gene_aliases()
    
    def pathways(self, with_ids=None):
        if with_ids is not None:
            return self.api.get_pathways_by_genes(with_ids)
        else:
            return [p.entry_id for p in self.api.list_pathways(self.org_code)]
    
    def list_pathways(self):
        return self.pathways()
    
    def get_linked_pathways(self, pathway_id):
        self.api.get_linked_pathways(pathway_id)
        
    def enzymes(self, genes=None):
        raise NotImplementedError()
    
    def get_enriched_pathways(self, genes, reference=None, prob=obiProb.Binomial(), callback=None):
        """ Return a dictionary with enriched pathways ids as keys
        and (list_of_genes, p_value, num_of_reference_genes) tuples 
        as items.
        
        """
        allPathways = defaultdict(lambda :[[], 1.0, []])
        from Orange.orng import orngMisc
        milestones = orngMisc.progressBarMilestones(len(genes), 100)
        pathways_db = KEGGPathways()
        
        pathways_for_gene = []
        for i, gene in enumerate(genes):
            pathways_for_gene.append(self.pathways([gene]))
            if callback and i in milestones:
                callback(i*50.0/len(genes))
                
        # precache for speed 
        pathways_db.pre_cache([pid for pfg in pathways_for_gene for pid in pfg]) 
        for i, (gene, pathways) in enumerate(zip(genes, pathways_for_gene)):
            for pathway in pathways:
                if pathways_db.get_entry(pathway).gene: 
                    allPathways[pathway][0].append(gene)
            if callback and i in milestones:
                callback(50.0 + i*50.0/len(genes))
        reference = set(reference if reference is not None else self.genes.keys())
        
        pItems = allPathways.items()
        
        for i, (p_id, entry) in enumerate(pItems):
            pathway = pathways_db.get_entry(p_id)
            entry[2].extend(reference.intersection(pathway.gene or []))
            entry[1] = prob.p_value(len(entry[0]), len(reference), len(entry[2]), len(genes))
        return dict([(pid, (genes, p, len(ref))) for pid, (genes, p, ref) in allPathways.items()])
        
    def get_genes_by_enzyme(self, enzyme):
        enzyme = Enzymes().get_entry(enzyme)
        return enzyme.genes.get(self.org_code, []) if enzyme.genes else []
    
    def get_genes_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).genes()
    
    def get_enzymes_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).enzymes()
    
    def get_compounds_by_pathway(self, pathway_id):
        return KEGGPathway(pathway_id).compounds()
    
    def get_pathways_by_genes(self, gene_ids):
        return self.api.get_pathways_by_genes(gene_ids)
        gene_ids = set(gene_ids)
        pathways = [self.genes[id].pathway for id in gene_ids if self.genes[id].pathway]
        pathways = reduce(set.union, pathways, set())
        return [id for id in pathways if gene_ids.issubset(KEGGPathway(id).genes())] 
    
    def get_pathways_by_enzymes(self, enzyme_ids):
        enzyme_ids = set(enzyme_ids)
        pathways = [KEGGEnzymes()[id].pathway for id in enzyme_ids]
        pathwats = reduce(set.union, pathways, set())
        return [id for id in pathways if enzyme_ids.issubset(KEGGPathway(id).enzymes())]
    
    def get_pathways_by_compounds(self, compound_ids):
        compound_ids = set(compound_ids)
        pathways = [KEGGCompounds()[id].pathway for id in compound_ids]
        pathwats = reduce(set.union, pathways, set())
        return [id for id in pathways if compound_ids.issubset(KEGGPathway(id).compounds())]
    
    def get_enzymes_by_compound(self, compound_id):
        return KEGGCompound()[compound_id].enzyme
    
    def get_enzymes_by_gene(self, gene_id):
        return self.genes[gene_id].enzymes
    
    def get_compounds_by_enzyme(self, enzyme_id):
        return self._enzymes_to_compounds.get(enzyme_id)
    
    @deprecated_keywords({"caseSensitive": "case_sensitive"})
    def get_unique_gene_ids(self, genes, case_sensitive=True):
        """Return a tuple with three elements. The first is a dictionary mapping from unique gene
        ids to gene names in genes, the second is a list of conflicting gene names and the third is a list
        of unknown genes.
        """
        unique, conflicting, unknown = {}, [], []
        for gene in genes:
            names = self.genematcher.match(gene)
            if len(names) == 1:
                unique[names[0]] = gene
            elif len(names) == 0:
                unknown.append(gene)
            else:
                conflicting.append(gene)
        return unique, conflicting, unknown
    
    def get_genes(self):
        return self.genes
    
    @classmethod
    def organism_name_search(cls, name):
        genome = KEGGGenome()
        if name not in genome:
            ids = genome.search(name)
            if not ids:
                from .. import obiTaxonomy
                ids = obiTaxonomy.search(name)
                ids = [id for id in ids if genome.search(id)]
            name = ids.pop(0) if ids else name
            
        try:
            return genome[name].entry_key
        except KeyError:
            raise ValueError("Organism with name='%s' not found in KEGG." % name)
        
    @classmethod
    def organism_version(cls, name):
        name = cls.organism_name_search(name)
        genome = KEGGGenome()
        info = genome.api.binfo(name)
        return info.release
#        orngServerFiles.localpath_download("KEGG", "kegg_genes_%s.tar.gz" % name)
#        return orngServerFiles.info("KEGG", "kegg_genes_%s.tar.gz" % name)["datetime"]
    
    def _set_genematcher(self, genematcher):
        setattr(self, "_genematcher", genematcher)
        
    def _get_genematcher(self):
        if getattr(self, "_genematcher", None) == None:
            from .. import obiGene
            if self.org_code == "ddi":
                self._genematcher = obiGene.matcher([obiGene.GMKEGG(self.org_code), obiGene.GMDicty(),
                                                     [obiGene.GMKEGG(self.org_code), obiGene.GMDicty()]])
            else:
                self._genematcher = obiGene.matcher([obiGene.GMKEGG(self.org_code)])
            self._genematcher.set_targets(self.genes.keys())
        return self._genematcher
    
    genematcher = property(_get_genematcher, _set_genematcher)
    
KEGGOrganism = Organism
    
def organism_name_search(name):
    return KEGGOrganism.organism_name_search(name)

def pathways(org):
    return KEGGPathway.list(org)

def organisms():
    return KEGGOrganism.organisms()

def from_taxid(taxid):
    genome = KEGGGenome()
    res = genome.search(taxid)
    print taxid, res
    for r in res:
        e = genome[r]
        
        if e.taxid in [taxid,  genome.TAXID_MAP.get(taxid, taxid)]:
            return e.org_code()
        
    return None

def to_taxid(name):
    genome = KEGGGenome()
    if name in genome:
        return genome[name].taxid
    
    keys = genome.search(name)
    if keys:
        return genome[keys[0]].taxid
    else:
        return None

def create_gene_sets():
    pass

from .. import obiGene
from Orange.utils import ConsoleProgressBar

class MatcherAliasesKEGG(obiGene.MatcherAliasesPickled):
    DOMAIN = "KEGG"
    VERSION = "v3.0"
    def create_aliases(self):
        import cPickle
        files = set(serverfiles.ServerFiles().listfiles(self.DOMAIN))
        ids_filename = "kegg_gene_id_aliases_" + self.organism + ".pickle"
        if ids_filename in files:
            filename = serverfiles.localpath_download(self.DOMAIN, ids_filename)
            
            aliases = cPickle.load(open(filename, "rb"))
        else:
            pb = ConsoleProgressBar("Retriving KEGG ids:")
            kegg_org = KEGGOrganism(self.organism)
            genes = kegg_org.genes
            genes.pre_cache(progress_callback=pb.set_state)
            aliases = []
            for key, entry in genes.iteritems():
                aliases.append(set([key]) | set(entry.alt_names))
            filename = serverfiles.localpath_download(self.DOMAIN, ids_filename)
            cPickle.dump(aliases, open(filename, "wb"))
            
        return aliases
    
    def filename(self):
        return "kegg3_" + self.organism
    
    def aliases_path(self):
        ids_filename = "kegg_gene_id_aliases_" + self.organism + ".pickle"
        return serverfiles.localpath(self.DOMAIN, ids_filename)
    
    def create_aliases_version(self):
        files = set(serverfiles.listfiles(self.DOMAIN))
        ids_filename = "kegg_gene_id_aliases_" + self.organism + ".pickle"
        if ids_filename in files:
            version = serverfiles.info(self.DOMAIN, ids_filename)["datetime"]
        else:
            kegg_org = KEGGOrganism(self.organism)
            genes = kegg_org.genes
            version = genes.info.release
        return version
        
    def __init__(self, organism, **kwargs):
        self.organism = organism
        sf = serverfiles.ServerFiles()
        files = set(sf.listfiles(self.DOMAIN))
        ids_filename = "kegg_gene_id_aliases_" + self.organism + ".pickle"
        if ids_filename in files:
            serverfiles.update(self.DOMAIN, ids_filename)
            
        obiGene.MatcherAliasesPickled.__init__(self, **kwargs)

def main():
    KEGGGenome()
    import doctest
    extraglobs = {"api": KeggApi()}
    doctest.testmod(optionflags=doctest.ELLIPSIS, extraglobs=extraglobs)

if __name__ == "__main__":
    sys.exit(main())
