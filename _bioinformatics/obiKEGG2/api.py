"""
KEGG api interface.

"""
from __future__ import absolute_import

from contextlib import closing

from .service import web_service
from .types import *

class KeggApi(object):
    """ KEGG API """
    
    def __init__(self):
        self.service = web_service()
        
    ##################
    # Meta information
    ##################
    
    def list_databases(self):
        """ Returns a list of available databases.
        
        >>> api.list_databases()
        [Definition(entry_id='nt',...
         
        """
        return map(Definition.from_items, self.service.list_databases())
    
    def list_organisms(self):
        """ Return a list of all available organisms
        
        >>> api.list_organisms()
        [Definition(entry_id='hsa',...
        
        """
        return map(Definition.from_items, self.service.list_organisms())
    
    def list_pathways(self, organism):
        """ Return a list of all available pathways for `organism`
        
        >>> api.list_pathways("hsa")
        [Definition(entry_id=',...
        
        """
        return map(Definition.from_items, self.service.list_pathways(organism))
        
    #######
    # DBGET
    #######
     
    def binfo(self, db):
        """ Return info for database `db`
        
        >>> print api.dbinfo("gb")
        genbank          GenBank nucleic acid sequence database
        gb               Release 186.0, Oct 11
                         National Center for Biotechnology Information
                         144,458,648 entries, 132,067,413,372 bases
                         Last update: 11/10/24
                         <dbget> <fasta> <blast>
                         
        """
        result = self.service.binfo(db)
        if result is not None:
            return BInfo.from_text(str(result))
        else:
            return result
    
    def bfind(self, db, keywords):
        """ Search database 'db' for keywords
        """
        result = self.service.bfind(" ".join([db, keywords]))
        if result is not None:
            return str(result)
        else:
            return result
    
    def bget(self, ids):
        """
        """
        if not isinstance(ids, basestring):
            # Sequence of ids
            ids = " ".join(ids)
        result = self.service.bget(ids)
        if result is not None:
            return str(result)
        else:
            return result
    
    def btit(self, ids):
        """
        """
        if not isinstance(ids, basestring):
            ids = " ".join(ids)
            
        result = self.service.btit(ids)
        if result is not None:
            return str(result)
        else:
            return result
    
    def bconv(self, ids):
        if not isinstance(ids, basestring):
            ids = " ".join(ids)
            
        result = self.service.bconv(ids)
        if result is not None:
            return str(result)
        else:
            return result
    
    ########
    # LinkDB
    ########
    
    def get_linkdb_by_entry(self, entry_id, db, offset, limit):
        links = self.service.get_linkdb_by_entry(entry_id, db, offset, limit)
        return [LinkDBRelation(**d) for d in \
                map(dict, links)]
        
    def get_linkdb_between_databases(self, from_db, to_db, offset, limit):
        links = self.service.get_linkdb_between_databases(from_db, to_db, offset, limit)
        return [LinkDBRelation(**d) for d in \
                map(dict, links)]
        
    def get_genes_by_enzyme(self, enzyme_id, org):
        return self.service.get_genes_by_enzyme(enzyme_id, org)
    
    def get_enzymes_by_gene(self, genes_id):
        return self.service.get_enzymes_by_gene(genes_id)
    
    def get_enzymes_by_compound(self, compound_id):
        return self.service.get_enzymes_by_compound(compound_id)
    
    def get_enzymes_by_glycan(self, glycan_id):
        return self.service.get_enzymes_by_glycan(glycan_id)
    
    def get_enzymes_by_reaction(self, reaction_id):
        return self.service.get_enzymes_by_reaction(reaction_id)
    
    def get_compounds_by_enzyme(self, enzyme_id):
        return self.service.get_compounds_by_enzyme(enzyme_id)
    
    def get_compounds_by_reaction(self, reaction_id):
        return self.service.get_compounds_by_reaction(reaction_id)
    
    def get_glycans_by_enzyme(self, enzyme_id):
        return self.service.get_glycans_by_enzyme(enzyme_id)
    
    def get_glycans_by_reaction(self, reaction_id):
        return self.service.get_glycans_by_reaction(reaction_id)
    
    def get_reactions_by_enzyme(self, enzyme_id):
        return self.service.get_reactions_by_enzyme(enzyme_id)
    
    def get_reactions_by_compound(self, compound_id):
        return self.service.get_reactions_by_compound(compound_id)
    
    def get_reactions_by_glycan(self, glycan_id):
        return self.service.get_reactions_by_glycan(glycan_id)
    
    ######
    # SSDB
    ######
    
    def get_best_best_neighbors_by_gene(self, genes_id, offset, limit):
        ssr = self.service.get_best_best_neighbors_by_gene(genes_id, offset, limit)
        return [SSDBRelation(**d) for d in \
                map(dict, ssr)]
    
    def get_best_neighbors_by_gene(self, genes_id, offset, limit):
        ssr = self.service.get_best_neighbors_by_gene(genes_id, offset, limit)
        return [SSDBRelation(**d) for d in \
                map(dict, ssr)]
    
    def get_reverse_best_neighbors_by_gene(self, genes_id, offset, limit):
        ssr = self.service.get_reverse_best_neighbors_by_gene(genes_id, offset, limit)
        return [SSDBRelation(**d) for d in \
                map(dict, ssr)]
    
    def get_paralogs_by_gene(self, genes_id, offset, limit):
        ssr =  self.service.get_paralogs_by_gene(genes_id, offset, limit)
        return [SSDBRelation(**d) for d in \
                map(dict, ssr)]
    
    #######
    # Motif
    #######
    
    def get_motifs_by_gene(self, genes_id, db):
        motif = self.service.get_motifs_by_gene(genes_id, db)
        return [MotifResult(**d) for d in \
                map(dict, motif)]
    
    def get_genes_by_motifs(self, motif_id_list, offset, limit):
        genes = self.service.get_genes_by_motifs(motif_id_list, offset, limit)
        return [Definition(**d) for d in \
                map(dict, genes)]
    
    ####
    # KO
    ####
    
    def get_ko_by_gene(self, genes_id):
        return self.service.get_ko_by_gene(genes_id)
    
    def get_ko_by_ko_class(self, ko_class_id):
        return self.service.get_ko_by_ko_class(ko_class_id)
    
    def get_genes_by_ko_class(self, ko_class_id, org, offset, limit):
        return self.service.get_genes_by_ko_class(ko_class_id, org, offset, limit)
    
    def get_genes_by_ko(self, ko_id, org):
        return self.service.get_genes_by_ko(ko_id, org)
    
    #########
    # Pathway
    #########
    
    def mark_pathway_by_objects(self, pathway_id, object_id_list):
        return self.service.mark_pathway_by_objects(pathway_id, object_id_list)
    
    def color_pathway_by_objects(self, pathway_id, object_id_list, fg_color_list, bg_color_list):
        return self.service.color_pathway_by_objects(pathway_id, object_id_list, fg_color_list, bg_color_list)
    
    def color_pathway_by_elements(self, pathway_id, element_id_list, fg_color_list, bg_color_list):
        return self.service.color_pathway_by_elements(pathway_id, element_id_list, fg_color_list, bg_color_list)
    
    def get_html_of_marked_pathway_by_objects(self, pathway_id, object_id_list):
        return self.service.get_html_of_marked_pathway_by_objects(pathway_id, object_id_list)
    
    def get_html_of_colored_pathway_by_objects(self, pathway_id, object_id_list, fg_color_list, bg_color_list):
        return self.service.get_html_of_colored_pathway_by_objects(pathway_id, object_id_list, fg_color_list, bg_color_list)
    
    def get_html_of_colored_pathway_by_elements(self, pathway_id, element_id_list, fg_color_list, bg_color_list):
        return self.service.get_html_of_colored_pathway_by_elements(pathway_id, element_id_list, fg_color_list, bg_color_list)
    
    def get_references_by_pathway(self, pathway_id):
        return self.service.get_references_by_pathway(pathway_id)
    
    def get_element_relations_by_pathway(self, pathway_id):
        return self.service.get_element_relations_by_pathway(pathway_id)
    
    
    
    def get_genes_by_organism(self, organism, offset=None, limit=None):
        if offset is None and limit is None:
            offset = 0
            limit = self.get_number_of_genes_by_organism(organism)
            
        return self.service.get_genes_by_organism(organism, offset, limit)
    
    def get_number_of_genes_by_organism(self, organism):
        return self.service.get_number_of_genes_by_organism(organism)
    
    ####################
    # Objects by pathway
    ####################
    
    def get_elements_by_pathway(self, pathway_id):
        return self.service.get_elements_by_pathway(pathway_id)
    
    def get_genes_by_pathway(self, pathway_id):
        return self.service.get_genes_by_pathway(pathway_id)
    
    def get_enzymes_by_pathway(self, pathway_id):
        return self.service.get_enzymes_by_pathway(pathway_id)
    
    def get_compounds_by_pathway(self, pathway_id):
        return self.service.get_compounds_by_pathway(pathway_id)
    
    def get_drugs_by_pathway(self, pathway_id):
        return self.service.get_drugs_by_pathway(pathway_id)
    
    def get_glycans_by_pathway(self, pathway_id):
        return self.service.get_glycans_by_pathway(pathway_id)
    
    def get_reactions_by_pathway(self, pathway_id):
        return self.get_reactions_by_pathway(pathway_id)
    
    def get_kos_by_pathway(self, pathway_id):
        return self.service.get_kos_by_pathway(pathway_id)
    
    #####################
    # Pathways by objects
    #####################
    
    def get_pathways_by_genes(self, gene_list):
        return map(str, self.service.get_pathways_by_genes(gene_list))
    
    def get_pathways_by_enzymes(self, enzyme_list):
        return map(str, self.service.get_pathways_by_enzymes(enzyme_list))
    
    def get_pathways_by_compounds(self, compound_list):
        return map(str, self.service.get_pathways_by_compounds(compound_list))
    
    def get_pathways_by_drugs(self, drug_list):
        return map(str, self.service.get_pathways_by_drugs(drug_list))
    
    def get_pathways_by_glycans(self, glycan_list):
        return map(str, self.service.get_pathways_by_glycans(glycan_list))
    
    def get_pathways_by_reactions(self, reaction_list):
        return map(str, self.service.get_pathways_by_reactions(reaction_list))
    
    def get_pathways_by_kos(self, ko_list):
        return map(str, self.service.get_pathways_by_kos(ko_list))
    
    ##########################
    # Relations among pathways
    ##########################
    
    def get_linked_pathways(self, pathway_id):
        if not pathway_id.startswith("path:"):
            pathway_id = "path:" + pathway_id
        return map(str, self.service.get_linked_pathways(pathway_id))
    
    
"""
KEGG api with caching
"""

import os

from . import caching
from .caching import cached_method, cache_entry, touch_dir

try:
    from functools import lru_cache
except ImportError:
    # TODO: move a copy of lru_cache in .caching if distributing this as a
    # standalone package
    from Orange.utils import lru_cache

    
class CachedKeggApi(KeggApi):
    def __init__(self, store=None):
        KeggApi.__init__(self)
        if store is None:
            self.store = {}
    
    # Needed API for cached decorator.
    def cache_store(self):
        from . import conf
        path = conf.params["cache.path"]
        touch_dir(path)
        return caching.Sqlite3Store(os.path.join(path,
                                                 "kegg_api_cache.sqlite3"))
    
    def last_modified(self, args, kwargs=None):
        return getattr(self, "default_release", "")
    
    def set_default_release(self, release):
        self.default_release = release
        
    
    ##################
    # Meta information
    ##################
    
    @lru_cache() # not persistently cached
    def list_databases(self):
        return KeggApi.list_databases(self)
    
    @cached_method
    def list_organisms(self):
        return KeggApi.list_organisms(self)
    
    @cached_method
    def list_pathways(self, organism):
        return KeggApi.list_pathways(self, organism)
    
    #######
    # DBGET
    #######
    
    @lru_cache() # not persistently cached
    def binfo(self, db):
        return KeggApi.binfo(self, db)
    
    @cached_method
    def bfind(self, db, keywords):
        return KeggApi.bfind(self, db, keywords)
    
    @cached_method
    def bget(self, ids):
        rval = KeggApi.bget(self, ids)
        return rval
    
    @cached_method
    def bget(self, ids):
        if not isinstance(ids, basestring):
            return self._batch_bget(ids)
        else:
            return KeggApi.bget(self, ids)
        
    def _batch_bget(self, ids):
        if len(ids) > 100:
            raise ValueError("Can batch at most 100 ids at a time.")
        
        bget = self.bget
        uncached = []
        with closing(bget.cache_store()) as store:
            # Which ids are already cached
            # TODO: Invalidate entries by release string.
            for id in ids:
                key = bget.key_from_args((id,))
                if key not in store:
                    uncached.append(id)
                
        if uncached:
            # in case there are duplicate ids
            uncached = sorted(set(uncached))
            rval = KeggApi.bget(self, uncached)
            if rval is not None:
                entrys = rval.split("///\n")
            else:
                entrys = []
                
            if entrys and not entrys[-1].strip():
                # Delete the last newline if present
                del entrys[-1]
            
            if len(entrys) == len(uncached):
                with closing(bget.cache_store()) as store:
                    for id, entry in zip(uncached, entrys):
                        key = bget.key_from_args((id,))
                        if entry is not None:
                            entry = entry + "///\n"
                        store[key] = cache_entry(entry, mtime=datetime.now())
                        
            else:
                # Try to bisect the uncached list
                if len(uncached) > 1 and len(uncached) - len(entrys) < 4:
                    split = len(uncached) / 2
                    self._batch_bget(uncached[:split])
                    self._batch_bget(uncached[split:])
                else:
                    import warnings
                    warnings.warn("Batch contains invalid ids", UserWarning)
        
        # Finally join all the results, but drop all None objects
        entries = filter(lambda e: e is not None, map(bget, ids))
        
        rval = "".join(entries)
        return rval
    
    @cached_method
    def btit(self, ids):
        return KeggApi.btit(self, ids)
    
    @cached_method
    def bconv(self, ids):
        return KeggApi.bconv(self, ids)
    
    ########
    # LinkDB
    ########
    
    @cached_method
    def get_linkdb_by_entry(self, entry_id, db, offset, limit):
       return KeggApi.get_linkdb_by_entry(self, entry_id, db, offset, limit)
        
    @cached_method
    def get_linkdb_between_databases(self, from_db, to_db, offset, limit):
        return KeggApi.get_linkdb_between_databases(self, from_db, to_db, offset, limit)
            
    @cached_method
    def get_genes_by_enzyme(self, enzyme_id, org):
        return KeggApi.get_genes_by_enzyme(self, enzyme_id, org)
    
    @cached_method
    def get_enzymes_by_gene(self, genes_id):
        return KeggApi.get_enzymes_by_gene(self, genes_id)
    
    @cached_method
    def get_enzymes_by_compound(self, compound_id):
        return KeggApi.get_enzymes_by_compound(self, compound_id)
    
    @cached_method
    def get_enzymes_by_glycan(self, glycan_id):
        return KeggApi.get_enzymes_by_glycan(self, glycan_id)
    
    @cached_method
    def get_enzymes_by_reaction(self, reaction_id):
        return KeggApi.get_enzymes_by_reaction(self, reaction_id)
    
    @cached_method
    def get_compounds_by_enzyme(self, enzyme_id):
        return KeggApi.get_compounds_by_enzyme(self, enzyme_id)
    
    @cached_method
    def get_compounds_by_reaction(self, reaction_id):
        return KeggApi.get_compounds_by_reaction(self, reaction_id)
    
    @cached_method
    def get_glycans_by_enzyme(self, enzyme_id):
        return KeggApi.get_glycans_by_enzyme(self, enzyme_id)
    
    @cached_method
    def get_glycans_by_reaction(self, reaction_id):
        return KeggApi.get_glycans_by_reaction(self, reaction_id)
    
    @cached_method
    def get_reactions_by_enzyme(self, enzyme_id):
        return KeggApi.get_reactions_by_enzyme(self, enzyme_id)
    
    @cached_method
    def get_reactions_by_compound(self, compound_id):
        return KeggApi.get_reactions_by_compound(self, compound_id)
    
    @cached_method
    def get_reactions_by_glycan(self, glycan_id):
        return KeggApi.get_reactions_by_glycan(self, glycan_id)
    
    ######
    # SSDB
    ######
    
    @cached_method
    def get_best_best_neighbors_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_best_best_neighbors_by_gene(self, genes_id, offset, limit)
    
    @cached_method
    def get_best_neighbors_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_best_neighbors_by_gene(self, genes_id, offset, limit)
    
    @cached_method
    def get_reverse_best_neighbors_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_reverse_best_neighbors_by_gene(self, genes_id, offset, limit)
    
    @cached_method
    def get_paralogs_by_gene(self, genes_id, offset, limit):
        return KeggApi.get_paralogs_by_gene(self, genes_id, offset, limit)
    
    #######
    # Motif
    #######
    
    @cached_method
    def get_motifs_by_gene(self, genes_id, db):
        return KeggApi.get_motifs_by_gene(self, genes_id, db)
    
    @cached_method
    def get_genes_by_motifs(self, motif_id_list, offset, limit):
        return KeggApi.get_genes_by_motifs(self, motif_id_list, offset, limit)

    ####
    # KO
    ####
    
    @cached_method
    def get_ko_by_gene(self, genes_id):
        return KeggApi.get_ko_by_gene(self, genes_id)
    
    @cached_method
    def get_ko_by_ko_class(self, ko_class_id):
        return KeggApi.service.get_ko_by_ko_class(self, ko_class_id)
    
    @cached_method
    def get_genes_by_ko_class(self, ko_class_id, org, offset, limit):
        return KeggApi.get_genes_by_ko_class(self, ko_class_id, org, offset, limit)
    
    @cached_method
    def get_genes_by_ko(self, ko_id, org):
        return KeggApi.get_genes_by_ko(self, ko_id, org)
    
    #########
    # Pathway
    #########
    
    # TODO
    
    
    
    @cached_method
    def get_genes_by_organism(self, organism, offset=None, limit=None):
        return KeggApi.get_genes_by_organism(self, organism, offset=offset, limit=limit)
    
    @cached_method
    def get_number_of_genes_by_organism(self, organism):
        return KeggApi.get_number_of_genes_by_organism(self, organism)
     
    @cached_method
    def get_pathways_by_genes(self, gene_list):
        return KeggApi.get_pathways_by_genes(self, gene_list)
    
    @cached_method
    def get_pathways_by_enzymes(self, enzyme_list):
        return KeggApi.get_pathways_by_enzymes(self, enzyme_list)
    
    @cached_method
    def get_pathways_by_compounds(self, compound_list):
        return KeggApi.get_pathways_by_compounds(self, compound_list)
    
    @cached_method
    def get_pathways_by_drugs(self, drug_list):
        return KeggApi.get_pathways_by_drugs(self, drug_list)
    
    @cached_method
    def get_pathways_by_glycans(self, glycan_list):
        return KeggApi.get_pathways_by_glycans(self, glycan_list)
    
    @cached_method
    def get_pathways_by_reactions(self, reaction_list):
        return KeggApi.get_pathways_by_reactions(self, reaction_list)
    
    @cached_method
    def get_pathways_by_kos(self, ko_list):
        return KeggApi.get_pathways_by_kos(self, ko_list)
    
    @cached_method
    def get_elements_by_pathway(self, pathway_id):
        return KeggApi.get_elements_by_pathway(self, pathway_id)
    
    @cached_method
    def get_genes_by_pathway(self, pathway_id):
        return KeggApi.get_genes_by_pathway(self, pathway_id)
    
    @cached_method
    def get_enzymes_by_pathway(self, pathway_id):
        return KeggApi.get_enzymes_by_pathway(self, pathway_id)
    
    @cached_method
    def get_compounds_by_pathway(self, pathway_id):
        return KeggApi.get_compounds_by_pathway(self, pathway_id)
    
    @cached_method
    def get_drugs_by_pathway(self, pathway_id):
        return KeggApi.get_drugs_by_pathway(self, pathway_id)
    
    @cached_method
    def get_glycans_by_pathway(self, pathway_id):
        return KeggApi.get_glycans_by_pathway(self, pathway_id)
    
    @cached_method
    def get_reactions_by_pathway(self, pathway_id):
        return KeggApi.get_reactions_by_pathway(self, pathway_id)
    
    @cached_method
    def get_kos_by_pathway(self, pathway_id):
        return KeggApi.get_kos_by_pathway(self, pathway_id)
    
