"""\
==============================================
KEGG - Kyoto Encyclopedia of Genes and Genomes
==============================================

This is a python module for access to `KEGG`_ using its web services.

To use this module you need to have `slumber` and `requests` package
installed.

.. _`KEGG`: http://www.genome.jp/kegg/


"""
from __future__ import absolute_import

import os
import sys
import urllib2

from collections import defaultdict
from itertools import chain
from datetime import datetime

from Orange.utils import lru_cache
from Orange.utils import progress_bar_milestones
from Orange.utils import (
    deprecated_keywords, deprecated_attribute, deprecated_function_name
)

from .. import obiProb

from . import databases
from . import entry

from .brite import BriteEntry, Brite

from . import api
from . import conf
from . import pathway

KEGGGenome = databases.Genome
KEGGGenes = databases.Genes
KEGGEnzyme = databases.Enzyme
KEGGReaction = databases.Reaction
KEGGPathways = databases.Pathway
KEGGCompound = databases.Compound

KEGGBrite = Brite
KEGGBriteEntry = BriteEntry

KEGGPathway = pathway.Pathway

DEFAULT_CACHE_DIR = conf.params["cache.path"]


class OrganismNotFoundError(Exception):
    pass


class Organism(object):
    """
    A convenience class for retrieving information regarding an
    organism in the KEGG Genes database.

    :param org: KEGGG organism code (e.g. "hsa", "sce")
    :type org: str

    """
    def __init__(self, org, genematcher=None):
        self.org_code = self.organism_name_search(org)
        self.genematcher = genematcher
        self.api = api.CachedKeggApi()

    @property
    def org(self):
        """
        KEGG organism code.
        """
        return self.org_code

    @property
    def genes(self):
        """
        An :class:`Genes` database instance for this organism.
        """
        # TODO: This should not be a property but a method.
        # I think it was only put here as back compatibility with old obiKEGG.
        if not hasattr(self, "_genes"):
            genes = KEGGGenes(self.org_code)
            self._genes = genes
        return self._genes

    def gene_aliases(self):
        """
        Return a list of sets of equal genes (synonyms) in KEGG for
        this organism.

        .. note::

            This only includes 'ncbi-geneid' and 'ncbi-gi' records
            from the KEGG Genes DBLINKS entries.

        """
        definitions = self.api.list(self.org_code)
        ncbi_geneid = self.api.conv(self.org_code, "ncbi-geneid")
        ncbi_gi = self.api.conv(self.org_code, "ncbi-gi")

        aliases = defaultdict(set)

        for entry_id, definition in definitions:
            # genes entry id without the organism code
            aliases[entry_id].add(entry_id.split(":", 1)[1])
            # all names in the NAME field (KEGG API list returns
            # 'NAME; DEFINITION') fields for genes
            names = definition.split(";")[0].split(",")
            aliases[entry_id].update([name.strip() for name in names])

        for source_id, target_id in chain(ncbi_geneid, ncbi_gi):
            aliases[target_id].add(source_id.split(":", 1)[1])

        return [set([entry_id]).union(names)
                for entry_id, names in aliases.iteritems()]

    def pathways(self, with_ids=None):
        """
        Return a list of all pathways for this organism.
        """
        if with_ids is not None:
            return self.api.get_pathways_by_genes(with_ids)
        else:
            return [p.entry_id for p in self.api.list_pathways(self.org_code)]

    def list_pathways(self):
        """
        List all pathways.
        """
        # NOTE: remove/deprecate and use pathways()
        return self.pathways()

    def get_linked_pathways(self, pathway_id):
        self.api.get_linked_pathways(pathway_id)

    def enzymes(self, genes=None):
        raise NotImplementedError()

    def get_enriched_pathways(self, genes, reference=None,
                              prob=obiProb.Binomial(), callback=None):
        """
        Return a dictionary with enriched pathways ids as keys
        and (list_of_genes, p_value, num_of_reference_genes) tuples
        as items.

        """
        if reference is None:
            reference = self.genes.keys()
        reference = set(reference)

        allPathways = defaultdict(lambda: [[], 1.0, []])
        milestones = progress_bar_milestones(len(genes), 100)
        pathways_db = KEGGPathways()

        pathways_for_gene = []
        for i, gene in enumerate(genes):
            pathways_for_gene.append(self.pathways([gene]))
            if callback and i in milestones:
                callback(i * 50.0 / len(genes))

        # pre-cache for speed
        pathways_db.pre_cache([pid for pfg in pathways_for_gene
                               for pid in pfg])
        for i, (gene, pathways) in enumerate(zip(genes, pathways_for_gene)):
            for pathway in pathways:
                if pathways_db.get_entry(pathway).gene:
                    allPathways[pathway][0].append(gene)
            if callback and i in milestones:
                callback(50.0 + i * 50.0 / len(genes))

        pItems = allPathways.items()

        for i, (p_id, entry) in enumerate(pItems):
            pathway = pathways_db.get_entry(p_id)
            entry[2].extend(reference.intersection(pathway.gene or []))
            entry[1] = prob.p_value(len(entry[0]), len(reference),
                                    len(entry[2]), len(genes))
        return dict([(pid, (genes, p, len(ref)))
                     for pid, (genes, p, ref) in allPathways.items()])

    def get_genes_by_enzyme(self, enzyme):
        enzyme = KEGGEnzyme().get_entry(enzyme)
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
        pathways = [self.genes[id].pathway for id in gene_ids
                    if self.genes[id].pathway]
        pathways = reduce(set.union, pathways, set())
        return [id for id in pathways
                if gene_ids.issubset(KEGGPathway(id).genes())]

    def get_pathways_by_enzymes(self, enzyme_ids):
        enzyme_ids = set(enzyme_ids)
        pathways = [KEGGEnzyme()[id].pathway for id in enzyme_ids]
        pathways = reduce(set.union, pathways, set())
        return [id for id in pathways
                if enzyme_ids.issubset(KEGGPathway(id).enzymes())]

    def get_pathways_by_compounds(self, compound_ids):
        compound_ids = set(compound_ids)
        pathways = [KEGGCompound()[id].pathway for id in compound_ids]
        pathways = reduce(set.union, pathways, set())
        return [id for id in pathways
                if compound_ids.issubset(KEGGPathway(id).compounds())]

    def get_enzymes_by_compound(self, compound_id):
        return KEGGCompound()[compound_id].enzyme

    def get_enzymes_by_gene(self, gene_id):
        return self.genes[gene_id].enzymes

    def get_compounds_by_enzyme(self, enzyme_id):
        return self._enzymes_to_compounds.get(enzyme_id)

    @deprecated_keywords({"caseSensitive": "case_sensitive"})
    def get_unique_gene_ids(self, genes, case_sensitive=True):
        """
        Return a tuple with three elements. The first is a dictionary
        mapping from unique geneids to gene names in genes, the second
        is a list of conflicting gene names and the third is a list of
        unknown genes.

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

    @deprecated_function_name
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
            return genome[name].organism_code
        except KeyError:
            raise OrganismNotFoundError(name)

    @classmethod
    def organism_version(cls, name):
        name = cls.organism_name_search(name)
        genome = KEGGGenome()
        info = genome.api.info(name)
        return info.release

    def _set_genematcher(self, genematcher):
        setattr(self, "_genematcher", genematcher)

    def _get_genematcher(self):
        if getattr(self, "_genematcher", None) is None:
            from .. import obiGene
            if self.org_code == "ddi":
                self._genematcher = obiGene.matcher(
                    [obiGene.GMKEGG(self.org_code), obiGene.GMDicty(),
                     [obiGene.GMKEGG(self.org_code), obiGene.GMDicty()]]
                )
            else:
                self._genematcher = obiGene.matcher(
                    [obiGene.GMKEGG(self.org_code)])

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
    for r in res:
        e = genome[r]

        if e.taxid in [taxid, genome.TAXID_MAP.get(taxid, taxid)]:
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


def main():
    KEGGGenome()
    import doctest
    extraglobs = {"api": KeggApi()}
    doctest.testmod(optionflags=doctest.ELLIPSIS, extraglobs=extraglobs)


if __name__ == "__main__":
    sys.exit(main())
