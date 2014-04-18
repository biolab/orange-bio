from collections import namedtuple

from PyQt4.QtCore import QTimer, pyqtSlot as Slot

import Orange.data
import Orange.network
import Orange.feature

from Orange.OrangeWidgets import OWWidget
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets import OWItemModels
from Orange.OrangeWidgets.OWConcurrent import ThreadExecutor, Future

from .. import ppi, taxonomy, gene

NAME = "Gene Network"
DESCRIPTION = "Extract a gene network for a set of genes."
ICON = "icons/GeneNetwork.svg"

INPUTS = [("Data", Orange.data.Table, "set_data")]
OUTPUTS = [("Network", Orange.network.Graph),
           ("Items", Orange.data.Table)]


Source = namedtuple(
    "Source",
    ["name", "constructor", "sf_domain", "sf_filename"]
)

SOURCES = [
    Source("BioGRID", ppi.BioGRID, "PPI", ppi.BioGRID.SERVER_FILE),
#     Source("STRING", ppi.STRING, "PPI", ppi.STRING.FILENAME)
]


class OWGeneNetwork(OWWidget.OWWidget):
    settingsList = ["taxid", "genes_in_row", "network_source",
                    "include_neighborhood"]
    contextHandlers = {
        "": OWWidget.DomainContextHandler(
            "", ["taxid", "gene_var_index", "genes_in_row"]
        )
    }

    def __init__(self, parent=None, signalManager=None, title="Gene Network"):
        super(OWGeneNetwork, self).__init__(parent, signalManager, title,
                                            wantMainArea=False)

        self.taxid = "9606"
        self.gene_var_index = -1
        self.genes_in_rows = False
        self.network_source = 0
        self.include_neighborhood = True
        self.autocommit = False

        self.loadSettings()

        self.taxids = taxonomy.common_taxids()
        self.current_taxid_index = self.taxids.index(self.taxid)

        self.db = None
        self.data = None
        self.geneinfo = None
        self._invalidated = False

        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info = OWGUI.widgetLabel(box, "No input\n")

        box = OWGUI.widgetBox(self.controlArea, "Network Source")
        OWGUI.comboBox(
            box, self, "network_source",
            items=[s.name for s in SOURCES],
            callback=self._update_source_db
        )
        box = OWGUI.widgetBox(self.controlArea, "Organism")
        self.organism_cb = OWGUI.comboBox(
            box, self, "current_taxid_index",
            items=map(taxonomy.name, self.taxids),
            callback=self._update_organism
        )
        box = OWGUI.widgetBox(self.controlArea, "Genes")
        self.genes_cb = OWGUI.comboBox(
            box, self, "gene_var_index", callback=self._update_query_genes
        )
        self.varmodel = OWItemModels.VariableListModel()
        self.genes_cb.setModel(self.varmodel)

        OWGUI.checkBox(
            box, self, "genes_in_rows",
            "Use attribute names",
            callback=self._update_query_genes
        )
        box = OWGUI.widgetBox(self.controlArea, "Neighborhood")
        OWGUI.checkBox(
            box, self, "include_neighborhood",
            "Include immediate neighborers",
            callback=self.invalidate
        )
        box = OWGUI.widgetBox(self.controlArea, "Commit")
#         cb = OWGUI.checkBox(box, self, "autocommit")
        OWGUI.button(box, self, "Commit", callback=self.commit, default=True)

        self.executor = ThreadExecutor()

#         QTimer.singleShot(0, self.initialize)
        self._update_source_db()
        self._update_organism()

    @Slot()
    def initialize(self):
        if self.db is None:
            self.db = ppi.BioGRID()

    def set_data(self, data):
        self.closeContext()
        self.data = data
        if data is not None:
            self.varmodel[:] = string_variables(data.domain)
            self.openContext("", data)
#             self._update_info()
            self.commit()
        else:
            self.varmodel[:] = []

    def set_source_db(self, dbindex):
        self.network_source = dbindex
        # by taxid??
        self.db = SOURCES[self.network_source]()
        self.invalidate()

    def set_organism(self, index):
        self.current_taxid_index = index
        self.taxid = self.taxids[index]
        self.invalidate()

    def set_gene_var(self, index):
        self.gene_var_index = index
        self.invalidate()

    def query_genes(self):
        if self.gene_var_index >= 0:
            var = self.varmodel[self.gene_var_index]
            genes = [unicode(inst[var]) for inst in self.data
                     if not inst[var].isSpecial()]
            return list(unique(genes))
        else:
            return []

    def invalidate(self):
        self._invalidated = True
        if self.autocommit:
            QTimer.singleShot(0, self._maybe_commit)

    @Slot()
    def _maybe_commit(self):
        if self._invalidated:
            self.commit()

    def commit(self):
        gene_var = self.varmodel[self.gene_var_index]
        query_genes = self.query_genes()
        geneinfo = self.geneinfo.result()

        ppidb = self._db
        tax_map = self._tax_mapping
        taxid = tax_map.get(self.taxid, None)
        if taxid is None:
            raise ValueError

        # Normalize the names to ppidb primary keys
        matcher = geneinfo.matcher
        query_genes = zip(query_genes, map(matcher.umatch, query_genes))
        synonyms = ppidb_synonym_mapping(ppidb, taxid)

        query_genes = [(query_gene, geneid, synonyms.get(geneid, None))
                        for query_gene, geneid in query_genes]

        query = [(syn[0], query_gene)
                 for query_gene, _, syn in query_genes if syn]
        net = extract_network(ppidb, dict(query), geneinfo,
                              self.include_neighborhood)

        # Items has the same domain as the input with query instances
        # corresponding to nodes where applicable (if include_neighborhood
        # then nodes not contained in the query will have instances with all
        # unknown values

        data_by_query_name = {str(inst[gene_var]): inst
                              for inst in self.data}
        items = Orange.data.Table(self.data.domain)
        for nodeid, node in sorted(net.nodes(data=True)):
            if "query_name" in node:
                items.append(data_by_query_name[node["query_name"]])
            else:
                items.append(Orange.data.Instance(self.data.domain))
        self.send("Network", net)
        self.send("Items", items)
        self._invalidated = False

    def fetch_db(self, which):
        pass

    def _update_source_db(self):
        source = ppi.BioGRID()
        orgs = map(str, source.organisms())

        tax_mapping = taxonomy_match(orgs, self.taxids)
        tax_mapping = {common_tax: db_tax
                       for db_tax, common_tax in tax_mapping.items()
                       if common_tax is not None}

        # TODO: Disable the items in Organism combo box that are not
        # in this list.
#         model = self.organism_cb.model()
#         for i, taxid in self.taxids:
#             enabled = taxid in tax_mapping

        self._db = source
        self._tax_mapping = tax_mapping
        self.invalidate()

    def _update_organism(self):
        self.taxid = self.taxids[self.current_taxid_index]
        if self.geneinfo is not None:
            self.geneinfo.cancel()

        self.geneinfo = self.executor.submit(
            lambda: gene.NCBIGeneInfo(self.taxid)
        )
        self.invalidate()

    def _update_query_genes(self):
        self.invalidate()


def unique(seq):
    seen = set()
    for el in seq:
        if el not in seen:
            seen.add(el)
            yield el


def string_variables(domain):
    variables = domain.variables + domain.getmetas().values()
    return [v for v in variables if isinstance(v, Orange.feature.String)]


def multimap_inverse(multimap):
    """
    Return a multimap inverse relation.
    """
    d = {}
    for key, values in multimap.iteritems():
        for v in values:
            d.setdefault(v, []).append(key)
    return d


def ppidb_synonym_mapping(ppidb, taxid):
    keys = ppidb.ids(taxid)
    mapping = {key: ppidb.synonyms(key) for key in keys}
    return multimap_inverse(mapping)


def taxonomy_match(query_taxids, target_taxids):
    taxid_mapping = {}
    target_taxids = set(target_taxids)

    for taxid in query_taxids:
        mapped = taxid_map(taxid, target_taxids)
        taxid_mapping[taxid] = mapped

    return taxid_mapping


def taxid_map(query, targets):
    if query in targets:
        return query

    lineage = taxonomy.lineage(query)
    if any(tid in targets for tid in lineage):
        return set(lineage).intersection(targets).pop()
    else:
        return None

from Orange.utils import serverfiles as sf


def fetch_ppidb(ppidb, progress=None):
    sf.localpath_download(
        ppidb.DOMAIN, ppidb.SERVER_FILE,
        callback=progress, verbose=True
    )
    db = ppidb()

    mapping = taxonomy_match(db.organisms(), taxonomy.common_taxids())
    return db, mapping


def extract_network(ppidb, query, geneinfo, include_neighborhood=True):
    """
    include neighborhood
    """
    from functools import partial
    from collections import defaultdict
    from itertools import count

    if not isinstance(query, dict):
        query = {name: name for name in query}

    graph = Orange.network.Graph()
    # node ids in Orange.network.Graph need to be in [0 .. n-1]
    nodeids = defaultdict(partial(next, count()))

    for key, query_name in query.items():
        nodeid = nodeids[key]
        graph.add_node(
            nodeid,
            key=key,
            synonyms=ppidb.synonyms(key),
            query_name=query_name
        )

    for key in query:
        edges = ppidb.edges(key)

        for id1, id2, score in edges:
            if include_neighborhood or id1 in query and id2 in query:
                nodeid1 = nodeids[id1]
                nodeid2 = nodeids[id2]
                if nodeid1 not in graph:
                    graph.add_node(
                        nodeid1, key=id1, synonyms=ppidb.synonyms(id1))
                if nodeid2 not in graph:
                    graph.add_node(
                        nodeid2, key=id1, synonyms=ppidb.synonyms(id2))

                if score is not None:
                    graph.add_edge(nodeid1, nodeid2, weight=score)
                else:
                    graph.add_edge(nodeid1, nodeid2)

    # The "items" should table should contain the query name (if
    nodedomain = Orange.data.Domain(
        [Orange.feature.String("Query name"),  # if applicable
         Orange.feature.String("id"),          # ppidb primary key
         Orange.feature.String("Synonyms"),    # ppidb synonyms
         Orange.feature.String("Symbol"),      # ncbi gene name ??
         Orange.feature.Discrete("source", values=["false", "true"])
         ],
        None
    )

    node_items = sorted(graph.node.items(), key=lambda t: nodeids[t[0]])

    nodeitems = Orange.data.Table(
        nodedomain,
        [[str(node.get("query_name", "")),
          str(node.get("key", "")),
          str(", ".join(node.get("synonyms", []))),
          str(node.get("query_name", nodeid)),
          "true" if "query_name" in node else "false"]
         for nodeid, node in node_items]
    )
    graph.set_items(nodeitems)
    return graph


def main():
    from PyQt4.QtGui import QApplication
    app = QApplication([])
    w = OWGeneNetwork()
    brown = Orange.data.Table("brown-selected")
    w.set_data(Orange.data.Table(brown[:5]))
    w.show()
    app.exec_()
    w.saveSettings()
    return 0

if __name__ == "__main__":
    main()
