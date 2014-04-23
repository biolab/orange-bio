from collections import namedtuple

from PyQt4.QtCore import QTimer, QThread, pyqtSlot as Slot

import Orange.data
import Orange.feature
import Orange.network
from Orange.orng.orngDataCaching import data_hints

from Orange.OrangeWidgets import OWWidget
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets import OWItemModels
from Orange.OrangeWidgets.OWConcurrent import ThreadExecutor, Task, methodinvoke

from .. import ppi, taxonomy, gene

NAME = "Gene Network"
DESCRIPTION = "Extract a gene network for a set of genes."
ICON = "icons/GeneNetwork.svg"

INPUTS = [("Data", Orange.data.Table, "set_data")]
OUTPUTS = [("Network", Orange.network.Graph),
           ("Items", Orange.data.Table)]


TAX_MAPPING_BIOGRID = {
    "4530": "39947",
    "4932": "559292",
    "562": "511145",
    "5833": "36329",
    "2104": None,
    "4754": None,
    "31033": None
}

Source = namedtuple(
    "Source",
    ["name", "constructor", "tax_mapping", "sf_domain", "sf_filename"]
)

SOURCES = [
    Source("BioGRID", ppi.BioGRID, ppi.BioGRID.TAXID_MAP,
           "PPI", ppi.BioGRID.SERVER_FILE),
    Source("STRING", ppi.STRING, ppi.STRING.TAXID_MAP,
           "PPI", "string.protein.{taxid}.sqlite")
]


class OWGeneNetwork(OWWidget.OWWidget):
    settingsList = ["taxid", "use_attr_names", "network_source",
                    "include_neighborhood"]
    contextHandlers = {
        "": OWWidget.DomainContextHandler(
            "", ["taxid", "gene_var_index", "use_attr_names"]
        )
    }

    def __init__(self, parent=None, signalManager=None, title="Gene Network"):
        super(OWGeneNetwork, self).__init__(parent, signalManager, title,
                                            wantMainArea=False)

        self.taxid = "9606"
        self.gene_var_index = -1
        self.use_attr_names = False
        self.network_source = 0
        self.include_neighborhood = True
        self.autocommit = False

        self.loadSettings()

        self.taxids = taxonomy.common_taxids()
        self.current_taxid_index = self.taxids.index(self.taxid)

        self.data = None
        self.geneinfo = None
        self.nettask = None
        self._invalidated = False

        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info = OWGUI.widgetLabel(box, "No data on input\n")

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
            box, self, "use_attr_names",
            "Use attribute names",
            callback=self._update_query_genes
        )

        box = OWGUI.widgetBox(self.controlArea, "Network")
        OWGUI.comboBox(
            box, self, "network_source",
            items=[s.name for s in SOURCES],
            callback=self._update_source_db
        )
        OWGUI.checkBox(
            box, self, "include_neighborhood",
            "Include immediate gene neighbors",
            callback=self.invalidate
        )
        box = OWGUI.widgetBox(self.controlArea, "Commit")
#         cb = OWGUI.checkBox(box, self, "autocommit")
        OWGUI.button(box, self, "Commit", callback=self.commit_, default=True)

        self.executor = ThreadExecutor()

    def set_data(self, data):
        self.closeContext()
        self.data = data
        if data is not None:
            self.varmodel[:] = string_variables(data.domain)
            taxid = data_hints.get_hint(data, "taxid", default=self.taxid)

            if taxid in self.taxids:
                self.set_organism(self.taxids.index(taxid))

            self.use_attr_names = data_hints.get_hint(
                data, "genesinrows", default=self.use_attr_names
            )

            self.openContext("", data)
            self._update_info()
            self.invalidate()
            self.commit_()
        else:
            self.varmodel[:] = []

    def set_source_db(self, dbindex):
        self.network_source = dbindex
        self.invalidate()

    def set_organism(self, index):
        self.current_taxid_index = index
        self.taxid = self.taxids[index]
        self.invalidate()

    def set_gene_var(self, index):
        self.gene_var_index = index
        self.invalidate()

    def query_genes(self):
        if self.use_attr_names:
            if self.data is not None:
                return [var.name for var in self.data.domain.attributes]
            else:
                return []
        elif self.gene_var_index >= 0:
            var = self.varmodel[self.gene_var_index]
            genes = [str(inst[var]) for inst in self.data
                     if not inst[var].isSpecial()]
            return list(unique(genes))
        else:
            return []

    def invalidate(self):
        self._invalidated = True

        if self.nettask is not None:
            self.nettask.finished.disconnect(self._on_result_ready)
            self.nettask.future().cancel()
            self.nettask = None

        if self.autocommit:
            QTimer.singleShot(10, self._maybe_commit)

    @Slot()
    def _maybe_commit(self):
        if self._invalidated:
            self.commit_()

    def _commit(self):
        query_genes = self.query_genes()
        geneinfo = self.geneinfo.result()
        ppidb = self.ppidb.result()
        source = SOURCES[self.network_source]
        ppidb = source.constructor()

#         ppidb = self._db
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

        gene_var = self.varmodel[self.gene_var_index]
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

    def commit_(self):
        include_neighborhood = self.include_neighborhood
        query_genes = self.query_genes()
        source = SOURCES[self.network_source]
        taxid = self.taxid

        if self.geneinfo is None:
            self.geneinfo = self.executor.submit(gene.NCBIGeneInfo, taxid)

        geneinfo_f = self.geneinfo
        taxmap = source.tax_mapping
        db_taxid = taxmap.get(taxid, taxid)

        if db_taxid is None:
            raise ValueError("invalid taxid")

        def takes_taxid(s):
            return "{taxid}" in source.sf_filename

        if takes_taxid(source):
            constructor = lambda: source.constructor(db_taxid)
        else:
            constructor = source.constructor

        self.nettask = Task(function=
            lambda: get_gene_network(
                constructor(), geneinfo_f.result(),
                db_taxid, query_genes,
                include_neighborhood=include_neighborhood
            )
        )
        self.nettask.finished.connect(self._on_result_ready)
        self.executor.submit(self.nettask)

        self.setBlocking(True)
        self.progressBarInit()
        self._invalidated = False
        self._update_info()

    @Slot()
    def _on_result_ready(self):
        net = self.nettask.result()
        self.progressBarFinished()
        self.setBlocking(False)
        self.send("Network", net)
        self._update_info()

    def _update_source_db(self):
        self.invalidate()

    def _update_organism(self):
        self.taxid = self.taxids[self.current_taxid_index]
        if self.geneinfo is not None:
            self.geneinfo.cancel()
        self.geneinfo = None
        self.invalidate()

    def _update_query_genes(self):
        self.invalidate()

    def _update_info(self):
        if self.data is None:
            self.info.setText("No data on input\n")
        else:
            names = self.query_genes()
            lines = ["%i unique genes on input" % len(set(names))]
            if self.nettask is not None:
                if not self.nettask.future().done():
                    lines.append("Retrieving ...")
                else:
                    net = self.nettask.result()
                    lines.append("%i nodes %i edges" %
                                 (len(net.nodes()), len(net.edges())))
            else:
                lines.append("")

            self.info.setText("\n".join(lines))


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


def fetch_ppidb(ppidb, taxid, progress=None):
    fname = ppidb.SERVER_FILE
    if "{taxid}" in fname:
        fname = fname.format(taxid=taxid)

    sf.localpath_download(
        ppidb.DOMAIN, fname,
        callback=progress, verbose=True
    )


def get_gene_network(ppidb, geneinfo, taxid, query_genes,
                     include_neighborhood=True, progress=None):
    # Normalize the names to ppidb primary keys

    matcher = geneinfo.matcher
    query_genes = zip(query_genes, map(matcher.umatch, query_genes))
    synonyms = ppidb_synonym_mapping(ppidb, taxid)

    query_genes = [(query_gene, geneid,
                    synonyms.get(query_gene, synonyms.get(geneid)))
                    for query_gene, geneid in query_genes]

    query = [(syn[0], query_gene)
             for query_gene, _, syn in query_genes if syn]
    net = extract_network(ppidb, dict(query), geneinfo, include_neighborhood)
    return net


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
         Orange.feature.Discrete("source", values=["false", "true"])],
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
