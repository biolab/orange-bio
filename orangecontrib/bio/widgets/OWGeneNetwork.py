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
OUTPUTS = [("Network", Orange.network.Graph)]

Source = namedtuple(
    "Source",
    ["name", "constructor", "tax_mapping", "sf_domain", "sf_filename",
     "score_filter"]
)

SOURCES = [
    Source("BioGRID", ppi.BioGRID, ppi.BioGRID.TAXID_MAP,
           "PPI", ppi.BioGRID.SERVER_FILE, False),
    Source("STRING", ppi.STRING, ppi.STRING.TAXID_MAP,
           "PPI", ppi.STRING.FILENAME, True)
]


class OWGeneNetwork(OWWidget.OWWidget):
    settingsList = ["taxid", "use_attr_names", "network_source",
                    "include_neighborhood", "min_score"]
    contextHandlers = {
        "": OWWidget.DomainContextHandler(
            "", ["taxid", "gene_var_index", "use_attr_names"]
        )
    }

    def __init__(self, parent=None, signalManager=None, title="Gene Network"):
        super(OWGeneNetwork, self).__init__(
            parent, signalManager, title, wantMainArea=False,
            resizingEnabled=False
        )

        self.taxid = "9606"
        self.gene_var_index = -1
        self.use_attr_names = False
        self.network_source = 1
        self.include_neighborhood = True
        self.autocommit = False
        self.min_score = 0.9
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
            callback=self._on_source_db_changed
        )
        OWGUI.checkBox(
            box, self, "include_neighborhood",
            "Include immediate gene neighbors",
            callback=self.invalidate
        )
        self.score_spin = OWGUI.doubleSpin(
            box, self, "min_score", 0.0, 1.0, step=0.001,
            label="Minimal edge score",
            callback=self.invalidate
        )
        self.score_spin.setEnabled(SOURCES[self.network_source].score_filter)

        box = OWGUI.widgetBox(self.controlArea, "Commit")
        OWGUI.button(box, self, "Commit", callback=self.commit, default=True)

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

            if not (0 <= self.gene_var_index < len(self.varmodel)):
                self.gene_var_index = len(self.varmodel) - 1

            self.openContext("", data)
            self.invalidate()
            self.commit()
        else:
            self.varmodel[:] = []
            self.send("Network", None)

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
            self.commit()

    @Slot()
    def advance(self):
        self.progressBarValue = (self.progressBarValue + 1) % 100

    @Slot(float)
    def set_progress(self, value):
        self.progressBarValue = value

    def commit(self):
        include_neighborhood = self.include_neighborhood
        query_genes = self.query_genes()
        source = SOURCES[self.network_source]
        if source.score_filter:
            min_score = self.min_score
            assert source.name == "STRING"
            min_score = min_score * 1000
        else:
            min_score = None

        taxid = self.taxid
        progress = methodinvoke(self, "advance")
        if self.geneinfo is None:
            self.geneinfo = self.executor.submit(
                fetch_ncbi_geneinfo, taxid, progress
            )

        geneinfo_f = self.geneinfo
        taxmap = source.tax_mapping
        db_taxid = taxmap.get(taxid, taxid)

        if db_taxid is None:
            raise ValueError("invalid taxid for this network")

        def fetch_network():
            geneinfo = geneinfo_f.result()
            ppidb = fetch_ppidb(source, db_taxid, progress)
            return get_gene_network(ppidb, geneinfo, db_taxid, query_genes,
                                    include_neighborhood=include_neighborhood,
                                    min_score=min_score,
                                    progress=methodinvoke(self, "set_progress", (float,)))

        self.nettask = Task(function=fetch_network)
        self.nettask.finished.connect(self._on_result_ready)

        self.executor.submit(self.nettask)

        self.setBlocking(True)
        self.setEnabled(False)
        self.progressBarInit()
        self._invalidated = False
        self._update_info()

    @Slot(object)
    def _on_result_ready(self,):
        self.progressBarFinished()
        self.setBlocking(False)
        self.setEnabled(True)
        net = self.nettask.result()
        self._update_info()
        self.send("Network", net)

    def _on_source_db_changed(self):
        source = SOURCES[self.network_source]
        self.score_spin.setEnabled(source.score_filter)
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


def fetch_ppidb(ppisource, taxid, progress=None):
    fname = ppisource.sf_filename

    if "{taxid}" in fname:
        if taxid in ppisource.tax_mapping:
            taxid_m = ppisource.tax_mapping[taxid]
            if taxid_m is None:
                raise ValueError(taxid)
            taxid = taxid_m
        fname = fname.format(taxid=taxid)
        constructor = lambda: ppisource.constructor(taxid)
    else:
        constructor = ppisource.constructor

    sf.localpath_download(
        ppisource.sf_domain, fname, callback=progress, verbose=True
    )
    return constructor()


def fetch_ncbi_geneinfo(taxid, progress=None):
    taxid = gene.NCBIGeneInfo.TAX_MAP.get(taxid, taxid)
    sf.localpath_download(
        "NCBI_geneinfo", "gene_info.{taxid}.db".format(taxid=taxid),
        callback=progress, verbose=True,
    )
    return gene.NCBIGeneInfo(taxid)


def get_gene_network(ppidb, geneinfo, taxid, query_genes,
                     include_neighborhood=True, min_score=None,
                     progress=None):
    if progress is not None:
        progress(1.0)

    # Normalize the names to ppidb primary keys
    matcher = geneinfo.matcher
    query_genes = zip(query_genes, map(matcher.umatch, query_genes))
    synonyms = ppidb_synonym_mapping(ppidb, taxid)

    query_genes = [(query_gene, geneid,
                    synonyms.get(query_gene, synonyms.get(geneid)))
                    for query_gene, geneid in query_genes]

    query = [(syn[0], query_gene)
             for query_gene, _, syn in query_genes if syn]

    net = extract_network(ppidb, dict(query), geneinfo, include_neighborhood,
                          min_score, progress=progress)
    return net


def extract_network(ppidb, query, geneinfo, include_neighborhood=True,
                    min_score=None, progress=None):
    """
    include neighborhood
    """
    from functools import partial
    from collections import defaultdict
    from itertools import count

    if not isinstance(query, dict):
        query = {name: name for name in query}

    report_weights = True
    if isinstance(ppidb, ppi.BioGRID):
        # BioGRID scores are not comparable (they can be p values,
        # confidence scores, ..., i.e. whatever was reported in the source
        # publication)
        report_weights = False
        if min_score is not None:
            raise ValueError("min_score used with BioGrid")

    graph = Orange.network.Graph()
    # node ids in Orange.network.Graph need to be in [0 .. n-1]
    nodeids = defaultdict(partial(next, count()))

    def gi_info(names):
        mapping = [(name, geneinfo.matcher.umatch(name)) for name in names]
        mapping = [(name, match) for name, match in mapping if match]
        entries = [(name, geneinfo[match]) for name, match in mapping]

        if len(entries) > 1:
            # try to resolve conflicts by prioritizing entries whose
            # symbol/gene_id/locus_tag exactly matches the synonym name.
            entries_ = [(name, entry) for name, entry in entries
                        if name in [entry.gene_id, entry.symbol, entry.locus_tag]]
            if len(entries_) == 1:
                entries = entries_

        if len(entries) == 0:
            return None
        elif len(entries) >= 1:
            # Need to report multiple mappings
            return entries[0][1]

    # Add query nodes.
    for key, query_name in query.items():
        nodeid = nodeids[key]
        synonyms = ppidb.synonyms(key)
        entry = gi_info(synonyms)
        graph.add_node(
            nodeid,
            key=key,
            synonyms=synonyms,
            query_name=query_name,
            symbol=entry.symbol if entry is not None else ""
        )

    if include_neighborhood:
        # extend the set of nodes in the network with immediate neighborers
        edges_iter = (edge for key in query for edge in ppidb.edges(key))
        for id1, id2, score in edges_iter:
            if min_score is None or score >= min_score:
                nodeid1 = nodeids[id1]
                nodeid2 = nodeids[id2]
                if nodeid1 not in graph:
                    synonyms1 = ppidb.synonyms(id1)
                    entry1 = gi_info(synonyms1)
                    symbol1 = entry1.symbol if entry1 is not None else ""
                    graph.add_node(
                        nodeid1, key=id1, synonyms=synonyms1,
                        symbol=symbol1
                    )

                if nodeid2 not in graph:
                    synonyms2 = ppidb.synonyms(id2)
                    entry2 = gi_info(synonyms2)
                    symbol2 = entry2.symbol if entry2 is not None else ""
                    graph.add_node(
                        nodeid2, key=id2, synonyms=synonyms2,
                        symbol=symbol2
                    )

    # add edges between nodes
    for i, id1 in enumerate(nodeids.keys()):
        if progress is not None:
            progress(100.0 * i / len(nodeids))

        for _, id2, score in ppidb.edges(id1):
            if id2 in nodeids and (min_score is None or score >= min_score):
                nodeid1 = nodeids[id1]
                nodeid2 = nodeids[id2]
                assert nodeid1 in graph and nodeid2 in graph
                if score is not None and report_weights:
                    graph.add_edge(nodeid1, nodeid2, weight=score)
                else:
                    graph.add_edge(nodeid1, nodeid2)

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
          str(node.get("symbol", nodeid)),
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
