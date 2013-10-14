"""
<name>Gene Atlas Tissue Expression</name>
<description></description>
<icon>icons/GeneExpressionAtlas.svg</icon>
"""

from __future__ import absolute_import

import os, sys
from collections import defaultdict
import shelve

import Orange
from Orange.utils import lru_cache
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from Orange.orng import orngDataCaching, orngServerFiles

from .. import obiArrayExpress, obiGene, obiGeneAtlas

NAME = "Gene Atlas Tissue Expression"
DESCRIPTION = ""
ICON = "icons/GeneExpressionAtlas.svg"
PRIORITY = 5000

INPUTS = [("Example Table", Orange.data.Table, "set_data")]
OUTPUTS = [("Selected Genes", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWGeneAtlasTissueExpression.OWGeneAtlasTissueExpression"]


TAXID_TO_ORG = obiGeneAtlas.TAXID_TO_ORG

class OWGeneAtlasTissueExpression(OWWidget):
    contextHandlers = {"": DomainContextHandler("", ["selected_organism",
                                                     "selected_gene_attr",
                                                     "genes_in_columns",
                                                     "selected_ef",
                                                     "selected_ef_value"])}
    settingsList = ["selected_organism", "selected_ef", "selected_ef_value"]
    
    ORGANISMS = obiArrayExpress.ATLAS_ORGANISMS
    FACTORS = ["Organism part", "Disease state", "Cell type"]
    
    def __init__(self, parent=None, signalManager=None, title="Gene Atlas Tissue Expression"):
        OWWidget.__init__(self, parent, signalManager, title)
        
        self.inputs = [("Example Table", Orange.data.Table, self.set_data)]
        self.outputs = [("Selected Genes", Orange.data.Table)]
        
        self.selected_organism = "Homo sapiens"
        self.selected_gene_attr = 0
        self.genes_in_columns = False
        self.selected_ef = 0
        self.selected_ef_value = 0
        
        self.loadSettings()
        
        #####
        # GUI
        #####
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.info_label = OWGUI.widgetLabel(box, "No data on input.\n")
        
        box = OWGUI.widgetBox(self.controlArea, "Organism", addSpace=True)
        cb = OWGUI.comboBox(box, self, "selected_organism",
                            items=self.ORGANISMS,
                            tooltip="Organism name",
                            callback=self.on_organism_change,
                            sendSelectedValue=True,
                            valueType=str
                            )
        cb.setMaximumWidth(250)
        
        box = OWGUI.widgetBox(self.controlArea, "Gene Attribute", addSpace=True)
        self.gene_attr_cb = OWGUI.comboBox(box, self, "selected_gene_attr",
                              tooltip="Attribute (column) containing the gene names.",
                              callback=self.on_gene_attr_change,
                              )
        self.gene_attr_cb.setMaximumWidth(250)
        
        cb = OWGUI.checkBox(box, self, "genes_in_columns", "Use attribute names",
                       tooltip="Gene names in columns.",
                       callback=self.on_genes_change,)
        cb.disables.append((-1, self.gene_attr_cb))
        cb.makeConsistent()
        
        box = OWGUI.widgetBox(self.controlArea, "Tissues", addSpace=True)
        self.categories_cb = OWGUI.comboBox(box, self, "selected_ef",
                                box="Categories", 
                                items=self.FACTORS,
                                tooltip="Experimental factor.",
                                callback=self.on_ef_change,
                                )
        self.categories_cb.box.setFlat(True)
        
        self.values_cb = OWGUI.comboBox(box, self, "selected_ef_value",
                                box="Values",
                                tooltip="Experimental factor value.",
                                callback=self.on_ef_value_change
                                )
        self.values_cb.setMaximumWidth(250)
        self.values_cb.box.setFlat(True)
        
        box = OWGUI.widgetBox(self.controlArea, "Cache", addSpace=True)
        OWGUI.button(box, self, "Clear cache",
                     callback=self.on_cache_clear,
                     tooltip="Clear Gene Atlas cache.")
        
        OWGUI.rubber(self.controlArea)
        
        OWGUI.button(self.controlArea, self, label="Commit",
                     callback=self.commit,
                     tooltip="Send selected genes")
        
        self.report_view = QTreeView(self.mainArea)
        self.report_view.setSelectionMode(QTreeView.ExtendedSelection)
        self.report_view.setSortingEnabled(True)
        self.report_view.setRootIsDecorated(False)
        self.report_view.setAlternatingRowColors(True)
        self.report_view.setEditTriggers(QTreeView.NoEditTriggers)
        self.mainArea.layout().addWidget(self.report_view)
        self.report_header = ["Gene symbol", "Up", "Down"]
        
        model = QStandardItemModel()
        model.setHorizontalHeaderLabels(self.report_header)
        self.report_view.setModel(model)
        
        self.data = None
        self.candidate_vars = []
        self.candidate_var_names = []
        self.results = {}, {}, {}
        
        self.ensembl_info = None
        self.gene_matcher = obiGene.GMDirect()
        self.loaded_matcher_taxid = None
        self.unknown_genes = []
        self.query_genes = []
        
#        self.set_organism(self.selected_organism, update_results=False)
        
        self.get_atlas_summary = obiGeneAtlas.get_atlas_summary
        
        #Cached construct_matcher
        @lru_cache(maxsize=3)
        def my_cached_matcher(org):
            return obiGeneAtlas.default_gene_matcher(org)
        self.construct_matcher = my_cached_matcher
        
    def set_data(self, data=None):
        """ Set the input example table with gene names. 
        """
        self.closeContext("")
        self.clear()
        self.data = data
        if data is not None:
            self.init_gene_box(data)
            taxid = orngDataCaching.data_hints.get_hint(data, key="taxid", default=None)
            if taxid is not None:
                if taxid in TAXID_TO_ORG:
                    self.selected_organism = TAXID_TO_ORG[taxid]
                    
            genesinrows = orngDataCaching.data_hints.get_hint(data, key="genesinrows", default=None)
            if genesinrows is not None:
                self.genes_in_columns = genesinrows
                
            self.openContext("", data)
            self.update_gene_matcher()
            self.update_info_box()
            self.run_query()
            self.update_ef_values_box()
            self.update_report_view()
        else:
            pass
            
    def clear(self):
        """ Clear the widget state.
        """
        self.data = None
        self.gene_attr_cb.clear()
        self.values_cb.clear()
        self.selected_gene_attr = 0
        self.selected_ef_value = 0
        self.candidate_vars = []
        self.candidate_var_names = []
        self.results = {}, {}, {}
        model = QStandardItemModel()
        model.setHorizontalHeaderLabels(self.report_header)
        self.report_view.setModel(model)

    def init_gene_box(self, data):
        """ Populate the 'Gene Attribute' box
        """
        vars = data.domain.variables + data.domain.get_metas().values()
        vars = [var for var in vars if isinstance(var, (Orange.feature.String))]
        self.candidate_vars = vars
        self.gene_attr_cb.addItems([var.name for var in vars])
        # See if any var name contains 'gene'
        gene_in_name = [v for v in vars if "gene" in v.name.lower()]
        if gene_in_name:
            self.selected_gene_attr = vars.index(gene_in_name[-1])
        self.candidate_var_names = [a.name for a in data.domain.attributes]
        if not self.candidate_vars:
            self.genes_in_columns = True            
            
    def on_organism_change(self):
        """ Called when the user changes the organism.
        """
        self.update_gene_matcher()
        self.run_query()
        self.update_ef_values_box()
        self.update_report_view()
    
    def on_gene_attr_change(self):
        """ Called  when the user changes source gene attribute.
        """
        self.update_info_box()
        self.run_query()
        self.update_ef_values_box()
        self.update_report_view()
    
    def on_genes_change(self):
        """ Called when the user toogles the "Genes in columns' check box. 
        """
        self.update_info_box()
        self.run_query()
        self.update_ef_values_box()
        self.update_report_view()
    
    def on_ef_change(self):
        """ Called when the user changes the EF factor.
        """
        self.update_ef_values_box()
        # TODO: restore previous selected value for factor (if value in ef_values)
        self.update_report_view()
        
    def on_ef_value_change(self):
        """ Called when the user changes the EF value.
        """
        results = self.results[self.selected_ef]
        ef_value = self.ef_values[self.selected_ef_value]
        self.update_report_view()
        
    def on_cache_clear(self):
        from contextlib import closing
        with closing(obiGeneAtlas._cache()) as cache:
            cache.clear()
        
    def input_genes(self, full_res=False):
        """ Extract input gene names from ``self.data``. 
        """
        if self.genes_in_columns:
            names = self.candidate_var_names
        else:
            names = []
            attr = self.candidate_vars[self.selected_gene_attr]
            for ex in self.data:
                if not ex[attr].is_special():
                    names.append(str(ex[attr]))
            
        return sorted(set(names))
    
    def match_genes(self, genes, full_res=True):
        """ Match the genes to ensembl gene ids.
        """
        if self.gene_matcher and self.ensembl_info:
            self.gene_matcher.set_targets(self.ensembl_info.keys())
            match_genes = map(self.gene_matcher.match, genes)
            genes = map(self.gene_matcher.umatch, genes)
        else:
            match_genes = None
        genes = sorted(set(genes) - set([None]))
        
        if full_res:
            return genes, match_genes
        else:
            return genes
    
    def run_query(self):
        """ Query the Gene Atlas.
        """
        if self.data is not None:
            genes = self.input_genes()
            matched, match_res = self.match_genes(genes, True)
            all_matched, unknown = [], []
            for gene, match in zip(genes, match_res):
                if len(match) == 1:
                    all_matched.append(match[0])
                elif len(match) > 1:
                    all_matched.append(match[0])
                else:
                    unknown.append(gene)
            self.unknown_genes = unknown
            self.warning(2)
            self.information(2)
            if float(len(unknown)) / len(genes) > 0.5:
                self.warning(2, "%i unknown genes" % len(unknown))
            elif unknown:
                self.information(2, "%i unknown genes" % len(unknown))
                
            gm = obiGene.GMDirect()
            gm.set_targets(all_matched)
            
            # Non threaded
#            self.results = self.get_atlas_summary(all_matched, self.selected_organism, gm)
            # Threaded
            self.error(0)
            self.update_info_box(query_running=True)
            self.progressBarInit()
            self.controlArea.setEnabled(False)
            self.results = None
            try:
                call = self.asyncCall(self.get_atlas_summary, (all_matched,
                                                               self.selected_organism,
                                                               gm),
                                      name="Query Gene Expression Atlas",
                                      onError=self.handle_assync_error)
                QObject.connect(call, SIGNAL("progressChanged(float)"), self.progressBarSet)
                call(progress_callback=call.emitProgressChanged)
                self.results = call.get_result(processEvents=True)
            except obiArrayExpress.GeneAtlasError, ex:
                self.error(0, str(ex))
            finally:
                self.controlArea.setEnabled(True)
                self.update_info_box(query_running=False)
                self.progressBarFinished()
            
            self.query_genes = genes
        
    def update_gene_matcher(self):
        taxid = dict((v,k) for k,v in TAXID_TO_ORG.items()).get(self.selected_organism)
        self.warning(0)
        if taxid != self.loaded_matcher_taxid:
            if taxid:
                self.ensembl_info = obiGene.EnsembleGeneInfo(taxid)
                self.gene_matcher = self.construct_matcher(taxid)
            else:
                self.warning(0, "Cannot match gene names!")
                self.ensembl_info = None
                self.gene_matcher = obiGene.GMDirect()
            self.loded_matcher_taxid = taxid
        
    def update_ef_values_box(self):
        """ Update the "Values" box.
        """
        ef_values = set() 
        results = self.results[self.selected_ef]
        efv_count = defaultdict(int)
        if self.results:
            for gene, val in  results.iteritems():
                for ef_val, (up, down) in val.iteritems():
                   ef_values.add(ef_val)
                   efv_count[ef_val] += up + down
        
        # Sort by up and down count sum of all genes
        self.ef_values = sorted(ef_values, key=efv_count.get, reverse=True)
        self.values_cb.clear()
        self.values_cb.addItems(self.ef_values)
        if self.ef_values:
            self.selected_ef_value = min(len(self.ef_values) - 1, self.selected_ef_value)
        
    def update_report_view(self):
        results = self.results[self.selected_ef]
        
        model = QStandardItemModel()
        
                
        def standard_item(val):
            item = StandardPyItem()
            item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            item.setData(QVariant(val), Qt.DisplayRole)
            return item
        
        if results and self.ef_values:
            ef_value = self.ef_values[self.selected_ef_value]
            for gene, efvs in results.iteritems():
                if ef_value in efvs:
                    up, down = efvs[ef_value]
                    model.appendRow([standard_item(gene),
                                     standard_item(up),
                                     standard_item(down)])
        
        model.setHorizontalHeaderLabels(self.report_header)
        self.report_view.setModel(model)
        self.report_view.resizeColumnToContents(0)
        self.report_view.resizeColumnToContents(1)
        self.report_view.resizeColumnToContents(2)
        
    def update_info_box(self, query_running=False):
        text = ""
        if self.data is not None:
            genes = self.input_genes()
#            unique, matched = self.match_genes(genes, full_res=True)
            # TODO: How many genes are matched (should update after the query is run).
            if self.results and not query_running:
#                matched = min(max([len(r) for r in self.results]), len(genes))
                matched = len(genes) - len(self.unknown_genes)
                with_results = len(reduce(set.union, [d.keys() for d in self.results], set()))
                text = "{0} gene names on input.\n{1} genes matched.\n{2} genes annotated.".format(len(genes), matched, with_results)
            else:
                text = "{0} gene names on input.\n\n".format(len(genes))
        else:
            text = "No data on input.\n\n"
        if query_running:
            text += "Querying Gene Atlas. Please wait."
        self.info_label.setText(text)
            
    def selected_report_genes(self):
        """ Return the gene names selected in the report view.
        """
        rows = self.report_view.selectionModel().selectedRows(0)
        genes = []
        for index in rows:
            name = str(index.data(Qt.DisplayRole).toString())
            genes.append(name)
        return genes
            
    def commit(self):
        genes = set(self.selected_report_genes())
        self.gene_matcher.set_targets(genes)
        
        genes = [gene.lower() for gene in genes]
        if self.genes_in_columns:
            attrs = [a for a in self.data.domain.attributes \
                     if a.name.lower() in genes or self.gene_matcher.umatch(a.name)]
            domain = Orange.data.Domain(attrs, self.data.domain.class_var)
            domain.add_metas(self.data.domain.get_metas())
            data = Orange.data.Table(domain, self.data)
        else:
            attr = self.candidate_vars[self.selected_gene_attr]
            examples = []
            for ex in self.data:
                if str(ex[attr]).lower() in genes or self.gene_matcher.umatch(str(ex[attr])):
                    examples.append(ex)
            if examples:
                data = Orange.data.Table(examples)
            else:
                data = None
        self.send("Selected Genes", data)
        
    def handle_assync_error(self, *args):
        pass
        
        
class StandardPyItem(QStandardItem):
    def __lt__(self, other):
        my = self.data(Qt.DisplayRole).toPyObject()
        other = other.data(Qt.DisplayRole).toPyObject()
        return my < other
        
    
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWGeneAtlasTissueExpression()

    data = Orange.data.Table(
        Orange.data.Domain([Orange.feature.String("Gene")], None),
        [["RUNX1"]]
    )

    w.show()
    w.set_data(data)
    app.exec_()
    w.saveSettings()
        
