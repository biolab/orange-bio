"""<name>Gene Atlas Tissue Expression</name>
<description></description>
"""
import os, sys
import obiArrayExpress
import Orange

import OWGUI
from OWWidget import *

from collections import defaultdict
from Orange.misc import lru_cache

import orngDataCaching
import obiGene

# Mapping for common taxids from obiTaxonomy
TAXID_TO_ORG = {"": "Anopheles gambiae",
                "3702": "Arabidopsis thaliana",
                "9913": "Bos taurus",
                "6239": "Caenorhabditis elegans",
                "7955": "Danio rerio",
                "7227": "Drosophila melanogaster",
                "": "Epstein barr virus",
                "": "Gallus gallus",
                "9606": "Homo sapiens",
                "": "Human cytomegalovirus",
                "": "Kaposi sarcoma-associated herpesvirus",
                "10090": "Mus musculus",
                "10116": "Rattus norvegicus",
                "4932": "Saccharomyces cerevisiae",
                "4896": "Schizosaccharomyces pombe",
                "8355": "Xenopus laevis"
     }

def construct_matcher(taxid):
    return obiGene.matcher([obiGene.GMEnsembl(taxid), obiGene.GMNCBI(taxid)])

class OWGeneAtlasTissueExpression(OWWidget):
    contextHandlers = {"": DomainContextHandler("", ["selected_organism",
                                                     "selected_gene_attr",
                                                     "genes_in_columns",
                                                     "selected_ef",
                                                     "selected_ef_value"])}
    settingsList = []
    
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
        
        self.gene_matcher = obiGene.GMDirect()
        self.query_genes = []
        
        # Cached get_atlas_summary 
        @lru_cache(maxsize=20)
        def get_atlas_summary(genes, organism):
            return obiArrayExpress.get_atlas_summary(list(genes), organism)
        
        self.get_atlas_summary = get_atlas_summary
        
        #Cached construct_matcher
        self.construct_matcher = lru_cache(maxsize=3)(construct_matcher)
        
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
        vars = [var for var in vars if isinstance(var, (Orange.data.variable.String))]
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
        
    def input_genes(self):
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
    
    def run_query(self):
        """ Query the Gene Atlas.
        """
        if self.data is not None:
            genes = self.input_genes()
            # Non threaded
            # self.results = self.get_atlas_summary(tuple(genes), self.selected_organism)
            # Threaded
            self.error(0)
            self.update_info_box(query_running=True)
            self.controlArea.setEnabled(False)
            try:
                call = self.asyncCall(self.get_atlas_summary, (tuple(genes),
                                                               self.selected_organism),
                                      name="Query Gene Expression Atlas",
                                      onError=self.handle_assync_error)
                
                call()
                self.results = call.get_result(processEvents=True)
            except obiArrayExpress.GeneAtlasError, ex:
                self.error(0, str(ex))
            
            finally:
                self.controlArea.setEnabled(True)
                self.update_info_box(query_running=False)
            
            self.query_genes = genes
        
    def update_gene_matcher(self):
        taxid = dict((v,k) for v,k in TAXID_TO_ORG.items()).get(self.selected_organism)
        if taxid:
            self.gene_matcher = self.construct_matcher(taxid)
        else:
            self.gene_matcher = obiGene.GMDirect()
            
        
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
            text = "{0} gene names on input.\n".format(len(genes))
        else:
            text = "No data on input.\n"
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
    data = Orange.data.Table("RUNX1.tab")
    w.show()
    w.set_data(data)
    app.exec_()
    w.saveSettings()
        