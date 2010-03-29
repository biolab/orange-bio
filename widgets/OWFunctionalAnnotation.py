"""<name>Gene Functional Annotation</name>
"""

from OWWidget import *

import obiGeneSets
import obiGene
import obiTaxonomy
import obiProb
import OWGUI

import math

from orngDataCaching import data_hints

from collections import defaultdict

class OWFunctionalAnnotation(OWWidget):
    settingsList = ["spiciesIndex", "genesinrows", "geneattr", "categoriesCheckState"]
    contextHandlers = {"":DomainContextHandler("", ["spiciesIndex", "genesinrows", "geneattr", "categoriesCheckState"])}
    
    def __init__(self, parent=None, signalManager=None, name="Gene Functional Annotation", **kwargs):
        OWWidget.__init__(self, parent, signalManager, name, **kwargs)
        self.inputs = [("Example Table", ExampleTable, self.setData, Default), ("Reference", ExampleTable, self.setReference)]
        self.outputs = [("Selected Examples", ExampleTable)]
        
        self.spiciesIndex = 0
        self.genesinrows = False
        self.geneattr = 0
        self.geneMatcherSettings = [False, False, True, False]
        self.useReferenceData = False
        self.useMinCountFilter = True
        self.useMaxPValFilter = True
        self.minClusterCount = 0
        self.maxPValue = 0.01
        
        self.categoriesCheckState = {}
        
        self.loadSettings()
        
        self.signalManager.setFreeze(1)
        QTimer.singleShot(50, self.updateHierarchy)
            
        self.infoBox = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Info"), "Info")
        self.infoBox.setText("No data on input")
        
        self.spiciesComboBox = OWGUI.comboBox(self.controlArea, self, "spiciesIndex", "Spicies", callback=lambda :self.data and self.updateAnnotations())
        
        box = OWGUI.widgetBox(self.controlArea, "Gene names")
        self.geneAttrComboBox = OWGUI.comboBox(box, self, "geneattr", "Gene attribute", sendSelectedValue=0, callback=self.updateAnnotations)
        OWGUI.checkBox(box, self, "genesinrows", "Use attribute names", callback=lambda :self.data and self.updateAnnotations(), disables=[(-1, self.geneAttrComboBox)])
        OWGUI.button(box, self, "Gene matcher settings", callback=self.updateGeneMatcherSettings, tooltip="Open gene matching settings dialog", debuggingEnabled=0)
        
        self.referenceRadioBox = OWGUI.radioButtonsInBox(self.controlArea, self, "useReferenceData", ["Entire genome", "Reference set (input)"],
                                                        tooltips=["Use entire genome for reference", "Use genes from Referece Examples input signal as reference"],
                                                        box="Reference", callback=self.updateAnnotations)
        
        box = OWGUI.widgetBox(self.controlArea, "Annotation Summary")
        self.groupsWidget = QTreeWidget(self)
        self.groupsWidget.setHeaderLabels(["Category"])
        box.layout().addWidget(self.groupsWidget)

        
        hLayout = QHBoxLayout()
        hLayout.setSpacing(5)
        sb, sbcb = OWGUI.spin(QWidget(self.mainArea), self, "minClusterCount", 0, 100, label="Min. Count", tooltip="Minimum gene count", callback=self.filterAnnotationsChartView, callbackOnReturn=True, checked="useMinCountFilter", checkCallback=self.filterAnnotationsChartView)
        dsp, dspcb = OWGUI.doubleSpin(QWidget(self.mainArea), self, "maxPValue", 0.0, 1.0, 0.0001, label="Max. P-Value", tooltip="Maximum P-Value", callback=self.filterAnnotationsChartView, callbackOnReturn=True, checked="useMaxPValFilter", checkCallback=self.filterAnnotationsChartView)
        
        hLayout.addWidget(sb)
        hLayout.addWidget(sbcb)
        hLayout.addWidget(dsp)
        hLayout.addWidget(dspcb)
        
#        self.controlArea.layout().addLayout(hLayout)
        
#        self.filterLineEdit = QLineEdit(self.mainArea)
        import OWGUIEx
        reload(OWGUIEx)
        self.filterLineEdit = OWGUIEx.QLineEditWithActions(self)
        self.filterLineEdit.setPlaceholderText("Filter ...")
        action = QAction(QIcon(os.path.join(orngEnviron.canvasDir, "icons", "delete.png")), "Clear", self)
#        action = QAction("Clear filter text", self)
        self.filterLineEdit.addAction(action, 0, Qt.AlignHCenter)
        self.connect(action, SIGNAL("triggered()"), self.filterLineEdit.clear)
        
        self.filterCompleter = QCompleter(self.filterLineEdit)
        self.filterLineEdit.setCompleter(self.filterCompleter)
        
        hLayout.addWidget(self.filterLineEdit)
        self.mainArea.layout().addLayout(hLayout)
        
        self.connect(self.filterLineEdit, SIGNAL("textChanged(QString)"), self.filterAnnotationsChartView)
        
        self.annotationsChartView = QTreeWidget(self)
        self.annotationsChartView.setHeaderLabels(["Category", "Term", "Count", "Reference count", "P-Value"])
        self.annotationsChartView.setAlternatingRowColors(True)
        self.annotationsChartView.setSortingEnabled(True)
        self.annotationsChartView.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.mainArea.layout().addWidget(self.annotationsChartView)
        self.taxid_list = []
        
        self.connect(self.groupsWidget, SIGNAL("itemClicked(QTreeWidgetItem *, int)"), self.subsetSelectionChanged)
        
        OWGUI.button(self.controlArea, self, "Commit", callback=self.commit)
        
        self.loadedGenematcher = "None"
        self.referenceData = None
        self.data = None
        
        self.treeItems = []
        
        self.resize(1024, 600)
        
    def updateHierarchy(self):
        try:
            all, local = obiGeneSets.list_all(), obiGeneSets.list_local()
            organisms = set(obiTaxonomy.essential_taxids() + [t[1] for t in all])
            self.taxid_list = list(organisms)
            self.spiciesComboBox.clear()
            self.spiciesComboBox.addItems([obiTaxonomy.name(id) for id in organisms])
            self.genesets = all
            self.emit(SIGNAL("hierarchyUpdated()"))
        finally:
            self.signalManager.setFreeze(0)
        
    def setData(self, data=None):
        self.data = data
        self.closeContext()
        self.geneAttrComboBox.clear()
        self.groupsWidget.clear()
        self.annotationsChartView.clear()
        
        if data:
            self.geneAttrs = [attr for attr in data.domain.variables + data.domain.getmetas().values() \
                              if attr.varType != orange.VarTypes.Continuous]
            self.geneAttrComboBox.addItems([attr.name for attr in self.geneAttrs])
            self.geneattr = min(self.geneattr, len(self.geneAttrs) - 1)
             
            taxid = data_hints.get_hint(data, "taxid", "")
            try:
                self.spiciesIndex = self.taxid_list.index(taxid)
            except Exception, ex:
                pass
            self.genesinrows = data_hints.get_hint(data, "genesinrows", self.genesinrows)
            
            self.openContext("", data)
            
            self.setHierarchy(self.getHierarchy(taxid=self.taxid_list[self.spiciesIndex]))
            
            self.loadedGenematcher = "None"
            self.updateAnnotations()
            
    def setReference(self, data=None):
        self.referenceData = data
        self.referenceRadioBox.setEnabled(bool(data))
        
    def getHierarchy(self, taxid):
        def recursive_dict():
            return defaultdict(recursive_dict)
        collection = recursive_dict()
        
        def collect(col, hier):
            if hier:
                collect(col[hier[0]], hier[1:])
                
        for hierarchy, taxid, _ in self.genesets:
            collect(collection[taxid], hierarchy)
        return collection[taxid]
        
    def setHierarchy(self, hierarchy):
        print self.categoriesCheckState
        def fill(col, parent, full=()):
            for key, value in sorted(col.items()):
                full_cat = full + (key,)
                item = QTreeWidgetItem(parent, [key])
                item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                item.setData(0, Qt.CheckStateRole, QVariant(self.categoriesCheckState.get(full_cat, Qt.Checked)))
                item.setExpanded(True)
                fill(value, item, full_cat)
                
        fill(hierarchy, self.groupsWidget)
        
    def selectedCategories(self):
#        def collect_selected(item):
#            if item.checkState(0) == Qt.Checked:
#                cat = str(item.data(0, Qt.DisplayRole).toString())
#                subcats = reduce(set.union, [collect_selected(item.child(i)) for i in range(item.childCount())], set())
#                return [(cat,) + subcat for subcat in sorted(subcats)] or [(cat,)]
#            else:
#                return []
#            
#        items = [self.groupsWidget.topLevelItem(i) for i in range(self.groupsWidget.topLevelItemCount())]
#        cats = reduce(set.union, [collect_selected(item) for item in items], set())
#        
#        return [(cat, taxid) for cat in cats]
        taxid = self.taxid_list[self.spiciesIndex]
        return [(key, taxid) for key, check in self.getHierarchyCheckState().items() if check == Qt.Checked]

    def getHierarchyCheckState(self):
        def collect(item, full=()):
            checked = item.checkState(0)
            name = str(item.data(0, Qt.DisplayRole).toString())
            full_cat = full + (name,)
            result = [(full_cat, checked)]
            for i in range(item.childCount()):
                result.extend(collect(item.child(i), full_cat))
            return result
            
        items = [self.groupsWidget.topLevelItem(i) for i in range(self.groupsWidget.topLevelItemCount())]
        states = reduce(list.__add__, [collect(item) for item in items], [])
        return dict(states)
            
    def subsetSelectionChanged(self, item, column):
        self.filterAnnotationsChartView()
        self.categoriesCheckState = self.getHierarchyCheckState()
        
    def updateGeneMatcherSettings(self):
        from OWGOEnrichmentAnalysis import GeneMatcherDialog
        dialog = GeneMatcherDialog(self, defaults=self.geneMatcherSettings, enabled=[True] * 4, modal=True)
        if dialog.exec_():
            self.geneMatcherSettings = [getattr(dialog, item[0]) for item in dialog.items]
            self.loadedGenematcher = "None"
            if self.data:
                self.updateAnnotations()
                
    def updateGenematcher(self):
        taxid = self.taxid_list[self.spiciesIndex]
        if taxid != self.loadedGenematcher:
            matchers = [obiGene.GMGO, obiGene.GMKEGG, obiGene.GMNCBI, obiGene.GMAffy]
            self.genematcher = obiGene.matcher([gm(taxid) for gm, use in zip(matchers, self.geneMatcherSettings) if use])
            self.genematcher.set_targets(self.referenceGenes())
            self.loadedGenematcher = taxid
            
    def genesFromExampleTable(self, table):
        if self.genesinrows:
            genes = [attr.name for attr in table.domain.attributes]
        else:
            geneattr = self.geneAttrs[self.geneattr]
            genes = [str(ex[geneattr]) for ex in table]
        return genes
    
    def clusterGenes(self):
        return self.genesFromExampleTable(self.data)
    
    def referenceGenes(self):
        if self.referenceData and self.useReferenceData:
            return self.genesFromExampleTable(self.referenceData)
        else:
            taxid = self.taxid_list[self.spiciesIndex]
            return obiGene.NCBIGeneInfo(taxid).keys()
    
    def _cached_name_lookup(self, func, cache):
        def f(name, cache=cache):
            if name not in cache:
                cache[name] = func(name)
            return cache[name]
        return f
    
    def mapGeneNames(self, names, cache=None):
        if cache is not None:
            umatch = self._cached_name_lookup(self.genematcher.umatch, cache)
        else:
            umatch = self.genematcher.umatch
        return [n for n in [umatch(name) for name in names] if n is not None]
    
    def enrichment(self, geneset, cluster, reference, pval=obiProb.Hypergeometric(), cache=None):
        genes = set(self.mapGeneNames(geneset.genes, cache))
        
        cmapped = genes.intersection(cluster)
        rmapped = genes.intersection(reference)
        return (cmapped, rmapped, pval.p_value(len(cmapped), len(reference), len(rmapped), len(cluster))) # TODO: compute all statistics here
    
    def updateAnnotations(self):
        self.annotationsChartView.clear()
        self.progressBarInit()
        self.updateGenematcher()
        categories = self.selectedCategories()
        collections = list(obiGeneSets.collections(*categories))
        clusterGenes, referenceGenes = self.clusterGenes(), self.referenceGenes()
        cache = {}
        countAll = len(clusterGenes)
        infoText = "%i genes on input\n" % countAll
        referenceGenes = set(self.mapGeneNames(referenceGenes, cache))
        self.progressBarSet(1)
        clusterGenes = set(self.mapGeneNames(clusterGenes, cache))
        self.progressBarSet(2)
        infoText += "%i (%.1f) gene names matched" % (len(clusterGenes), 100.0 * len(clusterGenes) / countAll)
        self.infoBox.setText(infoText)
        
        results = []
        from orngMisc import progressBarMilestones
        
        milestones = progressBarMilestones(len(collections), 100)
        for i, geneset in enumerate(collections):
            results.append((geneset, self.enrichment(geneset, clusterGenes, referenceGenes, cache=cache)))
            if i in milestones:
                self.progressBarSet(100.0 * i / len(collections))
#        results = [(geneset, self.enrichment(geneset, clusterGenes, referenceGenes, cache=cache)) for geneset in collections]
        fmt = lambda score, max_decimals=10: "%%.%if" % min(int(abs(math.log(max(score, 1e-10)))) + 2, max_decimals) if score > math.pow(10, -max_decimals) and score < 1 else "%.1f"
        self.annotationsChartView.clear()
        
        self.filterCompleter.setModel(None)
        
        self.treeItems = []
        for i, (geneset, (cmapped, rmapped, p_val)) in enumerate(results):
            if len(cmapped) > 0:
                item = QTreeWidgetItem(self.annotationsChartView, [" ".join(geneset.hierarchy), geneset.name, "", "", fmt(p_val) % p_val])
                item.setData(2, Qt.DisplayRole, QVariant(len(cmapped)))
                item.setData(3, Qt.DisplayRole, QVariant(len(rmapped)))
                item.setData(4, Qt.DisplayRole, QVariant(p_val))
                item.geneset= geneset
                self.treeItems.append(item)
                
        self.filterCompleter.setModel(self.annotationsChartView.model())
        self.filterCompleter.setCompletionColumn(1)
                
        for i in range(self.annotationsChartView.columnCount()):
            self.annotationsChartView.resizeColumnToContents(i)
            
        self.annotationsChartView.setColumnWidth(1, min(self.annotationsChartView.columnWidth(1), 300))
        self.progressBarFinished()
    
    def filterAnnotationsChartView(self, filterString=""):
        categories = set(" ".join(cat) for cat, taxid in self.selectedCategories())
#        print categories
        filterString = str(self.filterLineEdit.text())
        for item in self.treeItems:
            item_cat = str(item.data(0, Qt.EditRole).toString())
            (count, _), (pval, _) = item.data(2, Qt.EditRole).toInt(), item.data(4, Qt.EditRole).toDouble()
            geneset = item.geneset.name
#            print count, pval
            item.setHidden(item_cat not in categories or (self.useMinCountFilter and count < self.minClusterCount) or (self.useMaxPValFilter and pval > self.maxPValue) or filterString not in geneset)
                
#        filter_string = self.filterString.split()
        
    def commit(self):
        selected = self.annotationsChartView.selectedItems()
        genesets = [item.geneset for item in selected]
        cache = {}
        mappedNames = set(self.mapGeneNames(reduce(set.union, [geneset.genes for geneset in genesets], set()), cache))
        if self.genesinrows:
            mapped = [attr for attr in self.data.domain.attributes if self.genematcher.umatch(attr.name) in mappedNames]
            newdomain = orange.Domain(mapped, self.data.domain.classVar)
            newdomain.addmetas(self.data.domain.getmetas())
            data = orange.ExampleTable(newdomain, self.data)
        else:
            geneattr = self.geneAttrs[self.geneattr]
            selected = [1 if self.genematcher.umatch(str(ex[geneattr])) in mappedNames else 0 for ex in self.data]                
            data = self.data.select(selected)
            
#            if self.appendAnnotations:
#                meta = orange.StringVariable("Annotations")
#                data.domain.addmeta(orange.newmetaid(), meta)
#                for ex in data:
#                    geneattr = self.geneAttrs[self.geneattr]
#                    gene = str(ex[geneattr])
#                    annotations = getgene
        
        self.send("Selected Examples", data)
    
if __name__ == "__main__":
    import cProfile
    
    app = QApplication(sys.argv)
    w = OWFunctionalAnnotation()
    w.updateHierarchy()
    data = orange.ExampleTable("../../../doc/datasets/brown-selected")
#    data = orange.ExampleTable("../human")
#    print cProfile.runctx("w.setData(data)", globals(), locals())
    w.setData(data)
    w.show()
    app.exec_()
    w.saveSettings()
    