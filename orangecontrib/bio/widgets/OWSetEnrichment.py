"""<name>Set Enrichment</name>
<icon>icons/GeneSetEnrichment.svg</icon>
"""

from __future__ import absolute_import, with_statement

import math
import operator
from collections import defaultdict

from Orange.orng import orngEnviron, orngServerFiles
from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWGUI import LinkStyledItemDelegate, LinkRole
from Orange.OrangeWidgets.OWGUI import BarItemDelegate
from Orange.OrangeWidgets.OWWidget import *

from Orange.OrangeWidgets.OWConcurrent import ThreadExecutor, Task

from .utils.download import EnsureDownloaded

from .. import obiGene, obiGeneSets, obiProb, obiTaxonomy

NAME = "Set Enrichment"
DESCRIPTION = ""
ICON = "icons/GeneSetEnrichment.svg"
PRIORITY = 5000

INPUTS = [("Data", Orange.data.Table, "setData", Default),
          ("Reference", Orange.data.Table, "setReference")]
OUTPUTS = [("Data subset", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWSetEnrichment.OWSetEnrichment"]


def gsname(geneset):
    return geneset.name if geneset.name else geneset.id

fmtp = lambda score: "%0.5f" % score if score > 10e-4 else "%0.1e" % score
fmtpdet = lambda score: "%0.9f" % score if score > 10e-4 else "%0.5e" % score


def as_unicode(string):
    if isinstance(string, str):
        return string.decode("utf-8", errors="ignore")
    else:
        return string

# A translation table mapping punctuation, ... to spaces
_TR_TABLE = dict((ord(c), ord(" ")) for c in ".,!?()[]{}:;'\"<>")

def word_split(string):
    """
    Split a string into a list of words.
    """
    return as_unicode(string).translate(_TR_TABLE).split()


def _toPyObject(variant):
    val = variant.toPyObject()
    if isinstance(val, type(NotImplemented)): # PyQt 4.4 converts python int, floats ... to C types
        qtype = variant.type()
        if qtype == QVariant.Double:
            val, ok = variant.toDouble()
        elif qtype == QVariant.Int:
            val, ok = variant.toInt()
        elif qtype == QVariant.LongLong:
            val, ok = variant.toLongLong()
        elif qtype == QVariant.String:
            val = variant.toString()
    return val

class MyTreeWidget(QTreeWidget):
    def paintEvent(self, event):
        QTreeWidget.paintEvent(self, event)
        if getattr(self, "_userMessage", None):
            painter = QPainter(self.viewport())
            font = QFont(self.font())
            font.setPointSize(15)
            painter.setFont(font)
            painter.drawText(self.viewport().geometry(), Qt.AlignCenter, self._userMessage)
            painter.end()

class MyTreeWidgetItem(QTreeWidgetItem):
    def __lt__(self, other):
        if not self.treeWidget():
            return id(self) < id(other)
        column = self.treeWidget().sortColumn()
        if column in [4,5]:
            lhs = _toPyObject(self.data(column, 42))
            rhs = _toPyObject(other.data(column, 42))
        else:
            lhs = _toPyObject(self.data(column, Qt.DisplayRole))
            rhs = _toPyObject(other.data(column, Qt.DisplayRole))
        return lhs < rhs

def name_or_none(id):
    """Return organism name for ncbi taxid or None if not found.
    """
    try:
        return obiTaxonomy.name(id)
    except obiTaxonomy.UnknownSpeciesIdentifier:
        return None

class OWSetEnrichment(OWWidget):
    settingsList = ["speciesIndex", "genesinrows", "geneattr",
                    "categoriesCheckState", "useMinCountFilter",
                    "useMaxPValFilter", "useMaxFDRFilter", "minClusterCount",
                    "maxPValue", "maxFDR", "autocommit"]
    contextHandlers = {"":DomainContextHandler("", ["speciesIndex", "genesinrows", "geneattr", "categoriesCheckState"])}

    def refreshHierarchy(self):
        self.setHierarchy(*self.getHierarchy(taxid=self.taxid_list[self.speciesIndex]))

    def __init__(self, parent=None, signalManager=None, name="Set Enrichment", **kwargs):
        OWWidget.__init__(self, parent, signalManager, name, **kwargs)
        self.inputs = [("Data", ExampleTable, self.setData, Default), ("Reference", ExampleTable, self.setReference)]
        self.outputs = [("Data subset", ExampleTable)]

        self.speciesIndex = 0
        self.genesinrows = False
        self.geneattr = 0
        self.geneMatcherSettings = [False, False, True, False]
        self.useReferenceData = False
        self.useMinCountFilter = True
        self.useMaxPValFilter = True
        self.useMaxFDRFilter = True
        self.minClusterCount = 3
        self.maxPValue = 0.01
        self.maxFDR = 0.01
        self.autocommit = False
        self.categoriesCheckState = {}

        self.loadSettings()
        self._changed = False

        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.infoBox = OWGUI.widgetLabel(box, "Info")
        self.infoBox.setText("No data on input")

        self.speciesComboBox = OWGUI.comboBox(self.controlArea, self,
                      "speciesIndex", "Species",
                      callback=lambda: (self.refreshHierarchy(), self.data and self.updateAnnotations()),
                      debuggingEnabled=0)

        box = OWGUI.widgetBox(self.controlArea, "Entity names")
        self.geneAttrComboBox = OWGUI.comboBox(box, self, "geneattr",
                                "Entity feature",
                                sendSelectedValue=0,
                                callback=self.updateAnnotations)

        cb = OWGUI.checkBox(box, self, "genesinrows", "Use feature names",
                            callback=lambda :self.data and self.updateAnnotations(),
                            disables=[(-1, self.geneAttrComboBox)])
        cb.makeConsistent()

        OWGUI.button(box, self, "Gene matcher settings",
                     callback=self.updateGeneMatcherSettings,
                     tooltip="Open gene matching settings dialog",
                     debuggingEnabled=0)

        self.referenceRadioBox = OWGUI.radioButtonsInBox(self.controlArea,
                    self, "useReferenceData", ["All entities", "Reference set (input)"],
                    tooltips=["Use entire genome (for gene set enrichment) or all available entities for reference",
                              "Use entities from Reference Examples input signal as reference"],
                    box="Reference", callback=self.updateAnnotations)

        box = OWGUI.widgetBox(self.controlArea, "Entity Sets")
        self.groupsWidget = QTreeWidget(self)
        self.groupsWidget.setHeaderLabels(["Category"])
        box.layout().addWidget(self.groupsWidget)

        hLayout = QHBoxLayout()
        hLayout.setSpacing(10)
        hWidget = OWGUI.widgetBox(self.mainArea, orientation=hLayout)
        sb, sbcb = OWGUI.spin(hWidget, self, "minClusterCount",
                              0, 100, label="Entities",
                              tooltip="Minimum entity count",
                              callback=self.filterAnnotationsChartView,
                              callbackOnReturn=True,
                              checked="useMinCountFilter",
                              checkCallback=self.filterAnnotationsChartView)

        dsp, dspcb = OWGUI.doubleSpin(hWidget, self,
                        "maxPValue", 0.0, 1.0, 0.0001,
                        label="p-value",
                        tooltip="Maximum p-value",
                        callback=self.filterAnnotationsChartView,
                        callbackOnReturn=True,
                        checked="useMaxPValFilter",
                        checkCallback=self.filterAnnotationsChartView)

        dsfdr, dsfdrcb = OWGUI.doubleSpin(hWidget, self,
                        "maxFDR", 0.0, 1.0, 0.0001,
                        label="FDR",
                        tooltip="Maximum False discovery rate",
                        callback=self.filterAnnotationsChartView,
                        callbackOnReturn=True,
                        checked="useMaxFDRFilter",
                        checkCallback=self.filterAnnotationsChartView)

        from Orange.OrangeWidgets import OWGUIEx
        self.filterLineEdit = OWGUIEx.QLineEditWithActions(self)
        self.filterLineEdit.setPlaceholderText("Filter ...")
        action = QAction(QIcon(os.path.join(orngEnviron.canvasDir,
                        "icons", "delete_gray.png")), "Clear", self)

        self.filterLineEdit.addAction(action, 0, Qt.AlignHCenter)
        self.connect(action, SIGNAL("triggered()"), self.filterLineEdit.clear)

        self.filterCompleter = QCompleter(self.filterLineEdit)
        self.filterCompleter.setCaseSensitivity(Qt.CaseInsensitive)
        self.filterLineEdit.setCompleter(self.filterCompleter)

        hLayout.addWidget(self.filterLineEdit)
        self.mainArea.layout().addWidget(hWidget)

        self.connect(self.filterLineEdit, SIGNAL("textChanged(QString)"),
                     self.filterAnnotationsChartView)

        self.annotationsChartView = MyTreeWidget(self)
        self.annotationsChartView.setHeaderLabels(["Category", "Term",
                            "Count", "Reference count", "p-value", "FDR", "Enrichment"])
        self.annotationsChartView.setAlternatingRowColors(True)
        self.annotationsChartView.setSortingEnabled(True)
        self.annotationsChartView.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.annotationsChartView.setRootIsDecorated(False)
        self.annotationsChartView.viewport().setMouseTracking(True)
        self.annotationsChartView.itemSelectionChanged.connect(self.invalidate)
        self.mainArea.layout().addWidget(self.annotationsChartView)

        contextEventFilter = OWGUI.VisibleHeaderSectionContextEventFilter(self.annotationsChartView)
        self.annotationsChartView.header().installEventFilter(contextEventFilter)

        self.taxid_list = []

        self.connect(self.groupsWidget, SIGNAL("itemClicked(QTreeWidgetItem *, int)"), self.subsetSelectionChanged)

        box = OWGUI.widgetBox(self.controlArea, "Commit")
        cb = OWGUI.checkBox(box, self, "autocommit", "Commit on any change")
        b = OWGUI.button(box, self, "Commit", callback=self.commit,
                         default=True)
        OWGUI.setStopper(self, b, cb, "_changed", callback=self.commit)
        self.loadedGenematcher = "None"
        self.referenceData = None
        self.data = None

        self.treeItems = []

        self.resize(1024, 600)

        self.connect(self, SIGNAL("widgetStateChanged(QString, int, QString)"), self.onStateChange)

        self.updatingAnnotationsFlag = False
        self.currentAnnotatedCategories = []

        self.setBlocking(True)

        task = EnsureDownloaded(
            [("Taxonomy", "ncbi_taxonomy.tar.gz"),
             (obiGeneSets.sfdomain, "index.pck")]
        )

        def a1():
            oldi = self.speciesIndex
            self.updateHierarchy()
            self.speciesIndex = max(min(oldi, len(self.taxid_list)-1),0)

        task.finished.connect(a1)

        self._executor = ThreadExecutor()
        self._executor.submit(task)

    def no_gene_matching(self):
        #only direct gene matching
        return True if len(self.genematcher.matchers) == 1 else False

    def updateHierarchy(self):
        try:
            all, local = obiGeneSets.list_all(), obiGeneSets.list_local()
            organisms = set(obiTaxonomy.essential_taxids() + filter(None, [t[1] for t in all]))

            organism_names = map(name_or_none, organisms)
            organisms = [-1] +  [taxid for taxid, name in zip(organisms, organism_names) \
                         if name is not None]
            self.speciesComboBox.clear()
            self.taxid_list = list(organisms)
            self.speciesComboBox.addItems(["None" if id < 0 else obiTaxonomy.name(id) for id in self.taxid_list])
            self.genesets = all
        finally:
            self.setBlocking(False)

    def setData(self, data=None):
        self.data = data
        self.error(0)
        self.closeContext("")
        self.geneAttrComboBox.clear()
        self.groupsWidget.clear()
        self.annotationsChartView.clear()

        if not getattr(self,"taxid_list", None):
            QTimer.singleShot(100, lambda data=data: self.setData(data))
            return
        if data:
            self.geneAttrs = [attr for attr in data.domain.variables + data.domain.getmetas().values() \
                              if attr.varType != orange.VarTypes.Continuous]

            self.geneAttrComboBox.addItems([attr.name for attr in self.geneAttrs])
            self.geneattr = min(self.geneattr, len(self.geneAttrs) - 1)

            taxid = data_hints.get_hint(data, "taxid", "")
            try:
                self.speciesIndex = self.taxid_list.index(taxid)
            except ValueError, ex:
                pass
            self.genesinrows = data_hints.get_hint(data, "genesinrows", self.genesinrows)

            self.openContext("", data)
        
            self.refreshHierarchy()

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

        for hierarchy, t_id, _ in self.genesets:
            collect(collection[t_id], hierarchy)

        return (taxid, collection[taxid]), (None, collection[None])

    def setHierarchy(self, hierarchy, hierarchy_noorg):
        self.groupsWidgetItems = {}
        def fill(col, parent, full=(), org=""):
            for key, value in sorted(col.items()):
                full_cat = full + (key,)
                item = QTreeWidgetItem(parent, [key])
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                if value:
                    item.setFlags(item.flags() | Qt.ItemIsTristate)

                item.setData(0, Qt.CheckStateRole, QVariant(self.categoriesCheckState.get((full_cat, org), Qt.Checked)))
                item.setExpanded(True)
                item.category = full_cat
                item.organism = org
                self.groupsWidgetItems[full_cat] = item
                fill(value, item, full_cat, org=org)

        self.groupsWidget.clear()
        fill(hierarchy[1], self.groupsWidget, org=hierarchy[0])
        fill(hierarchy_noorg[1], self.groupsWidget, org=hierarchy_noorg[0])

#    def updateCategoryCounts(self):
#        for cat, item in self.groupWidgetItem:
#            item.setData(1, QVariant(), Qt.DisplayRole)

    def selectedCategories(self):
        return [(key, org) for (key, org), check in self.getHierarchyCheckState().items() if check == Qt.Checked]

    def getHierarchyCheckState(self):
        def collect(item, full=()):
            checked = item.checkState(0)
            name = str(item.data(0, Qt.DisplayRole).toString())
            full_cat = full + (name,)
            result = [((full_cat, item.organism), checked)]
            for i in range(item.childCount()):
                result.extend(collect(item.child(i), full_cat))
            return result

        items = [self.groupsWidget.topLevelItem(i) for i in range(self.groupsWidget.topLevelItemCount())]
        states = reduce(list.__add__, [collect(item) for item in items], [])
        return dict(states)

    def subsetSelectionChanged(self, item, column):
        self.categoriesCheckState = self.getHierarchyCheckState()
        categories = self.selectedCategories()

        if self.data is not None:
            if self.no_gene_matching() or not set(categories) <= set(self.currentAnnotatedCategories):
                self.updateAnnotations()
            else:
                self.filterAnnotationsChartView()

    def updateGeneMatcherSettings(self):
        from .OWGOEnrichmentAnalysis import GeneMatcherDialog
        dialog = GeneMatcherDialog(self, defaults=self.geneMatcherSettings, enabled=[True] * 4, modal=True)
        if dialog.exec_():
            self.geneMatcherSettings = [getattr(dialog, item[0]) for item in dialog.items]
            self.loadedGenematcher = "None"
            if self.data:
                self.updateAnnotations()

    def updateGenematcher(self):
        taxid = self.taxid_list[self.speciesIndex]
        if taxid != self.loadedGenematcher:
            self.progressBarInit()
            call = self.asyncCall(obiGene.matcher, name="Gene Matcher", blocking=True, thread=self.thread())
            call.connect(call, SIGNAL("progressChanged(float)"), self.progressBarSet)
            with orngServerFiles.DownloadProgress.setredirect(call.emitProgressChanged):
#            with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
                matchers = [obiGene.GMGO, obiGene.GMKEGG, obiGene.GMNCBI, obiGene.GMAffy]
                if any(self.geneMatcherSettings) and taxid != -1:
                    call.__call__([gm(taxid) for gm, use in zip(matchers, self.geneMatcherSettings) if use])
                    self.genematcher = call.get_result()
#                    self.genematcher = obiGene.matcher([gm(taxid) for gm, use in zip(matchers, self.geneMatcherSettings) if use])
                else:
                    call.__call__([])
                    self.genematcher = call.get_result()
                self.loadedGenematcher = taxid
            self.progressBarFinished()

    def genesFromExampleTable(self, table):
        if self.genesinrows:
            genes = [attr.name for attr in table.domain.attributes]
        else:
            geneattr = self.geneAttrs[self.geneattr]
            genes = [str(ex[geneattr]) for ex in table]
        return genes

    def clusterGenes(self):
        return self.genesFromExampleTable(self.data)

    def reference_from_ncbi(self):
        taxid = self.taxid_list[self.speciesIndex]
        call = self.asyncCall(obiGene.NCBIGeneInfo, (taxid,), name="Load reference genes", blocking=True, thread=self.thread())
        call.connect(call, SIGNAL("progressChanged(float)"), self.progressBarSet)
        with orngServerFiles.DownloadProgress.setredirect(call.emitProgressChanged):
            call.__call__()
            return call.get_result()

    def _cached_name_lookup(self, func, cache):
        def f(name, cache=cache):
            if name not in cache:
                cache[name] = func(name)
            return cache[name]
        return f

    def mapGeneNames(self, names, cache=None, passUnknown=False):
        if cache is not None:
            umatch = self._cached_name_lookup(self.genematcher.umatch, cache)
        else:
            umatch = self.genematcher.umatch
        if passUnknown:
            return [umatch(name) or name for name in names]
#            return [(mapped_name or name, mapped_name is not None) for mapped_name, name in zip(mapped, names)]
        return [n for n in [umatch(name) for name in names] if n is not None]

    def enrichment(self, geneset, cluster, reference, pval=obiProb.Hypergeometric(), cache=None):
        genes = set(self.mapGeneNames(geneset.genes, cache, passUnknown=False))

        cmapped = genes.intersection(cluster)
        rmapped = genes.intersection(reference)
        return (cmapped, rmapped, pval.p_value(len(cmapped), len(reference), len(rmapped), len(cluster)), float(len(cmapped)) / (len(cluster) or 1) / (float(len(rmapped) or 1) / (len(reference) or 1))) # TODO: compute all statistics here

    def updateAnnotations(self):
        if not self.taxid_list:
            return
        self.updatingAnnotationsFlag = True
        self.annotationsChartView.clear()
        self.error([0, 1])
        if not self.genesinrows and len(self.geneAttrs) == 0:
            self.error(0, "Input data contains no attributes with gene names")
            self.updatingAnnotationsFlag = False
            return

        self.updateGenematcher()

        self.progressBarInit()
        self.currentAnnotatedCategories = categories = self.selectedCategories()

        ## Load collections in a worker thread
        call = self.asyncCall(obiGeneSets.collections, categories, name="Loading collections", blocking=True, thread=self.thread())
        call.connect(call, SIGNAL("progressChanged(float)"), self.progressBarSet)
        with orngServerFiles.DownloadProgress.setredirect(call.emitProgressChanged):
            call.__call__()
            collections = list(call.get_result())

        clusterGenes = self.clusterGenes()
        cache = {}

        referenceGenes = self.genesFromExampleTable(self.referenceData) \
            if (self.referenceData and self.useReferenceData) else None

        if self.no_gene_matching():
            if referenceGenes == None:
                referenceGenes = set()
                for g in collections:
                    referenceGenes.update(set(g.genes))
            self.genematcher.set_targets(referenceGenes)
        else:
            if referenceGenes == None:
                referenceGenes = self.reference_from_ncbi()
            self.genematcher.set_targets(referenceGenes)
            referenceGenes = set(self.mapGeneNames(referenceGenes, cache, passUnknown=False))

        countAll = len(set(clusterGenes))
        infoText = "%i unique names on input\n" % countAll
        self.progressBarSet(1)
        clusterGenes = set(self.mapGeneNames(clusterGenes, cache, passUnknown=False))
        
        self.progressBarSet(2)
        if self.no_gene_matching():
            pass
        else:
            self.progressBarSet(2)
            infoText += "%i (%.1f) gene names matched" % (len(clusterGenes), 100.0 * len(clusterGenes) / countAll)
        self.infoBox.setText(infoText)

        results = []
        from Orange.orng.orngMisc import progressBarMilestones

        milestones = progressBarMilestones(len(collections), 100)
        for i, geneset in enumerate(collections):
            results.append((geneset, self.enrichment(geneset, clusterGenes, referenceGenes, cache=cache)))
            if i in milestones:
                self.progressBarSet(100.0 * i / len(collections))

        self.annotationsChartView.clear()

        maxCount = max([len(cm) for _, (cm, _, _, _) in results] + [1])
        maxRefCount = max([len(rc) for _, (_, rc, _, _) in results] + [1])
        countSpaces = int(math.ceil(math.log10(maxCount)))
        refSpaces = int(math.ceil(math.log(maxRefCount)))
        countFmt = "%"+str(countSpaces) + "s  (%.2f%%)"
        refFmt = "%"+str(refSpaces) + "s  (%.2f%%)"

        self.filterCompleter.setModel(None)
        linkFont = QFont(self.annotationsChartView.viewOptions().font)
        linkFont.setUnderline(True)
        self.treeItems = []
        for i, (geneset, (cmapped, rmapped, p_val, enrichment)) in enumerate(results):
            if len(cmapped) > 0:
                item = MyTreeWidgetItem(self.annotationsChartView, [", ".join(geneset.hierarchy), gsname(geneset)])
                item.setData(2, Qt.DisplayRole, QVariant(countFmt % (len(cmapped), 100.0*len(cmapped)/countAll)))
                item.setData(2, Qt.ToolTipRole, QVariant(len(cmapped))) # For filtering
                item.setData(3, Qt.DisplayRole, QVariant(refFmt % (len(rmapped), 100.0*len(rmapped)/len(referenceGenes))))
                item.setData(4, Qt.ToolTipRole, QVariant(fmtpdet(p_val)))
                item.setData(4, Qt.DisplayRole, QVariant(fmtp(p_val)))
                item.setData(4, 42, QVariant(p_val))
                #column 5, FDR, is computed in filterAnnotationsChartView
                item.setData(6, Qt.DisplayRole, QVariant(enrichment))
                item.setData(6, Qt.ToolTipRole, QVariant("%.3f" % enrichment))
                item.geneset= geneset
                self.treeItems.append(item)
                if geneset.link:
                    item.setData(1, LinkRole, QVariant(geneset.link))
                    item.setToolTip(1, geneset.link)
                    item.setFont(1, linkFont)
                    item.setForeground(1, QColor(Qt.blue))

        if not self.treeItems:
            self.warning(0, "No enriched sets found.")
        else:
            self.warning(0)

        allnames = set(as_unicode(gsname(geneset))
                       for geneset, (count, _, _, _) in results if count)

        allnames |= reduce(operator.ior,
                           (set(word_split(name)) for name in allnames),
                           set())

        self.completerModel = QStringListModel(sorted(allnames))
        self.filterCompleter.setModel(self.completerModel)

        
        if results:
            self.annotationsChartView.setItemDelegateForColumn(6, BarItemDelegate(self, scale=(0.0, max(t[1][3] for t in results))))
        self.annotationsChartView.setItemDelegateForColumn(1, LinkStyledItemDelegate(self.annotationsChartView))

        for i in range(self.annotationsChartView.columnCount()):
            self.annotationsChartView.resizeColumnToContents(i)

        self.annotationsChartView.setColumnWidth(1, min(self.annotationsChartView.columnWidth(1), 300))
        self.progressBarFinished()
        self.updatingAnnotationsFlag = False
        QTimer.singleShot(50, self.filterAnnotationsChartView)

    def filterAnnotationsChartView(self, filterString=""):
        if self.updatingAnnotationsFlag:
            return
        categories = set(", ".join(cat) for cat, taxid in self.selectedCategories())

    
        filterString = str(self.filterLineEdit.text()).lower()

        #hide categories
        itemsHiddenCat = []
        for item in self.treeItems:
            item_cat = str(item.data(0, Qt.EditRole).toString())
            geneset = gsname(item.geneset).lower()
            hidden = item_cat not in categories
            itemsHiddenCat.append(hidden)
        
        #compute FDR according the selected categories
        pvals = [ _toPyObject(item.data(4, 42)) for item, hidden in zip(self.treeItems, itemsHiddenCat) if not hidden ]
        fdrs = obiProb.FDR(pvals)

        #update FDR for the selected collections and apply filtering rules
        fdri = 0
        itemsHidden = []
        for item, hidden in zip(self.treeItems, itemsHiddenCat):
            if not hidden:
                fdr = fdrs[fdri]
                fdri += 1

                count, pval = _toPyObject(item.data(2, Qt.ToolTipRole)), _toPyObject(item.data(4, 42))

                hidden = (self.useMinCountFilter and count < self.minClusterCount) or \
                         (self.useMaxPValFilter and pval > self.maxPValue) or \
                         (self.useMaxFDRFilter and fdr > self.maxFDR)
                
                if not hidden:
                    item.setData(5, Qt.ToolTipRole, QVariant(fmtpdet(fdr)))
                    item.setData(5, Qt.DisplayRole, QVariant(fmtp(fdr)))
                    item.setData(5, 42, QVariant(fdr))

            item.setHidden(hidden)
            itemsHidden.append(hidden)
                

        if self.treeItems and all(itemsHidden):
            self.information(0, "All sets were filtered out.")
        else:
            self.information(0)

    def invalidate(self):
        if self.autocommit:
            self.commit()
        else:
            self._changed = True

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
            selected = [1 if self.genematcher.umatch(str(ex[geneattr])) in mappedNames else 0
                        for ex in self.data]
            data = self.data.select(selected)

        self.send("Data subset", data)
        self._changed = False

    def sendReport(self):
        self.reportSettings("Settings", [("Organism", obiTaxonomy.name(self.taxid_list[self.speciesIndex]))])
        self.reportSettings("Filter", [("Min cluster size", self.minClusterCount if self.useMinCountFilter else 0),
                                       ("Max p-value", self.maxPValue if self.useMaxPValFilter else 1.0),
                                       ("Max FDR", self.maxFDR if self.useMaxFDRFilter else 1.0)])

        self.reportSubsection("Annotations")
        self.reportRaw(reportItemView(self.annotationsChartView))

    def onStateChange(self, stateType, id, text):
        if stateType == "Warning" or stateType == "Info":
            self.annotationsChartView._userMessage = text
            self.annotationsChartView.viewport().update()


def reportItemView(view):
    model = view.model()
    return reportItemModel(view, model)

def reportItemModel(view, model, index=QModelIndex()):
    if not index.isValid() or model.hasChildren(index):
        columnCount, rowCount = model.columnCount(index), model.rowCount(index)
        if not index.isValid():
            text = '<table>\n<tr>' + ''.join('<th>%s</th>' % model.headerData(i, Qt.Horizontal, Qt.DisplayRole).toString() for i in range(columnCount)) +'</tr>\n'
        else:
#            variant = model.data(index, Qt.DisplayRole)
#            text = '<table' + (' caption="%s"' % variant.toString() if variant.isValid() else '') + '>\n'
            pass
        text += ''.join('<tr>' + ''.join('<td>' + reportItemModel(view, model, model.index(row, column, index)) + '</td>' for column in range(columnCount)) + '</tr>\n' for row in range(rowCount) if not view.isRowHidden(row, index))
        text += '</table>'
        return text
    else:
        variant = model.data(index, Qt.DisplayRole)
        return str(variant.toString()) if variant.isValid() else ""



if __name__ == "__main__":  
    import cProfile

    app = QApplication(sys.argv)
    w = OWSetEnrichment()
    #data = orange.ExampleTable("yeast-class-RPR.tab")
    #data = orange.ExampleTable("/home/marko/orange-pubchem-data/pug/chems_matrix.tab")
    data = orange.ExampleTable("/home/marko/comp81_master.tab")
    w.loadSettings()
#    data = orange.ExampleTable("../human")
#    print cProfile.runctx("w.setData(data)", globals(), locals())
    w.setData(data)
    #w.setReference(data)
    w.show()
    app.exec_()
    w.saveSettings()


