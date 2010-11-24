"""<name>GEO DataSets</name>
<description>Access to Gene Expression Omnibus data sets.</description>
<priority>20</priority>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact>
<icon>icons/GEODataSets.png</icon>
"""

from __future__ import with_statement

import sys, os, glob
from OWWidget import *
import OWGUI, OWGUIEx
import obiGEO
import orngServerFiles

from orngDataCaching import data_hints

from collections import defaultdict
from functools import partial 

LOCAL_GDS_COLOR = Qt.darkGreen

class TreeModel(QAbstractItemModel):
    def __init__(self, data, header, parent):
        QAbstractItemModel.__init__(self, parent)
        self._data = [[QVariant(s) for s in row] for row in data]
        self._dataDict = {}
        self._header = {Qt.Horizontal: dict([(i, {Qt.DisplayRole: h}) for i, h in enumerate(header)])}
        self._roleData = {Qt.DisplayRole:self._data}
        dataStore = partial(defaultdict, partial(defaultdict, partial(defaultdict, QVariant)))
        self._roleData = dataStore(self._roleData)
        self._header = dataStore(self._header)
    
    def setColumnLinks(self, column, links):
        font =QFont()
        font.setUnderline(True)
        font = QVariant(font)
        for i, link in enumerate(links):
            self._roleData[LinkRole][i][column] = QVariant(link)
            self._roleData[Qt.FontRole][i][column] = font
            self._roleData[Qt.ForegroundRole][i][column] = QVariant(QColor(Qt.blue))
    
    def setRoleData(self, role, row, col, data):
        self._roleData[role][row][col] = data
        
    def setData(self, index, value, role=Qt.EditRole):
        self._roleData[role][index.row()][index.column()] = value
        self.emit(SIGNAL("dataChanged(QModelIndex, QModelIndex)"), index, index)
        
    def data(self, index, role):
        row, col = index.row(), index.column()
        return self._roleData[role][row][col]
        
    def index(self, row, col, parent=QModelIndex()):
        return self.createIndex(row, col, 0)
    
    def parent(self, index):
        return QModelIndex()
    
    def rowCount(self, index=QModelIndex()):
        if index.isValid():
            return 0
        else:
            return len(self._data)
        
    def columnCount(self, index):
        return len(self._header[Qt.Horizontal])

    def headerData(self, section, orientation, role):
        try:
            return QVariant(self._header[orientation][section][role])
        except KeyError, er:
#            print >> sys.stderr, er
            return QVariant()
        
    def setHeaderData(self, section, orientation, value, role=Qt.EditRole):
        self._header[orientation][section][role] = value
    
from OWGUI import LinkStyledItemDelegate, LinkRole

class OWGEODatasets(OWWidget):
    settingsList = ["outputRows", "mergeSpots", "gdsSelectionStates", "splitterSettings", "currentGds", "autoCommit"]

    def __init__(self, parent=None ,signalManager=None, name=" GEO Data sets"):
        OWWidget.__init__(self, parent ,signalManager, name)

        self.outputs = [("Expression Data", ExampleTable)]

        ## Settings
#        self.selectedSubsets = []
#        self.sampleSubsets = []
        self.selectedAnnotation = 0
        self.includeIf = False
        self.minSamples = 3
        self.autoCommit = False
        self.outputRows = 0
        self.mergeSpots = True
        self.filterString = ""
        self.currentGds = None
        self.selectionChanged = False
        self.autoCommit = False
        self.gdsSelectionStates = {}
        self.splitterSettings = ['\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01\xea\x00\x00\x00\xd7\x01\x00\x00\x00\x07\x01\x00\x00\x00\x02',
                                 '\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01\xb5\x00\x00\x02\x10\x01\x00\x00\x00\x07\x01\x00\x00\x00\x01']

        self.loadSettings()

        ## GUI
        self.infoBox = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Info", addSpace=True), "\n\n")
#        box = OWGUI.widgetBox(self.controlArea, "Sample Subset")
#        OWGUI.listBox(box, self, "selectedSubsets", "sampleSubsets", selectionMode=QListWidget.ExtendedSelection)
#        box = OWGUI.widgetBox(self.controlArea, "Sample Annotations (Types)")
#        self.annotationCombo = OWGUI.comboBox(box, self, "selectedAnnotation", items=["Include all"])
        
##        OWGUI.button(box, self, "Clear selection", callback=self.clearSubsetSelection)
##        c = OWGUI.checkBox(box, self, "includeIf", "Include if at least", callback=self.commitIf)
##        OWGUI.spin(OWGUI.indentedBox(box), self, "minSamples", 2, 100, posttext="samples", callback=self.commitIf)

        box = OWGUI.widgetBox(self.controlArea, "Output", addSpace=True)
        OWGUI.radioButtonsInBox(box, self, "outputRows", ["Genes or spots", "Samples"], "Rows", callback=self.commitIf)
        OWGUI.checkBox(box, self, "mergeSpots", "Merge spots of same gene", callback=self.commitIf)

        box = OWGUI.widgetBox(self.controlArea, "Output", addSpace=True)
        self.commitButton = OWGUI.button(box, self, "Commit", callback=self.commit)
#        self.commitButton.setDisabled(True)
        cb = OWGUI.checkBox(box, self, "autoCommit", "Commit on any change")
        OWGUI.setStopper(self, self.commitButton, cb, "selectionChanged", self.commit)
##        OWGUI.checkBox(box, self, "autoCommit", "Commit automatically")
        OWGUI.rubber(self.controlArea)

        self.filterLineEdit = OWGUIEx.lineEditHint(self.mainArea, self, "filterString", "Filter", caseSensitive=False, matchAnywhere=True, listUpdateCallback=self.filter, callbackOnType=True, callback=self.filter, delimiters=" ")
        splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(splitter)
        self.treeWidget = QTreeView(splitter)
#        splitter.addWidget(self.treeWidget)
        
        self.treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSortingEnabled(True)
        self.treeWidget.setAlternatingRowColors(True)
        self.treeWidget.setItemDelegate(LinkStyledItemDelegate(self.treeWidget))
        self.treeWidget.setItemDelegateForColumn(0, OWGUI.IndicatorItemDelegate(self.treeWidget, role=Qt.DisplayRole))
        
#        self.mainArea.layout().addWidget(self.treeWidget)
        self.connect(self.treeWidget, SIGNAL("itemSelectionChanged ()"), self.updateSelection)
        self.treeWidget.viewport().setMouseTracking(True)
##        self.connect(self.treeWidget, SIGNAL("currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*))"), self.updateSelection_)
#        self.connect(self.treeWidget.model(), SIGNAL("layoutChanged()"), self.filter)
        splitterH = QSplitter(Qt.Horizontal, splitter) 
#        splitter.addWidget(splitterH)
        box = OWGUI.widgetBox(splitterH, "Description")
        self.infoGDS = OWGUI.widgetLabel(box, "")
        self.infoGDS.setWordWrap(True)
        OWGUI.rubber(box)
        
        box = OWGUI.widgetBox(splitterH, "Sample Annotations")
        self.annotationsTree = QTreeWidget(box)
        self.annotationsTree.setHeaderLabels(["Type (Sample annotations)", "Sample count"])
        self.annotationsTree.setRootIsDecorated(True)
        box.layout().addWidget(self.annotationsTree)
        self.connect(self.annotationsTree, SIGNAL("itemChanged(QTreeWidgetItem * , int)"), self.annotationSelectionChanged)
        self._annotationsUpdating = False
        self.splitters = splitter, splitterH
        self.connect(splitter, SIGNAL("splitterMoved(int, int)"), self.splitterMoved)
        self.connect(splitterH, SIGNAL("splitterMoved(int, int)"), self.splitterMoved)
        
        for sp, setting in zip(self.splitters, self.splitterSettings):
            sp.restoreState(setting)
            
        self.searchKeys = ["dataset_id", "title", "platform_organism", "description"]
        self.cells = []
#        self.currentGds = None
        QTimer.singleShot(50, self.updateTable)
        self.resize(1000, 600)

    def updateInfo(self):
        gds_info = obiGEO.GDSInfo()
        text = "%i datasets\n%i datasets cached\n" % (len(gds_info), len(glob.glob(orngServerFiles.localpath("GEO") + "/GDS*")))
        filtered = len([row for row in range(len(self.cells)) if self.rowFiltered(row)])
        if len(self.cells) != filtered:
            text += ("%i after filtering") % filtered
        self.infoBox.setText(text)
        
    def updateTable(self):
        self.treeItems = []
        self.progressBarInit()
        with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
            info = obiGEO.GDSInfo()
        milestones = set(range(0, len(info), max(len(info)/100, 1)))
        self.cells = cells = []
        gdsLinks = []
        pmLinks = []
        localGDS = []
        self.gds = []
        for i, (name, gds) in enumerate(info.items()):
            local = os.path.exists(orngServerFiles.localpath(obiGEO.DOMAIN, gds["dataset_id"] + ".soft.gz"))
            cells.append([" " if local else "", gds["dataset_id"], gds["title"], gds["platform_organism"], len(gds["samples"]), gds["feature_count"],
                          gds["gene_count"], len(gds["subsets"]), gds.get("pubmed_id", "")])

            gdsLinks.append("http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=%s" % gds["dataset_id"])
            pmLinks.append("http://www.ncbi.nlm.nih.gov/pubmed/%s" % gds.get("pubmed_id") if gds.get("pubmed_id") else QVariant())

            if local:
                localGDS.append(i)
            self.gds.append(gds)
            
            if i in milestones:
                self.progressBarSet(100.0*i/len(info))

        model = TreeModel(cells, ["", "ID", "Title", "Organism", "Samples", "Features", "Genes", "Subsets", "PubMedID"], self.treeWidget)
        model.setColumnLinks(1, gdsLinks)
        model.setColumnLinks(8, pmLinks)
#        for i in localGDS:
#            model._roleData[Qt.ForegroundRole][i].update(zip(range(2, 8), [QVariant(QColor(LOCAL_GDS_COLOR))] * 6))
#            mode.._roleData[OWGUI.IndicatorItemDelegate.]
        proxyModel = QSortFilterProxyModel(self.treeWidget)
        proxyModel.setSourceModel(model)
        self.treeWidget.setModel(proxyModel)
        self.connect(self.treeWidget.selectionModel(), SIGNAL("selectionChanged(QItemSelection , QItemSelection )"), self.updateSelection)
        filterItems = " ".join([self.gds[i][key] for i in range(len(self.gds)) for key in self.searchKeys])
        filterItems = reduce(lambda s, d: s.replace(d, " "), [",", ".", ":", "!", "?", "(", ")"], filterItems.lower())
        filterItems = sorted(set(filterItems.split(" ")))
        self.filterLineEdit.setItems(filterItems)
        
        for i in range(8):
            self.treeWidget.resizeColumnToContents(i)
        self.treeWidget.setColumnWidth(1, min(self.treeWidget.columnWidth(1), 300))
        self.treeWidget.setColumnWidth(2, min(self.treeWidget.columnWidth(2), 200))
        self.progressBarFinished()

        if self.currentGds:
            gdss = [(i, model.data(model.index(i,1), Qt.DisplayRole)) for i in range(model.rowCount())]
            current = [i for i, variant in gdss if variant.isValid() and str(variant.toString()) == self.currentGds["dataset_id"]]
            if current:
                mapFromSource = self.treeWidget.model().mapFromSource
                self.treeWidget.selectionModel().select(mapFromSource(model.index(current[0], 0)), QItemSelectionModel.Select | QItemSelectionModel.Rows)
            
        self.updateInfo()

    def updateSelection(self, *args):
        current = self.treeWidget.selectedIndexes()
        mapToSource = self.treeWidget.model().mapToSource
        current = [mapToSource(index).row() for index in current]
        if current:
            self.currentGds = self.gds[current[0]]
            self.setAnnotations(self.currentGds)
            self.infoGDS.setText(self.currentGds.get("description", ""))
        else:
            self.currentGds = None
#        self.commitButton.setDisabled(not bool(self.currentGds))
        self.commitIf()
        
    
    def setAnnotations(self, gds):
#        self.sampleSubsets = ["%s (%d)" % (s["description"], len(s["sample_id"])) for s in gds["subsets"]]
#        self.annotationCombo.clear()
#        self.annotationCombo.addItems(["Include all"] + list(set([sampleinfo["type"] for sampleinfo in gds["subsets"]])))
        
        self._annotationsUpdating = True
        self.annotationsTree.clear()
        annotations = reduce(lambda d, info: d[info["type"]].add(info["description"]) or d, gds["subsets"], defaultdict(set))
        subsetscount = dict([(s["description"], str(len(s["sample_id"]))) for s in gds["subsets"]])
        for type, subsets in annotations.items():
            key = (gds["dataset_id"], type)
            subsetItem = QTreeWidgetItem(self.annotationsTree, [type])
            subsetItem.setFlags(subsetItem.flags() | Qt.ItemIsUserCheckable | Qt.ItemIsTristate)
            subsetItem.setCheckState(0, self.gdsSelectionStates.get(key, Qt.Checked))
            subsetItem.key = key
            for subset in subsets:
                key = (gds["dataset_id"], type, subset)
                item = QTreeWidgetItem(subsetItem, [subset, subsetscount.get(subset, "")])
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(0, self.gdsSelectionStates.get(key, Qt.Checked))
                item.key = key
        self._annotationsUpdating = False
        self.annotationsTree.expandAll()
        for i in range(self.annotationsTree.columnCount()):
            self.annotationsTree.resizeColumnToContents(i)
                
    def annotationSelectionChanged(self, item, column):
        if self._annotationsUpdating:
            return 
        for i in range(self.annotationsTree.topLevelItemCount()):
            item = self.annotationsTree.topLevelItem(i)
            self.gdsSelectionStates[item.key] = item.checkState(0)
            for j in range(item.childCount()):
                child = item.child(j)
                self.gdsSelectionStates[child.key] = child.checkState(0)
        
    def rowFiltered(self, row):
        filterStrings = self.filterString.lower().split()
        try:
            string = " ".join([self.gds[row].get(key, "").lower() for key in self.searchKeys])
            return not all([s in string for s in filterStrings])
        except UnicodeDecodeError:
            string = " ".join([unicode(self.gds[row].get(key, "").lower(), errors="ignore") for key in self.searchKeys])
            return not all([s in string for s in filterStrings])
    
    def filter(self):
        filterStrings = self.filterString.lower().split()
        mapFromSource = self.treeWidget.model().mapFromSource
        index = self.treeWidget.model().sourceModel().index
        
        for i, row in enumerate(self.cells):
            self.treeWidget.setRowHidden(mapFromSource(index(i, 0)).row(), QModelIndex(), self.rowFiltered(i))
        self.updateInfo()

    def commitIf(self):
        if self.autoCommit:
            self.commit()
        else:
            self.selectionChanged = True
            
    def commit(self):
        if self.currentGds:
#            classes = [s["description"] for s in self.currentGds["subsets"]]
#            classes = [classes[i] for i in self.selectedSubsets] or None
#            sample_type = self.annotationCombo.currentText()
            sample_type = None
            self.progressBarInit()
            self.progressBarSet(10)
            
            def getdata(gds_id, **kwargs):
                gds = obiGEO.GDS(gds_id)
                data = gds.getdata(**kwargs)
                return data
            from OWConcurrent import createTask
            self.setEnabled(False)
            qApp.processEvents()
            self.get_data_async = createTask(getdata, (self.currentGds["dataset_id"],), dict(report_genes=self.mergeSpots,
                                                                   transpose=self.outputRows,
                                                                   sample_type=sample_type if sample_type!="Include all" else None),
                                             onResult=self.onData, onFinished=lambda: self.setEnabled(True),
                                             threadPool=QThreadPool.globalInstance()
                                             )
            
#            data = gds.getdata(report_genes=self.mergeSpots, transpose=self.outputRows, sample_type=sample_type if sample_type!="Include all" else None)
    def onData(self, data):
        self.progressBarSet(50)
#            samples = set([self.annotationsTree.topLevelItem(i).child(j).key[2] for i in range(self.annotationsTree.topLevelItemCount()) for j in self.annotationsTree.topLevelItem(i).])
        class itemiter(QTreeWidgetItemIterator):
            def next(self):
                self.__iadd__(1)
                if self.value():
                    return self.value()
                else:
                    raise StopIteration
            def __iter__(self):
                return self
            
        samples = set([str(item.text(0)) for item in itemiter(self.annotationsTree) if self.gdsSelectionStates.get(item.key, True)])
        
        if self.outputRows:
            select = [1 if samples.intersection(str(ex.getclass()).split("|")) else 0 for ex in data]
            data = data.select(select)
            data.domain.classVar.values = ["|".join([cl for cl in val.split("|") if cl in samples]) for val in data.domain.classVar.values]
        else:
            domain = orange.Domain([attr for attr in data.domain.attributes if samples.intersection(attr.attributes.values())], data.domain.classVar)
            domain.addmetas(data.domain.getmetas())
            for attr in domain.attributes:
                attr.attributes = dict([(key, value) for key, value in attr.attributes.items() if value in samples])
            data = orange.ExampleTable(domain, data)

        data_hints.set_hint(data, "taxid", self.currentGds.get("taxid", ""), 10.0)
        data_hints.set_hint(data, "genesinrows", self.outputRows, 10.0)
        
        self.progressBarFinished()
        self.send("Expression Data", data)

        model = self.treeWidget.model().sourceModel()
        row = self.gds.index(self.currentGds)
#            model._roleData[Qt.ForegroundRole][row].update(zip(range(1, 7), [QVariant(QColor(LOCAL_GDS_COLOR))] * 6))
        model.setData(model.index(row, 0),  QVariant(" "), Qt.DisplayRole) 
#            model.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), model.index(row, 0), model.index(row, 0))
        self.updateInfo()
        self.selectionChanged = False
        
    def splitterMoved(self, *args):
        self.splitterSettings = [str(sp.saveState()) for sp in self.splitters]

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWGEODatasets()
    w.show()
    app.exec_()
    w.saveSettings()
