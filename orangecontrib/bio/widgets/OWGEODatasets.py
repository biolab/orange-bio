"""<name>GEO Data Sets</name>
<description>Access to Gene Expression Omnibus data sets.</description>
<priority>20</priority>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact>
<icon>icons/GEODataSets.svg</icon>
"""

from __future__ import absolute_import, with_statement

from collections import defaultdict
from functools import partial 
import sys, os, glob

from Orange.orng import orngServerFiles
from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets import OWGUI, OWGUIEx
from Orange.OrangeWidgets.OWWidget import *

from .. import obiGEO

NAME = "GEO Data Sets"
DESCRIPTION = "Access to Gene Expression Omnibus data sets."
ICON = "icons/GEODataSets.svg"
PRIORITY = 20

INPUTS = []
OUTPUTS = [("Expression Data", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWGEODatasets.OWGEODatasets"]


LOCAL_GDS_COLOR = Qt.darkGreen

TextFilterRole = OWGUI.OrangeUserRole.next()

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
        
from Orange.utils import lru_cache

class MySortFilterProxyModel(QSortFilterProxyModel):    
    def __init__(self, parent=None):
        QSortFilterProxyModel.__init__(self, parent)
        self._filter_strings = []
        self._cache = {}
        self._cache_fixed = {}
        self._cache_prefix = {}
        self._row_text = {}
        
        # Create a cached version of _filteredRows
        self._filteredRows = lru_cache(100)(self._filteredRows) 

    def setSourceModel(self, model):
        """ Set the source model for the filter
        """ 
        self._filter_strings = []
        self._cache = {}
        self._cache_fixed = {}
        self._cache_prefix = {}
        self._row_text = {}
        QSortFilterProxyModel.setSourceModel(self, model)
        
    def addFilterFixedString(self, string, invalidate=True):
        """ Add `string` filter to the list of filters. If invalidate is
        True the filter cache will be recomputed.
        """
        self._filter_strings.append(string)
        all_rows = range(self.sourceModel().rowCount())
        row_text = [self.rowFilterText(row) for row in all_rows]
        self._cache[string] = [string in text for text in row_text]
        if invalidate:
            self.updateCached()
            self.invalidateFilter()
        
    def removeFilterFixedString(self, index=-1, invalidate=True):
        """ Remove the `index`-th filter string. If invalidate is True the
        filter cache will be recomputed.
        """
        string = self._filter_strings.pop(index) 
        del self._cache[string] 
        if invalidate:
            self.updateCached()
            self.invalidateFilter()
            
    def setFilterFixedStrings(self, strings):
        """ Set a list of string to be the new filters.
        """
        s_time = time.time()
        to_remove = set(self._filter_strings) - set(strings)
        to_add = set(strings) - set(self._filter_strings)
        for str in to_remove:
            self.removeFilterFixedString(self._filter_strings.index(str), invalidate=False)
        
        for str in to_add:
            self.addFilterFixedString(str, invalidate=False)
        self.updateCached()
        self.invalidateFilter()
            
    def _filteredRows(self, filter_strings):
        """ Return a dictionary mapping row indexes to True False values.
        .. note:: This helper function is wrapped in the __init__ method. 
        """
        all_rows = range(self.sourceModel().rowCount())
        cache = self._cache
        return dict([(row, all([cache[str][row] for str in filter_strings])) for row in all_rows])
    
    def updateCached(self):
        """ Update the combined filter cache.
        """
        self._cache_fixed = self._filteredRows(tuple(sorted(self._filter_strings))) 
        
    def setFilterFixedString(self, string):
        """Should this raise an error? It is not being used.
        """
        QSortFilterProxyModel.setFilterFixedString(self, string)
        
    def rowFilterText(self, row):
        """ Return text for `row` to filter on. 
        """
        f_role = self.filterRole()
        f_column = self.filterKeyColumn()
        s_model = self.sourceModel()
        data = s_model.data(s_model.index(row, f_column), f_role)
        if isinstance(data, QVariant):
            data = unicode(data.toString(), errors="ignore")
        else:
            data = unicode(data, errors="ignore")
        return data
        
    def filterAcceptsRow(self, row, parent): 
        return self._cache_fixed.get(row, True)
    
    def lessThan(self, left, right):
        if left.column() == 1 and right.column(): # TODO: Remove fixed column handling
            left_gds = str(left.data(Qt.DisplayRole).toString())
            right_gds = str(right.data(Qt.DisplayRole).toString())
            left_gds = left_gds.lstrip("GDS")
            right_gds = right_gds.lstrip("GDS")
            try:
                return int(left_gds) < int(right_gds)
            except Exception, ex:
                pass
        return QSortFilterProxyModel.lessThan(self, left, right)
    
from Orange.OrangeWidgets.OWGUI import LinkStyledItemDelegate, LinkRole

def childiter(item):
    """ Iterate over the children of an QTreeWidgetItem instance.
    """
    for i in range(item.childCount()):
        yield item.child(i)
                
class OWGEODatasets(OWWidget):
    settingsList = ["outputRows", "mergeSpots", "gdsSelectionStates", "splitterSettings", "currentGds", "autoCommit"]

    def __init__(self, parent=None ,signalManager=None, name=" GEO Data Sets"):
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

        box = OWGUI.widgetBox(self.controlArea, "Output", addSpace=True)
        OWGUI.radioButtonsInBox(box, self, "outputRows", ["Genes or spots", "Samples"], "Rows", callback=self.commitIf)
        OWGUI.checkBox(box, self, "mergeSpots", "Merge spots of same gene", callback=self.commitIf)

        box = OWGUI.widgetBox(self.controlArea, "Output", addSpace=True)
        self.commitButton = OWGUI.button(box, self, "Commit", callback=self.commit)
        cb = OWGUI.checkBox(box, self, "autoCommit", "Commit on any change")
        OWGUI.setStopper(self, self.commitButton, cb, "selectionChanged", self.commit)
        OWGUI.rubber(self.controlArea)

        self.filterLineEdit = OWGUIEx.lineEditHint(self.mainArea, self, "filterString", "Filter",
                                   caseSensitive=False, matchAnywhere=True, 
                                   #listUpdateCallback=self.filter, callbackOnType=False, 
                                   callback=self.filter,  delimiters=" ")
        
        splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(splitter)
        self.treeWidget = QTreeView(splitter)
        
        self.treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSortingEnabled(True)
        self.treeWidget.setAlternatingRowColors(True)
        self.treeWidget.setUniformRowHeights(True)
        self.treeWidget.setItemDelegate(LinkStyledItemDelegate(self.treeWidget))
        self.treeWidget.setItemDelegateForColumn(0, OWGUI.IndicatorItemDelegate(self.treeWidget, role=Qt.DisplayRole))
        
        self.connect(self.treeWidget, SIGNAL("itemSelectionChanged ()"), self.updateSelection)
        self.treeWidget.viewport().setMouseTracking(True)
        
        splitterH = QSplitter(Qt.Horizontal, splitter) 
        
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
        
        QTimer.singleShot(50, self.updateTable)
        self.resize(1000, 600)

    def updateInfo(self):
        gds_info = self.gds_info #obiGEO.GDSInfo()
        text = "%i datasets\n%i datasets cached\n" % (len(gds_info), len(glob.glob(orngServerFiles.localpath("GEO") + "/GDS*")))
        filtered = self.treeWidget.model().rowCount()
        if len(self.cells) != filtered:
            text += ("%i after filtering") % filtered
        self.infoBox.setText(text)
        
    def updateTable(self):
        self.treeItems = []
        self.progressBarInit()
        with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
            self.gds_info = info = obiGEO.GDSInfo()

        self.cells = cells = []
        gdsLinks = []
        pmLinks = []
        localGDS = []
        full_text_search_data = []
        self.progressBarSet(10)
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

            full_text_search_data.append(unicode(" | ".join([gds.get(key, "").lower() for key in self.searchKeys]), errors="ignore"))

        self.progressBarSet(20)
        model = TreeModel(cells, ["", "ID", "Title", "Organism", "Samples", "Features", "Genes", "Subsets", "PubMedID"], self.treeWidget)
        model.setColumnLinks(1, gdsLinks)
        model.setColumnLinks(8, pmLinks)

        for i, text in enumerate(full_text_search_data):
            model.setData(model.index(i, 0), QVariant(text), TextFilterRole)

        proxyModel = MySortFilterProxyModel(self.treeWidget)
        proxyModel.setSourceModel(model)
        proxyModel.setFilterKeyColumn(0)
        proxyModel.setFilterRole(TextFilterRole)
        proxyModel.setFilterCaseSensitivity(False)
        proxyModel.setFilterFixedString(self.filterString)

        self.treeWidget.setModel(proxyModel)
        self.connect(self.treeWidget.selectionModel(), SIGNAL("selectionChanged(QItemSelection , QItemSelection )"), self.updateSelection)
        filterItems = " ".join([self.gds[i][key] for i in range(len(self.gds)) for key in self.searchKeys])
        filterItems = reduce(lambda s, d: s.replace(d, " "),
                             [",", ".", ":", ";", "!", "?", "(", ")", "{", "}"
                              "[", "]", "_", "-", "+", "\\", "|", "/", "%", "#",
                              "@", "$", "^", "&", "*", "<", ">", "~", "`"],
                             filterItems.lower())
        filterItems = sorted(set(filterItems.split(" ")))
        filterItems = [item for item in filterItems if len(filterItems) > 3]
        self.filterLineEdit.setItems(filterItems)

        self.progressBarSet(40)

        for i in range(8):
            self.treeWidget.resizeColumnToContents(i)

        self.progressBarSet(70)

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
        filter_string = unicode(self.filterLineEdit.text(), errors="ignore")
        proxyModel = self.treeWidget.model()
        if proxyModel:
            strings = filter_string.lower().strip().split()
            proxyModel.setFilterFixedStrings(strings)
            self.updateInfo()

    def selectedSamples(self):
        """ Return the currently selected sample annotations (list of
        sample type, sample value tuples).
        
        .. note:: if some Sample annotation type has no selected values.
                  this method will return all values for it.
        
        """
        samples = []
        unused_types = []
        used_types = []
        for stype in childiter(self.annotationsTree.invisibleRootItem()): 
            selected_values = []
            all_values = []
            for sval in childiter(stype):
                value = (str(stype.text(0)), str(sval.text(0)))
                if self.gdsSelectionStates.get(sval.key, True):
                    selected_values.append(value)
                all_values.append(value)
            if selected_values:
                samples.extend(selected_values)
                used_types.append(str(stype.text(0)))
            else:
                # If no sample of sample type is selected we don't filter on it.
                samples.extend(all_values)
                unused_types.append(str(stype.text(0)))
        
        return samples, used_types
    
    def commitIf(self):
        if self.autoCommit:
            self.commit()
        else:
            self.selectionChanged = True
    
    def commit(self):
        if self.currentGds:
            self.error(0) 
            sample_type = None
            self.progressBarInit()
            self.progressBarSet(10)
            
            def getdata(gds_id, **kwargs):
                gds = obiGEO.GDS(gds_id)
                data = gds.getdata(**kwargs)
                return data
            
            _, groups = self.selectedSamples()
            if len(groups) == 1 and self.outputRows:
                sample_type = groups[0]

            self.setEnabled(False)
            call = self.asyncCall(getdata, (self.currentGds["dataset_id"],), dict(report_genes=self.mergeSpots,
                                           transpose=self.outputRows,
                                           sample_type=sample_type),
                                  onResult=self.onData,
                                  onFinished=lambda: self.setEnabled(True),
                                  onError=self.onAsyncError,
                                  threadPool=QThreadPool.globalInstance()
                                 )
            call.__call__() #invoke

    def onAsyncError(self, (exctype, value, tb)):
        import ftplib
        if issubclass(exctype, ftplib.error_temp):
            self.error(0, "Can not download dataset from NCBI ftp server! Try again later.")
        elif issubclass(exctype, ftplib.all_errors):
            self.error(0, "Error while connecting to the NCBI ftp server! %s" % str(value))
        else:
            sys.excepthook(exctype, value, tb)
            
        self.progressBarFinished()

    def onData(self, data):
        self.progressBarSet(50)
        
        samples,_ = self.selectedSamples()
        
        self.warning(0)
        message = None
        if self.outputRows:
            def samplesinst(ex):
                out = []
                for i,a in data.domain.get_metas().items():
                    out.append((a.name, ex[i].value))
                if data.domain.class_var.name != 'class':
                    out.append((data.domain.class_var.name, ex[-1].value))
                return out
            samples = set(samples)

            select = [1 if samples.issuperset(samplesinst(ex)) else 0 for ex in data]
            data = data.select(select)
            if len(data) == 0:
                message = "No samples with selected sample annotations."
        else:
            samples = set(samples)
            domain = orange.Domain([attr for attr in data.domain.attributes if samples.issuperset(attr.attributes.items())], data.domain.classVar)
            domain.addmetas(data.domain.getmetas())
            if len(domain.attributes) == 0:
                message = "No samples with selected sample annotations."
            stypes = set(s[0] for s in samples)
            for attr in domain.attributes:
                attr.attributes = dict([(key, value) for key, value in attr.attributes.items() if key in stypes])
            data = orange.ExampleTable(domain, data)
        
        if message is not None:
            self.warning(0, message)
            
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
