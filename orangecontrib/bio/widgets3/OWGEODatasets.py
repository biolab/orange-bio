"""<name>GEO Data Sets</name>
<description>Access to Gene Expression Omnibus data sets.</description>
<priority>20</priority>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact>
<icon>icons/GEODataSets.svg</icon>
"""

from __future__ import absolute_import, with_statement

import sys
import os
import glob
import string
import urllib2
from collections import defaultdict
from functools import partial

from Orange.utils import lru_cache
from Orange.utils import serverfiles
from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets import OWGUI, OWGUIEx
from Orange.OrangeWidgets.OWWidget import *

from Orange.OrangeWidgets.OWConcurrent import (
    ThreadExecutor, Task, methodinvoke
)

from .. import geo

NAME = "GEO Data Sets"
DESCRIPTION = "Access to Gene Expression Omnibus data sets."
ICON = "icons/GEODataSets.svg"
PRIORITY = 20

INPUTS = []
OUTPUTS = [("Expression Data", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWGEODatasets.OWGEODatasets"]


TextFilterRole = OWGUI.OrangeUserRole.next()


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
        """Set the source model for the filter.
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
            self.invalidate()

    def setFilterFixedStrings(self, strings):
        """Set a list of string to be the new filters.
        """
        to_remove = set(self._filter_strings) - set(strings)
        to_add = set(strings) - set(self._filter_strings)
        for str in to_remove:
            self.removeFilterFixedString(
                self._filter_strings.index(str),
                invalidate=False)

        for str in to_add:
            self.addFilterFixedString(str, invalidate=False)
        self.updateCached()
        self.invalidate()

    def _filteredRows(self, filter_strings):
        """Return a dictionary mapping row indexes to True False values.

        .. note:: This helper function is wrapped in the __init__ method.

        """
        all_rows = range(self.sourceModel().rowCount())
        cache = self._cache
        return dict([(row, all([cache[str][row] for str in filter_strings]))
                     for row in all_rows])

    def updateCached(self):
        """Update the combined filter cache.
        """
        self._cache_fixed = self._filteredRows(
            tuple(sorted(self._filter_strings)))

    def setFilterFixedString(self, string):
        """Should this raise an error? It is not being used.
        """
        QSortFilterProxyModel.setFilterFixedString(self, string)

    def rowFilterText(self, row):
        """Return text for `row` to filter on.
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
        # TODO: Remove fixed column handling
        if left.column() == 1 and right.column():
            left_gds = str(left.data(Qt.DisplayRole).toString())
            right_gds = str(right.data(Qt.DisplayRole).toString())
            left_gds = left_gds.lstrip("GDS")
            right_gds = right_gds.lstrip("GDS")
            try:
                return int(left_gds) < int(right_gds)
            except ValueError:
                pass
        return QSortFilterProxyModel.lessThan(self, left, right)

from Orange.OrangeWidgets.OWGUI import LinkStyledItemDelegate, LinkRole


def childiter(item):
    """ Iterate over the children of an QTreeWidgetItem instance.
    """
    for i in range(item.childCount()):
        yield item.child(i)


class OWGEODatasets(OWWidget):
    settingsList = ["outputRows", "mergeSpots", "gdsSelectionStates",
                    "splitterSettings", "currentGds", "autoCommit",
                    "datasetNames"]

    def __init__(self, parent=None, signalManager=None, name=" GEO Data Sets"):
        OWWidget.__init__(self, parent, signalManager, name)

        self.outputs = [("Expression Data", ExampleTable)]

        ## Settings
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
        self.splitterSettings = [
            '\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01\xea\x00\x00\x00\xd7\x01\x00\x00\x00\x07\x01\x00\x00\x00\x02',
            '\x00\x00\x00\xff\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x01\xb5\x00\x00\x02\x10\x01\x00\x00\x00\x07\x01\x00\x00\x00\x01'
        ]
        self.datasetNames = {}
        self.loadSettings()

        self.datasetName = ""

        ## GUI
        self.infoBox = OWGUI.widgetLabel(
            OWGUI.widgetBox(self.controlArea, "Info", addSpace=True),
            "Initializing\n\n"
        )

        box = OWGUI.widgetBox(self.controlArea, "Output", addSpace=True)
        OWGUI.radioButtonsInBox(box, self, "outputRows",
                                ["Genes or spots", "Samples"], "Rows",
                                callback=self.commitIf)
        OWGUI.checkBox(box, self, "mergeSpots", "Merge spots of same gene",
                       callback=self.commitIf)

        OWGUI.separator(box)
        self.nameEdit = OWGUI.lineEdit(
            box, self, "datasetName", "Data set name",
            tooltip="Override the default output data set name",
            callback=self.onNameEdited
        )
        self.nameEdit.setPlaceholderText("")

        box = OWGUI.widgetBox(self.controlArea, "Commit", addSpace=True)
        self.commitButton = OWGUI.button(box, self, "Commit",
                                         callback=self.commit)
        cb = OWGUI.checkBox(box, self, "autoCommit", "Commit on any change")
        OWGUI.setStopper(self, self.commitButton, cb, "selectionChanged",
                         self.commit)
        OWGUI.rubber(self.controlArea)

        self.filterLineEdit = OWGUIEx.lineEditHint(
            self.mainArea, self, "filterString", "Filter",
            caseSensitive=False, matchAnywhere=True,
            callback=self.filter,  delimiters=" ")

        splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(splitter)
        self.treeWidget = QTreeView(splitter)

        self.treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSortingEnabled(True)
        self.treeWidget.setAlternatingRowColors(True)
        self.treeWidget.setUniformRowHeights(True)
        self.treeWidget.setEditTriggers(QTreeView.NoEditTriggers)

        linkdelegate = LinkStyledItemDelegate(self.treeWidget)
        self.treeWidget.setItemDelegateForColumn(1, linkdelegate)
        self.treeWidget.setItemDelegateForColumn(8, linkdelegate)
        self.treeWidget.setItemDelegateForColumn(
            0, OWGUI.IndicatorItemDelegate(self.treeWidget,
                                           role=Qt.DisplayRole))

        proxyModel = MySortFilterProxyModel(self.treeWidget)
        self.treeWidget.setModel(proxyModel)
        self.treeWidget.selectionModel().selectionChanged.connect(
            self.updateSelection
        )
        self.treeWidget.viewport().setMouseTracking(True)

        splitterH = QSplitter(Qt.Horizontal, splitter)

        box = OWGUI.widgetBox(splitterH, "Description")
        self.infoGDS = OWGUI.widgetLabel(box, "")
        self.infoGDS.setWordWrap(True)
        OWGUI.rubber(box)

        box = OWGUI.widgetBox(splitterH, "Sample Annotations")
        self.annotationsTree = QTreeWidget(box)
        self.annotationsTree.setHeaderLabels(
            ["Type (Sample annotations)", "Sample count"]
        )
        self.annotationsTree.setRootIsDecorated(True)
        box.layout().addWidget(self.annotationsTree)
        self.annotationsTree.itemChanged.connect(
            self.annotationSelectionChanged
        )
        self._annotationsUpdating = False
        self.splitters = splitter, splitterH

        for sp, setting in zip(self.splitters, self.splitterSettings):
            sp.splitterMoved.connect(self.splitterMoved)
            sp.restoreState(setting)

        self.searchKeys = ["dataset_id", "title", "platform_organism",
                           "description"]

        self.gds = []
        self.gds_info = None

        self.resize(1000, 600)

        self.setBlocking(True)
        self.setEnabled(False)
        self.progressBarInit()

        self._executor = ThreadExecutor()

        func = partial(get_gds_model,
                       methodinvoke(self, "_setProgress", (float,)))
        self._inittask = Task(function=func)
        self._inittask.finished.connect(self._initializemodel)
        self._executor.submit(self._inittask)

        self._datatask = None

    @pyqtSlot(float)
    def _setProgress(self, value):
        self.progressBarValue = value

    def _initializemodel(self):
        assert self.thread() is QThread.currentThread()
        model, self.gds_info, self.gds = self._inittask.result()
        model.setParent(self)

        proxy = self.treeWidget.model()
        proxy.setFilterKeyColumn(0)
        proxy.setFilterRole(TextFilterRole)
        proxy.setFilterCaseSensitivity(False)
        proxy.setFilterFixedString(self.filterString)

        proxy.setSourceModel(model)
        proxy.sort(0, Qt.DescendingOrder)

        self.progressBarFinished()
        self.setBlocking(False)
        self.setEnabled(True)

        filter_items = " ".join(
            gds[key] for gds in self.gds for key in self.searchKeys
        )
        tr_chars = ",.:;!?(){}[]_-+\\|/%#@$^&*<>~`"
        tr_table = string.maketrans(tr_chars, " " * len(tr_chars))
        filter_items = filter_items.translate(tr_table)

        filter_items = sorted(set(filter_items.split(" ")))
        filter_items = [item for item in filter_items if len(item) > 3]
        self.filterLineEdit.setItems(filter_items)

        if self.currentGds:
            gdss = [(i, proxy.data(proxy.index(i, 1), Qt.DisplayRole))
                    for i in range(proxy.rowCount())]
            current = [i for i, variant in gdss
                       if variant.isValid() and
                       str(variant.toString()) == self.currentGds["dataset_id"]]
            if current:
                current_index = proxy.index(current[0], 0)
                self.treeWidget.selectionModel().select(
                    current_index,
                    QItemSelectionModel.Select |
                    QItemSelectionModel.Rows
                )
                self.treeWidget.scrollTo(current_index)

        for i in range(8):
            self.treeWidget.resizeColumnToContents(i)

        self.treeWidget.setColumnWidth(
            1, min(self.treeWidget.columnWidth(1), 300))
        self.treeWidget.setColumnWidth(
            2, min(self.treeWidget.columnWidth(2), 200))

        self.updateInfo()

    def updateInfo(self):
        gds_info = self.gds_info
        text = ("%i datasets\n%i datasets cached\n" %
                (len(gds_info),
                 len(glob.glob(serverfiles.localpath("GEO") + "/GDS*"))))
        filtered = self.treeWidget.model().rowCount()
        if len(self.gds) != filtered:
            text += ("%i after filtering") % filtered
        self.infoBox.setText(text)

    def updateSelection(self, *args):
        current = self.treeWidget.selectedIndexes()
        mapToSource = self.treeWidget.model().mapToSource
        current = [mapToSource(index).row() for index in current]
        if current:
            self.currentGds = self.gds[current[0]]
            self.setAnnotations(self.currentGds)
            self.infoGDS.setText(self.currentGds.get("description", ""))
            self.nameEdit.setPlaceholderText(self.currentGds["title"])
            self.datasetName = \
                self.datasetNames.get(self.currentGds["dataset_id"], "")
        else:
            self.currentGds = None
            self.nameEdit.setPlaceholderText("")
            self.datasetName = ""

        self.commitIf()

    def setAnnotations(self, gds):
        self._annotationsUpdating = True
        self.annotationsTree.clear()

        annotations = defaultdict(set)
        subsetscount = {}
        for desc in gds["subsets"]:
            annotations[desc["type"]].add(desc["description"])
            subsetscount[desc["description"]] = str(len(desc["sample_id"]))

        for type, subsets in annotations.items():
            key = (gds["dataset_id"], type)
            subsetItem = QTreeWidgetItem(self.annotationsTree, [type])
            subsetItem.setFlags(subsetItem.flags() | Qt.ItemIsUserCheckable |
                                Qt.ItemIsTristate)
            subsetItem.setCheckState(
                0, self.gdsSelectionStates.get(key, Qt.Checked)
            )
            subsetItem.key = key
            for subset in subsets:
                key = (gds["dataset_id"], type, subset)
                item = QTreeWidgetItem(
                    subsetItem, [subset, subsetscount.get(subset, "")]
                )
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(
                    0, self.gdsSelectionStates.get(key, Qt.Checked)
                )
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

    def filter(self):
        filter_string = unicode(self.filterLineEdit.text(), errors="ignore")
        proxyModel = self.treeWidget.model()
        if proxyModel:
            strings = filter_string.lower().strip().split()
            proxyModel.setFilterFixedStrings(strings)
            self.updateInfo()

    def selectedSamples(self):
        """
        Return the currently selected sample annotations.

        The return value is a list of selected (sample type, sample value)
        tuples.

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
                # If no sample of sample type is selected we don't filter
                # on it.
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

            _, groups = self.selectedSamples()
            if len(groups) == 1 and self.outputRows:
                sample_type = groups[0]

            self.setEnabled(False)
            self.setBlocking(True)

            def get_data(gds_id, report_genes, transpose, sample_type, title):
                gds = geo.GDS(gds_id)
                data = gds.getdata(
                    report_genes=report_genes, transpose=transpose,
                    sample_type=sample_type
                )
                data.name = title
                return data

            get_data = partial(
                get_data, self.currentGds["dataset_id"],
                report_genes=self.mergeSpots,
                transpose=self.outputRows,
                sample_type=sample_type,
                title=self.datasetName or self.currentGds["title"]
            )
            self._datatask = Task(function=get_data)
            self._datatask.finished.connect(self._on_dataready)
            self._executor.submit(self._datatask)

    def _on_dataready(self):
        self.setEnabled(True)
        self.setBlocking(False)

        self.progressBarSet(50)

        try:
            data = self._datatask.result()
        except urllib2.URLError as error:
            self.error(0, "Error while connecting to the NCBI ftp server! %r" %
                       error)
            self._datatask = None
            self.progressBarFinished()
            return

        self._datatask = None

        data_name = data.name
        samples, _ = self.selectedSamples()

        self.warning(0)
        message = None
        if self.outputRows:
            def samplesinst(ex):
                out = []
                for i, a in data.domain.get_metas().items():
                    out.append((a.name, ex[i].value))
                if data.domain.class_var.name != 'class':
                    out.append((data.domain.class_var.name, ex[-1].value))
                return out
            samples = set(samples)

            select = [1 if samples.issuperset(samplesinst(ex)) else 0
                      for ex in data]
            data = data.select(select)
            if len(data) == 0:
                message = "No samples with selected sample annotations."
        else:
            samples = set(samples)
            domain = orange.Domain(
                [attr for attr in data.domain.attributes
                 if samples.issuperset(attr.attributes.items())],
                data.domain.classVar
            )
            domain.addmetas(data.domain.getmetas())
            if len(domain.attributes) == 0:
                message = "No samples with selected sample annotations."
            stypes = set(s[0] for s in samples)
            for attr in domain.attributes:
                attr.attributes = dict(
                    (key, value) for key, value in attr.attributes.items()
                    if key in stypes
                )
            data = orange.ExampleTable(domain, data)

        if message is not None:
            self.warning(0, message)

        data_hints.set_hint(data, "taxid", self.currentGds.get("taxid", ""),
                            10.0)
        data_hints.set_hint(data, "genesinrows", self.outputRows, 10.0)

        self.progressBarFinished()
        data.name = data_name
        self.send("Expression Data", data)

        model = self.treeWidget.model().sourceModel()
        row = self.gds.index(self.currentGds)

        model.setData(model.index(row, 0),  QVariant(" "), Qt.DisplayRole)

        self.updateInfo()
        self.selectionChanged = False

    def splitterMoved(self, *args):
        self.splitterSettings = [str(sp.saveState()) for sp in self.splitters]

    def onDeleteWidget(self):
        if self._inittask:
            self._inittask.future().cancel()
            self._inittask.finished.disconnect(self._initializemodel)
        if self._datatask:
            self._datatask.future().cancel()
            self._datatask.finished.disconnect(self._on_dataready)
        self._executor.shutdown(wait=False)

        super(OWGEODatasets, self).onDeleteWidget()

    def onNameEdited(self):
        if self.currentGds:
            gds_id = self.currentGds["dataset_id"]
            self.datasetNames[gds_id] = unicode(self.nameEdit.text())
            self.commitIf()

def get_gds_model(progress=lambda val: None):
    """
    Initialize and return a GDS datasets model.

    :param progress: A progress callback.
    :rval tuple:
        A tuple of (QStandardItemModel, geo.GDSInfo, [geo.GDS])

    .. note::
        The returned QStandardItemModel's thread affinity is set to
        the GUI thread.

    """
    progress(1)
    info = geo.GDSInfo()
    search_keys = ["dataset_id", "title", "platform_organism", "description"]
    cache_dir = serverfiles.localpath(geo.DOMAIN)
    gds_link = "http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc={0}"
    pm_link = "http://www.ncbi.nlm.nih.gov/pubmed/{0}"
    gds_list = []

    def is_cached(gds):
        return os.path.exists(os.path.join(cache_dir, gds["dataset_id"]) +
                              ".soft.gz")

    def item(displayvalue, item_values={}):
        item = QStandardItem()
        item.setData(displayvalue, Qt.DisplayRole)
        for role, value in item_values.iteritems():
            item.setData(value, role)
        return item

    def gds_to_row(gds):
        #: Text for easier full search.
        search_text = unicode(
            " | ".join([gds.get(key, "").lower()
                        for key in search_keys]),
            errors="ignore"
        )
        row = [
            item(" " if is_cached(gds) else "",
                 {TextFilterRole: search_text}),
            item(gds["dataset_id"],
                 {LinkRole: gds_link.format(gds["dataset_id"])}),
            item(gds["title"]),
            item(gds["platform_organism"]),
            item(len(gds["samples"])),
            item(gds["feature_count"]),
            item(gds["gene_count"]),
            item(len(gds["subsets"])),
            item(gds.get("pubmed_id", ""),
                 {LinkRole: pm_link.format(gds["pubmed_id"])
                            if gds.get("pubmed_id")
                            else QVariant()})
        ]
        return row

    model = QStandardItemModel()
    model.setHorizontalHeaderLabels(
        ["", "ID", "Title", "Organism", "Samples", "Features",
         "Genes", "Subsets", "PubMedID"]
    )
    progress(20)
    for gds in info.values():
        model.appendRow(gds_to_row(gds))

        gds_list.append(gds)

    progress(50)

    if QThread.currentThread() is not QCoreApplication.instance().thread():
        model.moveToThread(QCoreApplication.instance().thread())
    return model, info, gds_list


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWGEODatasets()
    w.show()
    app.exec_()
    w.saveSettings()
