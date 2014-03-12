"""
<name>Gene Info</name>
<description>Displays gene information from NCBI and other sources.</description>
<priority>2010</priority>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact>
<icon>icons/GeneInfo.svg</icon>
"""

from __future__ import absolute_import, with_statement

import sys
from collections import defaultdict
from functools import partial

from PyQt4.QtCore import pyqtSlot as Slot

import Orange

from Orange.utils import serverfiles
from Orange.utils import lru_cache

from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWGUI import LinkStyledItemDelegate, LinkRole

from Orange.OrangeWidgets.OWWidget import *

from Orange.OrangeWidgets.OWConcurrent import \
    ThreadExecutor, Task, methodinvoke


from .. import gene, taxonomy
from .utils import download


NAME = "Gene Info"
DESCRIPTION = "Displays gene information from NCBI and other sources."
ICON = "icons/GeneInfo.svg"
PRIORITY = 2010

INPUTS = [("Examples", Orange.data.Table, "setData")]
OUTPUTS = [("Selected Examples", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWGeneInfo.OWGeneInfo"]


class TreeModel(QAbstractItemModel):

    def __init__(self, data, header, parent):
        QAbstractItemModel.__init__(self, parent)
        self._data = [[QVariant(s) for s in row] for row in data]
        self._dataDict = {}
        self._header = header
        self._roleData = {Qt.DisplayRole: self._data}
        self._roleData = partial(
            defaultdict,
            partial(defaultdict,
                    partial(defaultdict, QVariant)))(self._roleData)

    def setColumnLinks(self, column, links):
        font = QFont()
        font.setUnderline(True)
        font = QVariant(font)
        for i, link in enumerate(links):
            self._roleData[LinkRole][i][column] = QVariant(link)
            self._roleData[Qt.FontRole][i][column] = font
            self._roleData[Qt.ForegroundRole][i][column] = \
                QVariant(QColor(Qt.blue))

    def setRoleData(self, role, row, col, data):
        self._roleData[role][row][col] = data

    def data(self, index, role=Qt.DisplayRole):
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

    def columnCount(self, index=QModelIndex()):
        return len(self._header)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return QVariant(self._header[section])
        return QVariant()


class LinkFmt(object):

    def __init__(self, link_fmt, name):
        self.link_fmt = link_fmt
        self.name = name

    def format(self, *args, **kwargs):
        return Link(self.link_fmt.format(*args, **kwargs), **kwargs)

    def __repr__(self):
        return "<LinkFmt " + repr(self.name) + " >"

    def __str__(self):
        return self.name


class Link(object):

    def __init__(self, link, text=None, **kwargs):
        self.link = link
        self.text = text if text is not None else "link"
        self.__dict__.update(kwargs)


@lru_cache(maxsize=2)
def get_ncbi_info(taxid):
    return gene.NCBIGeneInfo(taxid)


def ncbi_info(taxid, genes, advance=None):
    taxid = gene.NCBIGeneInfo.TAX_MAP.get(taxid, taxid)
    download.ensure_downloaded(
        "NCBI_geneinfo",
        "gene_info.%s.db" % taxid,
        advance
    )
    info = get_ncbi_info(taxid)

    schema_link = LinkFmt(
        "http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch={gene_id}",
        name="NCBI ID")

    schema = [schema_link, "Symbol", "Locus Tag", "Chromosome",
              "Description", "Synonyms", "Nomenclature"]
    ret = []
    for gene_name in genes:
        gi = info.get_info(gene_name)
        if gi:
            ret.append([schema_link.format(gene_id=gi.gene_id, text=gi.gene_id),
                        gi.symbol + " (%s)" % gene_name if gene_name != gi.symbol else gi.symbol,
                        gi.locus_tag or "",
                        gi.chromosome or "",
                        gi.description or "",
                        ", ".join(gi.synonyms),
                        gi.symbol_from_nomenclature_authority or ""
                        ])
        else:
            ret.append(None)
    return schema, ret


def dicty_info(taxid, genes, advance=None):
    from .. import dicty
    download.ensure_downloaded(
        dicty.DictyBase.domain,
        dicty.DictyBase.filename,
        advance
    )
    info = dicty.DictyBase()
    name_matcher = gene.GMDicty()
    name_matcher.set_targets(info.info.keys())
    schema_link = LinkFmt(
        "http://dictybase.org/db/cgi-bin/gene_page.pl?dictybaseid={gene_id}",
        name="Dicty Base Id")
    schema = [schema_link, "Name", "Synonyms", "Gene Products"]

    ret = []
    for gene_name in genes:
        gene_name = name_matcher.umatch(gene_name)
        gi = info.info.get(gene_name, None)
        if gi:
            ret.append([schema_link.format(gene_id=gene_name, text=gene_name),
                        gi[0] + " (%s)" % gene_name if gene_name != gi[0] else gi[0],  # Gene Name
                        ", ".join(gi[1]),  # Synonyms
                        gi[2] or "",  # Gene Products
                        ])

        else:
            ret.append(None)

    return schema, ret


INFO_SOURCES = {
    "default": [("NCBI Info", ncbi_info)],
    "352472": [("NCBI Info", ncbi_info),
               ("Dicty Base", dicty_info)]
}


class OWGeneInfo(OWWidget):
    settingsList = ["organismIndex", "geneAttr", "useAttr", "autoCommit",
                    "taxid"]
    contextHandlers = {
        "": DomainContextHandler(
            "", ["organismIndex", "geneAttr", "useAttr", "useAltSource",
                 "taxid"]
        )
    }

    def __init__(self, parent=None, signalManager=None, name="Gene Info"):
        OWWidget.__init__(self, parent, signalManager, name)

        self.inputs = [("Examples", Orange.data.Table, self.setData)]
        self.outputs = [("Selected Examples", Orange.data.Table)]

        self.organismIndex = 0
        self.taxid = None
        self.geneAttr = 0
        self.useAttr = False
        self.autoCommit = False
        self.searchString = ""
        self.selectionChangedFlag = False
        self.useAltSource = 0
        self.loadSettings()

        self.__initialized = False
        self.initfuture = None
        self.itemsfuture = None

        self.infoLabel = OWGUI.widgetLabel(
            OWGUI.widgetBox(self.controlArea, "Info", addSpace=True),
            "Initializing\n"
        )

        self.organisms = None
        self.organismBox = OWGUI.widgetBox(
            self.controlArea, "Organism", addSpace=True)

        self.organismComboBox = OWGUI.comboBox(
            self.organismBox, self, "organismIndex",
            callback=self._onSelectedOrganismChanged,
            debuggingEnabled=0)

        # For now only support one alt source, with a checkbox
        # In the future this can be extended to multiple selections
        self.altSourceCheck = OWGUI.checkBox(self.organismBox, self,
                            "useAltSource", "Show information from dictyBase",
                            callback=self.onAltSourceChange,
#                            debuggingEnabled=0,
                            )
        self.altSourceCheck.hide()

        box = OWGUI.widgetBox(self.controlArea, "Gene names", addSpace=True)
        self.geneAttrComboBox = OWGUI.comboBox(
            box, self, "geneAttr",
            "Gene atttibute", callback=self.updateInfoItems
        )
        OWGUI.checkBox(box, self, "useAttr", "Use attribute names",
                       callback=self.updateInfoItems,
                       disables=[(-1, self.geneAttrComboBox)])

        self.geneAttrComboBox.setDisabled(bool(self.useAttr))

        box = OWGUI.widgetBox(self.controlArea, "Commit", addSpace=True)
        b = OWGUI.button(box, self, "Commit", callback=self.commit)
        c = OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.setStopper(self, b, c, "selectionChangedFlag",
                         callback=self.commit)

        # A label for dictyExpress link
        self.dictyExpressBox = OWGUI.widgetBox(
            self.controlArea, "Dicty Express")
        self.linkLabel = OWGUI.widgetLabel(self.dictyExpressBox, "")
        self.linkLabel.setOpenExternalLinks(False)
        self.connect(self.linkLabel, SIGNAL("linkActivated(QString)"),
                     self.onDictyExpressLink)
        self.dictyExpressBox.hide()

        OWGUI.rubber(self.controlArea)

        OWGUI.lineEdit(self.mainArea, self, "searchString", "Filter",
                       callbackOnType=True, callback=self.searchUpdate)

        self.treeWidget = QTreeView(self.mainArea)
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSelectionMode(
            QAbstractItemView.ExtendedSelection)
        self.treeWidget.setItemDelegate(
            LinkStyledItemDelegate(self.treeWidget))
        self.treeWidget.setUniformRowHeights(True)
        self.treeWidget.viewport().setMouseTracking(True)
        self.treeWidget.setSortingEnabled(True)
        self.mainArea.layout().addWidget(self.treeWidget)

        box = OWGUI.widgetBox(self.mainArea, "",
                              orientation="horizontal")
        OWGUI.button(box, self, "Select Filtered",
                     callback=self.selectFiltered)
        OWGUI.button(box, self, "Clear Selection",
                     callback=self.treeWidget.clearSelection)

        self.resize(1000, 700)

        self.geneinfo = []
        self.cells = []
        self.row2geneinfo = {}
        self.data = None

        # : (# input genes, # matches genes)
        self.matchedInfo = 0, 0
        self.selectionUpdateInProgress = False

        self.setBlocking(True)
        self.executor = ThreadExecutor(self)

        self.progressBarInit()

        task = Task(
            function=partial(
                taxonomy.ensure_downloaded,
                callback=methodinvoke(self, "advance", ())
            )
        )

        task.resultReady.connect(self.initialize)
        task.exceptionReady.connect(self._onInitializeError)

        self.initfuture = self.executor.submit(task)

    @Slot()
    def advance(self):
        assert self.thread() is QThread.currentThread()

        self.progressBarSet(self.progressBarValue + 1,
                            processEventsFlags=None)

    def initialize(self):
        if self.__initialized:
            # Already initialized
            return

        self.progressBarFinished()

        self.organisms = sorted(
            set([name.split(".")[-2] for name in
                 serverfiles.listfiles("NCBI_geneinfo")] +
                gene.NCBIGeneInfo.essential_taxids())
        )

        self.organismComboBox.addItems(
            [taxonomy.name(tax_id) for tax_id in self.organisms]
        )
        if self.taxid in self.organisms:
            self.organismIndex = self.organisms.index(self.taxid)

        self.infoLabel.setText("No data on input\n")
        self.__initialized = True
        self.initfuture = None

        self.setBlocking(False)

    def _onInitializeError(self, exc):
        sys.excepthook(type(exc), exc.args, None)
        self.error(0, "Could not download the necessary files.")

    def _onSelectedOrganismChanged(self):
        self.taxid = self.organisms[self.organismIndex]
        if self.data is not None:
            self.updateInfoItems()

    def setData(self, data=None):
        if not self.__initialized:
            self.initfuture.result()
            self.initialize()

        if self.itemsfuture is not None:
            raise Exception("Already processing")

        self.closeContext()
        self.data = data

        if data:
            self.geneAttrComboBox.clear()
            self.attributes = \
                [attr for attr in (data.domain.variables +
                                   data.domain.getmetas().values())
                 if isinstance(attr, (Orange.feature.String,
                                      Orange.feature.Discrete))]

            self.geneAttrComboBox.addItems(
                [attr.name for attr in self.attributes]
            )

            self.taxid = data_hints.get_hint(self.data, "taxid", self.taxid)
            self.useAttr = data_hints.get_hint(
                self.data, "genesinrows", self.useAttr)

            self.openContext("", data)
            self.geneAttr = min(self.geneAttr, len(self.attributes) - 1)

            if self.taxid in self.organisms:
                self.organismIndex = self.organisms.index(self.taxid)

            self.updateInfoItems()
        else:
            self.clear()

    def infoSource(self):
        """ Return the current selected info source getter function from
        INFO_SOURCES
        """
        org = self.organisms[min(self.organismIndex, len(self.organisms) - 1)]
        if org not in INFO_SOURCES:
            org = "default"
        sources = INFO_SOURCES[org]
        name, func = sources[min(self.useAltSource, len(sources) - 1)]
        return name, func

    def inputGenes(self):
        if self.useAttr:
            genes = [attr.name for attr in self.data.domain.attributes]
        elif self.attributes:
            attr = self.attributes[self.geneAttr]
            genes = [str(ex[attr]) for ex in self.data
                     if not ex[attr].isSpecial()]
        else:
            genes = []
        return genes

    def updateInfoItems(self):
        self.warning(0)
        if not self.data:
            return

        genes = self.inputGenes()
        if self.useAttr:
            genes = [attr.name for attr in self.data.domain.attributes]
        elif self.attributes:
            attr = self.attributes[self.geneAttr]
            genes = [str(ex[attr]) for ex in self.data
                     if not ex[attr].isSpecial()]
        else:
            genes = []
        if not genes:
            self.warning(0, "Could not extract genes from input dataset.")

        self.warning(1)
        org = self.organisms[min(self.organismIndex, len(self.organisms) - 1)]
        source_name, info_getter = self.infoSource()

        self.error(0)

        self.updateDictyExpressLink(genes, show=org == "352472")
        self.altSourceCheck.setVisible(org == "352472")

        self.progressBarInit()
        self.setBlocking(True)
        self.setEnabled(False)
        self.infoLabel.setText("Retrieving info records.\n")

        self.genes = genes

        task = Task(
            function=partial(
                info_getter, org, genes,
                advance=methodinvoke(self, "advance", ()))
        )
        self.itemsfuture = self.executor.submit(task)
        task.finished.connect(self._onItemsCompleted)

    def _onItemsCompleted(self):
        self.setBlocking(False)
        self.progressBarFinished()
        self.setEnabled(True)

        try:
            schema, geneinfo = self.itemsfuture.result()
        finally:
            self.itemsfuture = None

        self.geneinfo = geneinfo = list(zip(self.genes, geneinfo))
        self.cells = cells = []
        self.row2geneinfo = {}
        links = []
        for i, (_, gi) in enumerate(geneinfo):
            if gi:
                row = []
                for _, item in zip(schema, gi):
                    if isinstance(item, Link):
                        # TODO: This should be handled by delegates
                        row.append(item.text)
                        links.append(item.link)
                    else:
                        row.append(item)
                cells.append(row)
                self.row2geneinfo[len(cells) - 1] = i

        model = TreeModel(cells, [str(col) for col in schema], None)

        model.setColumnLinks(0, links)
        proxyModel = QSortFilterProxyModel(self)
        proxyModel.setSourceModel(model)
        self.treeWidget.setModel(proxyModel)
        self.connect(self.treeWidget.selectionModel(),
                     SIGNAL("selectionChanged(QItemSelection , QItemSelection )"),
                     self.commitIf)

        for i in range(7):
            self.treeWidget.resizeColumnToContents(i)
            self.treeWidget.setColumnWidth(
                i, min(self.treeWidget.columnWidth(i), 200)
            )

        self.infoLabel.setText("%i genes\n%i matched NCBI's IDs" %
                               (len(self.genes), len(cells)))
        self.matchedInfo = len(self.genes), len(cells)

    def clear(self):
        self.infoLabel.setText("No data on input\n")
        self.treeWidget.setModel(
            TreeModel([], ["NCBI ID", "Symbol", "Locus Tag",
                           "Chromosome", "Description", "Synonyms",
                           "Nomenclature"], self.treeWidget))

        self.geneAttrComboBox.clear()
        self.send("Selected Examples", None)

    def commitIf(self, *args):
        if self.autoCommit and not self.selectionUpdateInProgress:
            self.commit()
        else:
            self.selectionChangedFlag = True

    def commit(self):
        if not self.data:
            return
        model = self.treeWidget.model()
        mapToSource = model.mapToSource
        selectedRows = self.treeWidget.selectedIndexes()
        selectedRows = [mapToSource(index).row() for index in selectedRows]
        model = model.sourceModel()

        selectedGeneids = [self.row2geneinfo[row] for row in selectedRows]
        selectedIds = [self.geneinfo[i][0] for i in selectedGeneids]
        selectedIds = set(selectedIds)
        gene2row = dict((self.geneinfo[self.row2geneinfo[row]][0], row)
                        for row in selectedRows)

        if self.useAttr:
            def is_selected(attr):
                return attr.name in selectedIds
            attrs = [attr for attr in self.data.domain.attributes
                     if is_selected(attr)]
            domain = Orange.data.Domain(attrs, self.data.domain.classVar)
            domain.addmetas(self.data.domain.getmetas())
            newdata = Orange.data.Table(domain, self.data)
            self.send("Selected Examples", newdata)
        elif self.attributes:
            attr = self.attributes[self.geneAttr]
            examples = [ex for ex in self.data if str(ex[attr]) in selectedIds]
            # Add gene info
            domain = Orange.data.Domain(
                self.data.domain, self.data.domain.classVar)
            domain.addmetas(self.data.domain.getmetas())
            n_columns = model.columnCount()

            headers = [str(model.headerData(i, Qt.Horizontal, Qt.DisplayRole)
                           .toString())
                       for i in range(n_columns)]
            new_meta_attrs = [(Orange.feature.Descriptor.new_meta_id(),
                               Orange.feature.String(name))
                              for name in headers]
            domain.addmetas(dict(new_meta_attrs))
            examples = [Orange.data.Instance(domain, ex) for ex in examples]
            for ex in examples:
                for i, (_, meta) in enumerate(new_meta_attrs):
                    index = model.index(gene2row[str(ex[attr])], i)
                    ex[meta] = str(
                        model.data(index, Qt.DisplayRole).toString()
                    )

            if examples:
                newdata = Orange.data.Table(examples)
            else:
                newdata = None
            self.send("Selected Examples", newdata)
        else:
            self.send("Selected Examples", None)

    def rowFiltered(self, row):
        searchStrings = self.searchString.lower().split()
        row = unicode(" ".join(self.cells[row]).lower(), errors="ignore")
        return not all([s in row for s in searchStrings])

    def searchUpdate(self):
        if not self.data:
            return
        searchStrings = self.searchString.lower().split()
        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            row = unicode(" ".join(row).lower(), errors="ignore")
            self.treeWidget.setRowHidden(
                mapFromSource(index(i, 0)).row(),
                QModelIndex(),
                not all([s in row for s in searchStrings]))

    def selectFiltered(self):
        if not self.data:
            return
        itemSelection = QItemSelection()

        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            if not self.rowFiltered(i):
                itemSelection.select(mapFromSource(index(i, 0)),
                                     mapFromSource(index(i, 0)))
        self.treeWidget.selectionModel().select(
            itemSelection,
            QItemSelectionModel.Select | QItemSelectionModel.Rows)

    def sendReport(self):
        from Orange.OrangeWidgets import OWReport
        genes, matched = self.matchedInfo

        if self.organisms:
            org = self.organisms[min(self.organismIndex,
                                     len(self.organisms) - 1)]
            org_name = taxonomy.name(org)
        else:
            org = None
            org_name = None
        if self.data is not None:
            self.reportRaw(
                "<p>Input: %i genes of which %i (%.1f%%) matched NCBI synonyms"
                "<br>"
                "Organism: %s"
                "<br>"
                "Filter: %s"
                "</p>" % (genes, matched, 100.0 * matched / genes, org_name,
                          self.searchString)
            )
            self.reportSubsection("Gene list")
            self.reportRaw(reportItemView(self.treeWidget))
        else:
            self.reportRaw("<p>No input</p>")

    def updateDictyExpressLink(self, genes, show=False):
        def fix(ddb):
            if ddb.startswith("DDB"):
                if not ddb.startswith("DDB_G"):
                    ddb = ddb.replace("DDB", "DDB_G")
                return ddb
            return None
        if show:
            genes = [fix(gene) for gene in genes if fix(gene)]
            link1 = '<a href="http://dictyexpress.biolab.si/run/index.php?gene=%s">Microarray profile</a>'
            link2 = '<a href="http://dictyexpress.biolab.si/run/index.php?gene=%s&db=rnaseq">RNA-Seq profile</a>'
            self.linkLabel.setText(link1 + "<br/>" + link2)

            show = any(genes)

        if show:
            self.dictyExpressBox.show()
        else:
            self.dictyExpressBox.hide()

    def onDictyExpressLink(self, link):
        if not self.data:
            return

        selectedIndexes = self.treeWidget.selectedIndexes()
        if not len(selectedIndexes):
            QMessageBox.information(
                self, "No gene ids selected",
                "Please select some genes and try again."
            )
            return
        model = self.treeWidget.model()
        mapToSource = model.mapToSource
        selectedRows = self.treeWidget.selectedIndexes()
        selectedRows = [mapToSource(index).row() for index in selectedRows]
        model = model.sourceModel()

        selectedGeneids = [self.row2geneinfo[row] for row in selectedRows]
        selectedIds = [self.geneinfo[i][0] for i in selectedGeneids]
        selectedIds = set(selectedIds)

        def fix(ddb):
            if ddb.startswith("DDB"):
                if not ddb.startswith("DDB_G"):
                    ddb = ddb.replace("DDB", "DDB_G")
                return ddb
            return None

        genes = [fix(gene) for gene in selectedIds if fix(gene)]
        url = str(link) % " ".join(genes)
        QDesktopServices.openUrl(QUrl(url))

    def onAltSourceChange(self):
        self.updateInfoItems()

    def onDeleteWidget(self):
        OWWidget.onDeleteWidget(self)

        # try to cancel pending tasks
        if self.initfuture:
            self.initfuture.cancel()
        if self.itemsfuture:
            self.itemsfuture.cancel()

        self.executor.shutdown()


def reportItemView(view):
    model = view.model()
    return reportItemModel(view, model)


def reportItemModel(view, model, index=QModelIndex()):
    if not index.isValid() or model.hasChildren(index):
        columnCount, rowCount = model.columnCount(index), model.rowCount(index)
        if not index.isValid():
            text = ('<table>\n<tr>' +
                    ''.join('<th>%s</th>' %
                            model.headerData(i, Qt.Horizontal, Qt.DisplayRole)
                            .toString()
                            for i in range(columnCount)) +
                    '</tr>\n')
        else:
            pass
        text += ''.join('<tr>' +
                        ''.join('<td>' + reportItemModel(view, model, model.index(row, column, index)) +
                                '</td>' for column in range(columnCount)) +
                        '</tr>\n'
                        for row in range(rowCount)
                        if not view.isRowHidden(row, index))
        text += '</table>'
        return text
    else:
        variant = model.data(index, Qt.DisplayRole)
        return str(variant.toString()) if variant.isValid() else ""


if __name__ == "__main__":
    app = QApplication(sys.argv)
    data = Orange.data.Table("brown-selected.tab")
    w = OWGeneInfo()
    w.show()

    w.setData(data)
    app.exec_()
    w.saveSettings()
