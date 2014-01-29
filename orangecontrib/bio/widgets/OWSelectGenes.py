import sys
import os
import re
import unicodedata
import operator
from collections import defaultdict, namedtuple
from operator import itemgetter
from xml.sax.saxutils import escape

from contextlib import contextmanager

from PyQt4.QtGui import (
    QFrame, QHBoxLayout, QPlainTextEdit, QSyntaxHighlighter, QTextCharFormat,
    QTextCursor, QCompleter, QStandardItemModel, QSortFilterProxyModel,
    QStandardItem, QListView, QTreeView, QStyle, QStyledItemDelegate,
    QStyleOptionViewItemV4, QPalette, QColor, QApplication, QAction,
    QToolButton, QItemSelectionModel, QPlainTextDocumentLayout, QTextDocument,
    QRadioButton, QButtonGroup, QStyleOptionButton, QMenu, QDialog,
    QStackedWidget, QComboBox, QFileDialog
)

from PyQt4.QtCore import Qt, QEvent, QVariant, QThread
from PyQt4.QtCore import pyqtSignal as Signal

import Orange

from Orange.OrangeWidgets.OWWidget import (
    OWWidget, DomainContextHandler, Default
)

from Orange.OrangeWidgets.OWConcurrent import (
    ThreadExecutor, Task, methodinvoke
)

from Orange.OrangeWidgets.OWItemModels import VariableListModel
from Orange.OrangeWidgets import OWGUI
from Orange.orng.orngDataCaching import data_hints
from Orange.bio import obiGene as geneinfo
from Orange.bio import obiTaxonomy as taxonomy


NAME = "Select Genes"
DESCRIPTION = "Select a specified subset of the input genes."
ICON = "icons/SelectGenes.svg"

INPUTS = [("Data", Orange.data.Table, "setData", Default),
          ("Gene Subset", Orange.data.Table, "setGeneSubset")]

OUTPUTS = [("Selected Data", Orange.data.Table)]


def toString(variant):
    if isinstance(variant, QVariant):
        return unicode(variant.toString())
    else:
        return unicode(variant)


def toBool(variant):
    if isinstance(variant, QVariant):
        return bool(variant.toPyObject())
    else:
        return bool(variant)


def toPyObject(variant):
    if isinstance(variant, QVariant):
        return variant.toPyObject()
    else:
        return variant


class SaveSlot(QStandardItem):
    ModifiedRole = next(OWGUI.OrangeUserRole)

    def __init__(self, name, savedata=None, modified=False):
        super(SaveSlot, self).__init__(name)

        self.savedata = savedata
        self.modified = modified
        self.document = None

    @property
    def name(self):
        return unicode(self.text())

    @property
    def modified(self):
        return toBool(self.data(SaveSlot.ModifiedRole))

    @modified.setter
    def modified(self, state):
        self.setData(bool(state), SaveSlot.ModifiedRole)


class SavedSlotDelegate(QStyledItemDelegate):

    def paint(self, painter, option, index):
        option = QStyleOptionViewItemV4(option)
        self.initStyleOption(option, index)

        modified = toBool(index.data(SaveSlot.ModifiedRole))
        if modified:
            option.palette.setColor(QPalette.Text, QColor(Qt.red))
            option.palette.setColor(QPalette.Highlight, QColor(Qt.darkRed))
            option.text = "*" + option.text

        if option.widget:
            widget = option.widget
            style = widget.style()
        else:
            widget = None
            style = QApplication.style()

        style.drawControl(QStyle.CE_ItemViewItem, option, painter, widget)


def radio_indicator_width(button):
    button.ensurePolished()
    style = button.style()
    option = QStyleOptionButton()
    button.initStyleOption(option)

    w = style.pixelMetric(QStyle.PM_ExclusiveIndicatorWidth, option, button)
    return w


class OWSelectGenes(OWWidget):

    contextHandlers = {
        "": DomainContextHandler(
            "", ["geneIndex", "taxid"]
        ),
        "subset": DomainContextHandler(
            "subset", ["subsetGeneIndex"]
        )
    }

    settingsList = ["autoCommit", "preserveOrder", "savedSelections",
                    "selectedSelectionIndex", "selectedSource",
                    "completeOnSymbols"]

    SelectInput, SelectCustom = 0, 1
    CompletionRole = Qt.UserRole + 1

    def __init__(self, parent=None, signalManager=None, title=NAME):
        OWWidget.__init__(self, parent, signalManager, title,
                          wantMainArea=False)

        self.geneIndex = None
        self.taxid = None
        self.autoCommit = False
        self.preserveOrder = True
        self.savedSelections = [
            ("Example", ["MRE11A", "RAD51", "MLH1", "MSH2", "DMC1"])
        ]

        self.selectedSelectionIndex = -1
        self.selectedSource = OWSelectGenes.SelectCustom
        self.completeOnSymbols = True

        self.loadSettings()

        # Input variables that could contain names
        self.variables = VariableListModel()
        # All gene names and their symbols
        self.geneNames = []
        # A list of (name, info) where name is from the input
        # (geneVar column) and info is the NCBIGeneInfo object if available
        # or None
        self.genes = []
        # Output changed flag
        self._changedFlag = False
        # Current gene names
        self.selection = []
        # Input data
        self.data = None
        self.subsetData = None
        # Input variables that could contain gene names from "Gene Subset"
        self.subsetVariables = VariableListModel()
        # Selected subset variable index
        self.subsetGeneIndex = -1
        self.organisms = []
        self.taxidindex = {}
        self.geneinfo = (None, None)
        self._executor = ThreadExecutor()

        self._infotask = None

        box = OWGUI.widgetBox(self.controlArea, "Gene Attribute")
        box.setToolTip("Column with gene names")
        self.attrsCombo = OWGUI.comboBox(
            box, self, "geneIndex",
            callback=self._onGeneIndexChanged,
        )
        self.attrsCombo.setModel(self.variables)

        box = OWGUI.widgetBox(self.controlArea, "Gene Selection")

        button1 = QRadioButton("Select genes from 'Gene Subset' input")
        button2 = QRadioButton("Select specified genes")

        box.layout().addWidget(button1)

        # Subset gene variable selection
        self.subsetbox = OWGUI.widgetBox(box, None)
        offset = radio_indicator_width(button1)
        self.subsetbox.layout().setContentsMargins(offset, 0, 0, 0)
        self.subsetbox.setEnabled(
            self.selectedSource == OWSelectGenes.SelectInput)

        box1 = OWGUI.widgetBox(self.subsetbox, "Gene Attribute", flat=True)
        self.subsetVarCombo = OWGUI.comboBox(
            box1, self, "subsetGeneIndex",
            callback=self._onSubsetGeneIndexChanged
        )
        self.subsetVarCombo.setModel(self.subsetVariables)
        self.subsetVarCombo.setToolTip(
            "Column with gene names in the 'Gene Subset' input"
        )
        OWGUI.button(box1, self, "Copy genes to saved subsets",
                     callback=self.copyToSaved)

        OWGUI.button(box1, self, "Append genes to current saved selection",
                     callback=self.appendToSelected)

        box.layout().addWidget(button2)

        self.selectedSourceButtons = group = QButtonGroup(box)
        group.addButton(button1, OWSelectGenes.SelectInput)
        group.addButton(button2, OWSelectGenes.SelectCustom)
        group.buttonClicked[int].connect(self._selectionSourceChanged)

        if self.selectedSource == OWSelectGenes.SelectInput:
            button1.setChecked(True)
        else:
            button2.setChecked(True)

        self.entrybox = OWGUI.widgetBox(box, None)
        offset = radio_indicator_width(button2)
        self.entrybox.layout().setContentsMargins(offset, 0, 0, 0)

        self.entrybox.setEnabled(
            self.selectedSource == OWSelectGenes.SelectCustom)

        box = OWGUI.widgetBox(self.entrybox, "Select Genes", flat=True)
        box.setToolTip("Enter gene names to select")
        box.layout().setSpacing(1)

        self.entryField = ListTextEdit(box)
        self.entryField.setTabChangesFocus(True)
        self.entryField.setDocument(self._createDocument())
        self.entryField.itemsChanged.connect(self._onItemsChanged)

        box.layout().addWidget(self.entryField)

        completer = ListCompleter()
        completer.setCompletionMode(QCompleter.PopupCompletion)
        completer.setCompletionRole(
            self.CompletionRole if self.completeOnSymbols else Qt.DisplayRole
        )
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        completer.setMaxVisibleItems(10)

        popup = QTreeView()
        popup.setSelectionMode(QTreeView.ExtendedSelection)
        popup.setEditTriggers(QTreeView.NoEditTriggers)
        popup.setRootIsDecorated(False)
        popup.setAlternatingRowColors(True)
        popup.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        popup.setMaximumWidth(500)
        popup.header().setStretchLastSection(False)
        popup.header().hide()

        completer.setPopup(popup)
        completer.setModel(SetFilterProxyModel(self))

        self.entryField.setCompleter(completer)

        toolbar = QFrame()
        toolbar.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Fixed)
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(1)

        addaction = QAction("+", self)
        addmenu = QMenu()

        action = addmenu.addAction("Import names from gene sets...")
        action.triggered.connect(self.importGeneSet)

        action = addmenu.addAction("Import names from text file...")
        action.triggered.connect(self.importFromFile)

        addaction.setMenu(addmenu)

        def button(action, popupMode=None):
            b = QToolButton()
            b.setDefaultAction(action)
            if popupMode is not None:
                b.setPopupMode(QToolButton.InstantPopup)
            return b

        b = button(addaction, popupMode=QToolButton.InstantPopup)
        layout.addWidget(b)

        moreaction = QAction("More", self)

        moremenu = QMenu()
        self.completeOnSymbolsAction = QAction(
            "Complete on gene symbol names", self,
            toolTip="Use symbol names for auto completion.",
            checkable=True,
            checked=self.completeOnSymbols
        )

        self.completeOnSymbolsAction.toggled[bool].connect(
            self._onToggleSymbolCompletion
        )

        moremenu.addAction(self.completeOnSymbolsAction)

        self.translateAction = QAction(
            "Translate all names to official symbol names", self,
            enabled=False
        )
        self.translateAction.triggered.connect(self._onTranslate)

        moremenu.addAction(self.translateAction)
        moreaction.setMenu(moremenu)

        b = button(moreaction, popupMode=QToolButton.InstantPopup)
        layout.addWidget(b)

        self.organismsCombo = ComboBox(focusPolicy=Qt.StrongFocus)
        self.organismsCombo.addItem("...")
        self.organismsCombo.model().item(0).setEnabled(False)
        self.organismsCombo.setMinimumWidth(200)
        self.organismsCombo.activated[int].connect(self._onOrganismActivated)

        layout.addSpacing(10)
        layout.addWidget(self.organismsCombo)
        self.pb = QProgressBar()
        self.pb.hide()
        self.pbstack = QStackedWidget()
        self.pbstack.addWidget(self.pb)

        layout.addSpacing(10)
        layout.addWidget(self.pbstack, 0, Qt.AlignRight | Qt.AlignVCenter)
        toolbar.setLayout(layout)

        box.layout().addWidget(toolbar)

        box = OWGUI.widgetBox(self.entrybox, "Saved Selections", flat=True)
        box.setToolTip("Save/Select/Update saved gene selections")
        box.layout().setSpacing(1)

        self.selectionsModel = QStandardItemModel()
        self.selectionsView = QListView()
        self.selectionsView.setAlternatingRowColors(True)
        self.selectionsView.setModel(self.selectionsModel)
        self.selectionsView.setItemDelegate(SavedSlotDelegate(self))
        self.selectionsView.selectionModel().selectionChanged.connect(
            self._onSelectedSaveSlotChanged
        )

        box.layout().addWidget(self.selectionsView)

        self.actionSave = QAction(
            "Save", self,
            toolTip="Save/Update the current selection")

        self.actionAdd = QAction(
            "+", self,
            toolTip="Create a new saved selection")

        self.actionRemove = QAction(
            u"\u2212", self,
            toolTip="Delete the current saved selection")

        toolbar = QFrame()
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(1)

        b = button(self.actionAdd)
        layout.addWidget(b)

        b = button(self.actionRemove)
        layout.addWidget(b)

        b = button(self.actionSave)
        layout.addWidget(b)

        layout.addStretch(100)
        toolbar.setLayout(layout)

        box.layout().addWidget(toolbar)

        self.actionSave.triggered.connect(self.saveSelection)
        self.actionAdd.triggered.connect(self.addSelection)
        self.actionRemove.triggered.connect(self.removeSelection)

        box = OWGUI.widgetBox(self.controlArea, "Output")
        OWGUI.checkBox(box, self, "preserveOrder", "Preserve input order",
                       tooltip="Preserve the order of the input data "
                               "instances.",
                       callback=self.invalidateOutput)
        cb = OWGUI.checkBox(box, self, "autoCommit", "Auto commit")
        button = OWGUI.button(box, self, "Commit", callback=self.commit)

        OWGUI.setStopper(self, button, cb, "_changedFlag", self.commit)

        # Gene set import dialog (initialized when required)
        self._genesetDialog = None

        # restore saved selections model.
        for name, names in self.savedSelections:
            item = SaveSlot(name, names)
            self.selectionsModel.appendRow([item])

        if self.selectedSelectionIndex != -1:
            self.selectionsView.selectionModel().select(
                self.selectionsModel.index(self.selectedSelectionIndex, 0),
                QItemSelectionModel.Select
            )

        self._updateActions()

        def list_available_organisms():
            return [(taxid, taxonomy.name(taxid))
                    for taxid in taxonomy.common_taxids()]

        self._inittask = Task(function=list_available_organisms)

        self._inittask.finished.connect(self._onInit)
        self._executor.submit(self._inittask)
        self.pb.show()
        self.pb.setRange(0, 0)

    def _onInit(self):
        self.organisms = self._inittask.result()
        self.organismsCombo.clear()
        self.organismsCombo.addItems([name for _, name in self.organisms])

        self.taxidindex = \
            {taxid: i for i, (taxid, _) in enumerate(self.organisms)}

        self.organismsCombo.setCurrentIndex(
            self.taxidindex.get(self.taxid, -1))

        self.pb.hide()

        if self.taxid in self.taxidindex:
            self._updateGeneInfo()

    def _updateGeneInfo(self):
        self.geneinfo = (self.taxid, None)
        self.translateAction.setEnabled(False)

        if self._infotask:
            # Try to cancel existing pending task
            self._infotask.future().cancel()
            self._infotask.resultReady.disconnect(self._onGeneInfoReady)
            self._infotask.exceptionReady.disconnect(self._onGeneInfoError)

        def ncbi_gene_info(taxid=self.taxid):
            try:
                return (taxid, geneinfo.NCBIGeneInfo(taxid))
            except BaseException:
                sys.excepthook(*sys.exc_info())
                raise

        task = Task(function=ncbi_gene_info)
        task.resultReady.connect(self._onGeneInfoReady)
        task.exceptionReady.connect(self._onGeneInfoError)
        self._infotask = task
        self._executor.submit(task)
        self.pb.show()

    def _onOrganismActivated(self, index):
        try:
            taxid, _ = self.organisms[index]
        except IndexError:
            pass
        else:
            if taxid != self.taxid:
                self.taxid = taxid
                self._updateGeneInfo()

    def setData(self, data):
        """
        Set the input data.
        """
        self.closeContext("")
        self.warning(0)
        self.data = data
        if data is not None:
            attrs = gene_candidates(data)
            self.variables[:] = attrs
            self.attrsCombo.setCurrentIndex(0)
            if attrs:
                self.geneIndex = 0
            else:
                self.geneIndex = -1
                self.warning(0, "No suitable columns for gene names.")
        else:
            self.variables[:] = []
            self.geneIndex = -1

        self._changedFlag = True

        oldtaxid = self.taxid
        if data is not None:
            self.taxid = data_hints.get_hint(data, "taxid", None)
        else:
            self.taxid = None

        self.openContext("", data)

        if self.taxid is None:
            self.geneinfo = (None, None)

        if oldtaxid != self.taxid:
            self.organismsCombo.setCurrentIndex(
                self.taxidindex.get(self.taxid, -1))

            if self.taxid in self.taxidindex:
                self._updateGeneInfo()

        self._updateCompletionModel()

        self.commit()

    def setGeneSubset(self, data):
        """
        Set the gene subset input.
        """
        self.closeContext("subset")
        self.warning(1)
        self.subsetData = data
        if data is not None:
            variables = gene_candidates(data)
            self.subsetVariables[:] = variables
            self.subsetVarCombo.setCurrentIndex(0)
            if variables:
                self.subsetGeneIndex = 0
            else:
                self.subsetGeneIndex = -1
                self.warning(1, "No suitable column for subset gene names.")
        else:
            self.subsetVariables[:] = []
            self.subsetGeneIndex = -1

        self.openContext("subset", data)

        if self.selectedSource == OWSelectGenes.SelectInput:
            self.commit()

    @property
    def geneVar(self):
        """
        Current gene attribute or None if none available.
        """
        index = self.attrsCombo.currentIndex()
        if self.data is not None and index >= 0:
            return self.variables[index]
        else:
            return None

    @property
    def subsetGeneVar(self):
        """
        Current subset gene attribute or None if not available.
        """
        index = self.subsetVarCombo.currentIndex()
        if self.subsetData is not None and index >= 0:
            return self.subsetVariables[index]
        else:
            return None

    def invalidateOutput(self):
        if self.autoCommit:
            self.commit()
        else:
            self._changedFlag = True

    def selectedGenes(self):
        """
        Return the names of the current selected genes.
        """
        selection = []
        if self.selectedSource == OWSelectGenes.SelectInput:
            var = self.subsetGeneVar
            if var is not None:
                values = [inst[var] for inst in self.subsetData]
                selection = [str(val) for val in values
                             if not val.is_special()]
        else:
            selection = self.selection
        return selection

    def commit(self):
        """
        Send the selected data subset to the output.
        """
        selection = self.selectedGenes()

        if self.geneinfo[1] is not None:
            backmap = dict((info.symbol, name) for name, info in self.genes
                           if info is not None)
            names = set([name for name, _ in self.genes])
            selection = [backmap.get(name, name) for name in selection]
    
        if self.geneVar is not None:
            data = select_by_genes(self.data, self.geneVar,
                                   gene_list=selection,
                                   preserve_order=self.preserveOrder)
        else:
            data = None

        self.send("Selected Data", data)
        self._changedFlag = False

    def setSelectionSource(self, source):
        if self.selectedSource != source:
            self.selectedSource = source
            self.subsetbox.setEnabled(source == OWSelectGenes.SelectInput)
            self.entrybox.setEnabled(source == OWSelectGenes.SelectCustom)
            b = self.selectedSourceButtons.button(source)
            b.setChecked(True)

    def _selectionSourceChanged(self, source):
        if self.selectedSource != source:
            self.selectedSource = source
            self.subsetbox.setEnabled(source == OWSelectGenes.SelectInput)
            self.entrybox.setEnabled(source == OWSelectGenes.SelectCustom)
            self.invalidateOutput()

    def _updateCompletionModel(self):
        var = self.geneVar
        if var is not None:
            names = [str(inst[var]) for inst in self.data
                     if not inst[var].is_special()]
        else:
            names = []

        infodict = {}

        if self.geneinfo[1] is not None:
            info = [(name, self.geneinfo[1].get_info(name, None))
                    for name in names]
            info = filter(itemgetter(1), info)
            infodict = dict(info)

        names = sorted(set(names))
        genes = zip(names, map(infodict.get, names))

        symbols = [info.symbol for _, info in genes if info is not None]

        model = QStandardItemModel()

        def make_row(name, info):
            if info is not None:
                col1 = QStandardItem(name)
                col1.setData(info.symbol, OWSelectGenes.CompletionRole)
                return [col1,
                        QStandardItem(info.symbol),
                        QStandardItem(info.description)]
            else:
                col1 = QStandardItem(name)
                col1.setData(name, OWSelectGenes.CompletionRole)
                return [col1]

        for name, info in genes:
            model.appendRow(make_row(name, info))

        self.geneNames = sorted(set(names) | set(symbols))
        self.genes = genes
        self.entryField.completer().model().setSourceModel(model)
        self.entryField.document().highlighter.setNames(names + symbols)

        self._updatePopupSections()

    def _onGeneIndexChanged(self):
        self._updateCompletionModel()
        self.invalidateOutput()

    def _onSubsetGeneIndexChanged(self):
        if self.selectedSource == OWSelectGenes.SelectInput:
            self.invalidateOutput()

    def _onItemsChanged(self, names):
        selection = set(names).intersection(self.geneNames)
        curr_selection = set(self.selection).intersection(self.geneNames)

        self.selection = names

        if selection != curr_selection:
            self.invalidateOutput()
            to_complete = sorted(set(self.geneNames) - set(names))
            self.entryField.completer().model().setFilterFixedSet(to_complete)

        item = self._selectedSaveSlot()
        if item:
            item.modified = item.savedata != names

    def _onToggleSymbolCompletion(self, state):
        completer = self.entryField.completer()
        completer.setCompletionRole(
            self.CompletionRole if state else Qt.DisplayRole
        )
        self.completeOnSymbols = state
        self._updatePopupSections()

    def _onTranslate(self):
        if self.geneinfo[1] is not None:
            items = self.entryField.items()
            entries = map(self.geneinfo[1].get_info, items)
            items = [info.symbol if info is not None else item
                     for item, info in zip(items, entries)]
            self.entryField.selectAll()
            self.entryField.insertPlainText(" ".join(items))

    def _onGeneInfoReady(self, geneinfo):
        assert QThread.currentThread() is self.thread()
        # Check if the gene info is for the correct (current requested)
        # organism (we might receive a late response from a previous
        # request)
        if self.geneinfo[0] == geneinfo[0]:
            self.geneinfo = geneinfo
            self.translateAction.setEnabled(True)
            self.pb.hide()
            self._updateCompletionModel()

    def _onGeneInfoError(self, exc):
        self.error(0, str(exc))

    def _updatePopupSections(self):
        completer = self.entryField.completer()
        popup = completer.popup()
        assert isinstance(popup, QTreeView)
        header = popup.header()
        # The column in which the symbols should be
        symbol_col = 0 if self.completeOnSymbols else 1
        if symbol_col != header.sectionPosition(1):
            header.moveSection(0, 1)

    def _selectedSaveSlot(self):
        """
        Return the current selected saved selection slot.
        """
        indexes = self.selectionsView.selectedIndexes()
        if indexes:
            return self.selectionsModel.item(indexes[0].row())
        else:
            return None

    def saveSelection(self):
        """
        Save (update) the items in the current selected selection.
        """
        item = self._selectedSaveSlot()
        if item:
            item.savedata = self.entryField.items()
            item.modified = False

    def copyToSaved(self):
        """
        Copy the current 'Gene Subset' names to saved selections.
        """
        if self.subsetGeneVar and \
                self.selectedSource == OWSelectGenes.SelectInput:
            names = self.selectedGenes()
            item = SaveSlot("New selection")
            item.savedata = names
            self.selectionsModel.appendRow([item])
            self.setSelectionSource(OWSelectGenes.SelectCustom)
            self.selectionsView.setCurrentIndex(item.index())
            self.selectionsView.edit(item.index())

    def appendToSelected(self):
        """
        Append the current 'Gene Subset' names to 'Select Genes' entry field.
        """
        if self.subsetGeneVar and \
                self.selectedSource == OWSelectGenes.SelectInput:
            names = self.selectedGenes()
            text = " ".join(names)
            self.entryField.appendPlainText(text)
            self.setSelectionSource(OWSelectGenes.SelectCustom)
            self.entryField.setFocus()
            self.entryField.moveCursor(QTextCursor.End)

    def addSelection(self, name=None):
        """
        Add a new saved selection entry initialized by the current items.

        The new slot will be selected.

        """
        item = SaveSlot(name or "New selection")
        item.savedata = self.entryField.items()
        self.selectionsModel.appendRow([item])
        self.selectionsView.setCurrentIndex(item.index())

        if not name:
            self.selectionsView.edit(item.index())

    def removeSelection(self):
        """
        Remove the current selected save slot.
        """
        item = self._selectedSaveSlot()
        if item:
            self.selectionsModel.removeRow(item.row())

    def importGeneSet(self):
        if self._genesetDialog is None:
            self._genesetDialog = GeneSetDialog(
                self, windowTitle="Import Gene Set Names")

        dialog = self._genesetDialog

        if self.taxid is not None:
            dialog.setCurrentOrganism(self.taxid)

        result = dialog.exec_()
        if result == QDialog.Accepted:
            gsets = dialog.selectedGeneSets()
            genes = reduce(operator.ior, (gs.genes for gs in gsets), set())
            text = " ".join(genes)
            self.entryField.appendPlainText(text)
            self.entryField.setFocus()
            self.entryField.moveCursor(QTextCursor.End)

    def importFromFile(self):
        filename = QFileDialog.getOpenFileName(
            self, "Open File", os.path.expanduser("~/"))

        if filename:
            filename = unicode(filename)
            with open(filename, "rU") as f:
                text = f.read()
            self.entryField.appendPlainText(text)
            self.entryField.setFocus()
            self.entryField.moveCursor(QTextCursor.End)

    def _onSelectedSaveSlotChanged(self):
        item = self._selectedSaveSlot()
        if item:
            if not item.document:
                item.document = self._createDocument()
                if item.savedata:
                    item.document.setPlainText(" ".join(item.savedata))

            item.document.highlighter.setNames(self.geneNames)

            self.entryField.setDocument(item.document)

        self._updateActions()

    def _createDocument(self):
        """
        Create and new QTextDocument instance for editing gene names.
        """
        doc = QTextDocument(self)
        doc.setDocumentLayout(QPlainTextDocumentLayout(doc))
        doc.highlighter = NameHighlight(doc)
        return doc

    def _updateActions(self):
        """
        Update the Save/remove action enabled state.
        """
        selected = bool(self._selectedSaveSlot())
        self.actionRemove.setEnabled(selected)
        self.actionSave.setEnabled(selected)

    def getSettings(self, *args, **kwargs):
        # copy the saved selections model back to widget settings.
        selections = []
        for i in range(self.selectionsModel.rowCount()):
            item = self.selectionsModel.item(i)
            selections.append((item.name, item.savedata))
        self.savedSelections = selections

        item = self._selectedSaveSlot()
        if item is None:
            self.selectedSelectionIndex = -1
        else:
            self.selectedSelectionIndex = item.row()

        return OWWidget.getSettings(self, *args, **kwargs)

    def sendReport(self):
        report = []
        if self.data is not None:
            report.append("%i instances on input." % len(self.data))
        else:
            report.append("No data on input.")

        if self.geneVar is not None:
            report.append("Gene names taken from %r attribute." %
                          escape(self.geneVar.name))

        self.reportSection("Input")
        self.startReportList()
        for item in report:
            self.addToReportList(item)
        self.finishReportList()
        report = []
        selection = self.selectedGenes()
        if self.selectedSource == OWSelectGenes.SelectInput:
            self.reportRaw(
                "<p>Gene Selection (from 'Gene Subset' input): %s</p>" %
                escape(" ".join(selection))
            )
        else:
            self.reportRaw(
                "<p>Gene Selection: %s</p>" %
                escape(" ".join(selection))
            )
        self.reportSettings(
            "Settings",
            [("Preserve order", self.preserveOrder)]
        )

    def onDeleteWidget(self):
        self._inittask.future().cancel()

        if self._infotask:
            self._infotask.future().cancel()

        self._executor.shutdown(wait=True)
        OWWidget.onDeleteWidget(self)


def is_string(feature):
    return isinstance(feature, Orange.feature.String)


def domain_variables(domain):
    """
    Return all feature descriptors from the domain.
    """
    vars = (domain.features +
            domain.class_vars +
            domain.getmetas().values())
    return vars


def gene_candidates(data):
    """
    Return features that could contain gene names.
    """
    vars = domain_variables(data.domain)
    vars = filter(is_string, vars)
    return vars


def select_by_genes(data, gene_feature, gene_list, preserve_order=True):
    if preserve_order:
        selection = set(gene_list)
        sel = [inst for inst in data
               if str(inst[gene_feature]) in selection]
    else:
        by_genes = defaultdict(list)
        for inst in data:
            by_genes[str(inst[gene_feature])].append(inst)

        sel = []
        for name in gene_list:
            sel.extend(by_genes.get(name, []))

    if sel:
        data = Orange.data.Table(data.domain, sel)
    else:
        data = Orange.data.Table(data.domain)

    return data


_CompletionState = namedtuple(
    "_CompletionState",
    ["start",  # completion prefix start position
     "pos",  # cursor position
     "anchor"]  # anchor position (inline completion end)
)


class ListTextEdit(QPlainTextEdit):
    """
    A text editor specialized for editing a list of items.
    """
    #: Emitted when the list items change.
    itemsChanged = Signal(list)

    def __init__(self, parent=None, **kwargs):
        QPlainTextEdit.__init__(self, parent, **kwargs)

        self._items = None
        self._completer = None
        self._completionState = _CompletionState(-1, -1, -1)

        self.cursorPositionChanged.connect(self._cursorPositionChanged)
        self.textChanged.connect(self._textChanged)

    def setCompleter(self, completer):
        """
        Set a completer for list items.
        """
        if self._completer is not None:
            self._completer.setWidget(None)
            self._completer.activated.disconnect(self._insertCompletion)

        self._completer = completer

        if self._completer:
            self._completer.setWidget(self)
            self._completer.activated.connect(self._insertCompletion)

    def completer(self):
        """
        Return the completer.
        """
        return self._completer

    def setItems(self, items):
        text = " ".join(items)
        self.setPlainText(text)

    def items(self):
        if self._items is None:
            self._items = self._getItems()
        return self._items

    def keyPressEvent(self, event):
        # TODO: in Qt 4.8 QPlainTextEdit uses inputMethodEvent for
        # non-ascii input

        if self._completer.popup().isVisible():
            if event.key() in [Qt.Key_Enter, Qt.Key_Return, Qt.Key_Escape,
                               Qt.Key_Tab, Qt.Key_Backtab]:
                # These need to propagate to the completer.
                event.ignore()
                return

        QPlainTextEdit.keyPressEvent(self, event)

        if not len(event.text()) or not is_printable(unicode(event.text())[0]):
            return

        text = unicode(self.toPlainText())
        cursor = self.textCursor()
        pos = cursor.position()

        if pos == len(text) or not text[pos].strip():
            # cursor is at end of text or whitespace
            # find the beginning of the current word
            whitespace = " \t\n\r\f\v"
            start = max([text.rfind(c, 0, pos) for c in whitespace]) + 1

            prefix = text[start:pos]

            if prefix:
                if self._completer.completionPrefix() != prefix:
                    self._completer.setCompletionPrefix(text[start:pos])

                rect = self.cursorRect()
                popup = self._completer.popup()
                if popup.isVisible():
                    rect.setWidth(popup.width())
                else:
                    view_adjust_size_to_contents(popup)
                    rect.setWidth(popup.width())

                # Popup the completer list
                self._completer.complete(rect)

                # Inline completion of a common prefix
                inline = self._commonCompletionPrefix()
                inline = inline[len(prefix):]

                self._completionState = \
                    _CompletionState(start, pos, pos + len(inline))

                cursor.insertText(inline)
                cursor.setPosition(pos, QTextCursor.KeepAnchor)
                self.setTextCursor(cursor)

            elif self._completer.popup().isVisible():
                self._stopCompletion()

    def _cursorPositionChanged(self):
        cursor = self.textCursor()
        pos = cursor.position()
        start, _, _ = self._completionState

        if start == -1:
            # completion not in progress
            return

        if pos <= start:
            # cursor moved before the start of the prefix
            self._stopCompletion()
            return

        text = unicode(self.toPlainText())
        # Find the end of the word started by completion prefix
        word_end = len(text)
        for i in range(start, len(text)):
            if text[i] in " \t\n\r\f\v":
                word_end = i
                break

        if pos > word_end:
            # cursor moved past the word boundary
            self._stopCompletion()

        # TODO: Update the prefix when moving the cursor
        # inside the word

    def _insertCompletion(self, item):
        if isinstance(item, list):
            completion = " ".join(item)
        else:
            completion = unicode(item)

        start, _, end = self._completionState

        self._stopCompletion()

        cursor = self.textCursor()
        # Replace the prefix+inline with the full completion
        # (correcting for the case-insensitive search).
        cursor.setPosition(min(end, self.document().characterCount()))
        cursor.setPosition(start, QTextCursor.KeepAnchor)

        cursor.insertText(completion + " ")

    def _commonCompletionPrefix(self):
        """
        Return the common prefix of items in the current completion model.
        """
        model = self._completer.completionModel()
        column = self._completer.completionColumn()
        role = self._completer.completionRole()
        items = [toString(model.index(i, column).data(role))
                 for i in range(model.rowCount())]

        if not items:
            return ""

        first = min(items)
        last = max(items)
        for i, c in enumerate(first):
            if c != last[i]:
                return first[:i]

        return first

    def _stopCompletion(self):
        self._completionState = _CompletionState(-1, -1, -1)
        if self._completer.popup().isVisible():
            self._completer.popup().hide()

    def _textChanged(self):
        items = self._getItems()
        if self._items != items:
            self._items = items
            self.itemsChanged.emit(items)

    def _getItems(self):
        """
        Return the current items (a list of strings).

        .. note:: The inline completion text is not included.

        """
        text = unicode(self.toPlainText())
        if self._completionState[0] != -1:
            # Remove the inline completion text
            _, pos, end = self._completionState
            text = text[:pos] + text[end:]
        return [item for item in text.split() if item.strip()]


def view_adjust_column_sizes(view, maxWidth=None):
    """
    Adjust view's column sizes to to contents.
    """
    if maxWidth is None:
        maxWidth = sys.maxint

    for col in range(view.model().columnCount()):
        width = min(view.sizeHintForColumn(col), maxWidth)
        view.setColumnWidth(col, width)


def view_adjust_size_to_contents(view):
    """
    Adjust the view to a reasonable size based in it's contents.
    """
    view_adjust_column_sizes(view)
    w = sum([view.columnWidth(col)
             for col in range(view.model().columnCount())])
    w += view.verticalScrollBar().sizeHint().width()

    h = view.sizeHintForRow(0) * 7
    h += view.horizontalScrollBar().sizeHint().height()
    view.resize(w, h)


class NameHighlight(QSyntaxHighlighter):
    def __init__(self, parent=None, **kwargs):
        super(NameHighlight, self).__init__(parent, **kwargs)

        self._names = set()

        self._format = QTextCharFormat()
        self._format.setForeground(Qt.blue)

        self._unrecognized_format = QTextCharFormat()
#         self._unrecognized_format.setFontStrikeOut(True)

    def setNames(self, names):
        self._names = set(names)
        self.rehighlight()

    def names(self):
        return set(self._names)

    def highlightBlock(self, text):
        text = unicode(text)
        pattern = re.compile(r"\S+")
        for match in pattern.finditer(text):
            name = text[match.start(): match.end()]
            match_len = match.end() - match.start()

            if not name.strip():
                continue

            if name in self._names:
                format = self._format
            else:
                format = self._unrecognized_format

            self.setFormat(match.start(), match_len, format)


@contextmanager
def signals_blocked(obj):
    blocked = obj.signalsBlocked()
    obj.blockSignals(True)
    try:
        yield
    finally:
        obj.blockSignals(blocked)


class ListCompleter(QCompleter):
    """
    A completer supporting selection of multiple list items.
    """
    activated = Signal(list)

    def __init__(self, *args, **kwargs):
        QCompleter.__init__(self, *args, **kwargs)

        popup = QListView()
        popup.setEditTriggers(QListView.NoEditTriggers)
        popup.setSelectionMode(QListView.ExtendedSelection)

        self.setPopup(popup)

    def setPopup(self, popup):
        QCompleter.setPopup(self, popup)

        popup.viewport().installEventFilter(self)
        popup.doubleClicked.connect(self._complete)

    def eventFilter(self, receiver, event):
        if event.type() == QEvent.KeyPress and receiver is self.popup():
            if event.key() in [Qt.Key_Enter, Qt.Key_Return, Qt.Key_Tab]:
                self._complete()
                return True

        elif event.type() == QEvent.MouseButtonRelease and \
                receiver is self.popup().viewport():
            # Process the event without emitting 'clicked', ... signal to
            # override the default QCompleter behavior
            with signals_blocked(self.popup()):
                QApplication.sendEvent(self.popup(), event)
                return True

        return QCompleter.eventFilter(self, receiver, event)

    def _complete(self):
        column = self.completionColumn()
        role = self.completionRole()
        indexes = self.popup().selectionModel().selectedRows(column)

        items = [toString(index.data(role)) for index in indexes]

        if self.popup().isVisible():
            self.popup().hide()

        if items:
            self.activated.emit(items)


class SetFilterProxyModel(QSortFilterProxyModel):
    def __init__(self, *args, **kwargs):
        super(SetFilterProxyModel, self).__init__(*args, **kwargs)
        self._filterFixedSet = None

    def setFilterFixedSet(self, items):
        if items is None:
            self._filterFixedSet = None
        else:
            self._filterFixedSet = set(items)
        self.invalidateFilter()

    def filterAcceptsRow(self, row, parent):
        if self._filterFixedSet is None:
            return True

        model = self.sourceModel()
        col = self.filterKeyColumn()
        var = model.data(model.index(row, col, parent),
                         self.filterRole())
        var = toString(var)
        return var in self._filterFixedSet


class ComboBox(QComboBox):
    """
    A combo box ignoring mouse wheel events.
    """
    def wheelEvent(self, event):
        event.ignore()
        return

# All control character categories.
_control = set(["Cc", "Cf", "Cs", "Co", "Cn"])


def is_printable(unichar):
    """
    Return True if the unicode character `unichar` is a printable character.
    """
    return unicodedata.category(unichar) not in _control


import sys
import itertools

from PyQt4.QtGui import (
    QVBoxLayout, QLineEdit, QDialogButtonBox,
    QProgressBar, QSizePolicy
)

from PyQt4.QtCore import QSize
from Orange.bio import obiGeneSets as genesets
from Orange.utils import serverfiles


class GeneSetView(QFrame):
    selectedOrganismChanged = Signal(str)
    selectionChanged = Signal()
    geneSetsLoaded = Signal()

    def __init__(self, *args, **kwargs):
        super(GeneSetView, self).__init__(*args, **kwargs)

        self._taxid = None

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        self._stack = QStackedWidget()
        self._stack.setContentsMargins(0, 0, 0, 0)
        self._stack.setSizePolicy(QSizePolicy.MinimumExpanding,
                                  QSizePolicy.Fixed)
        self.orgcombo = ComboBox(minimumWidth=150, focusPolicy=Qt.StrongFocus)
        self.orgcombo.activated[int].connect(self._on_organismSelected)
        self._stack.addWidget(self.orgcombo)

        self.progressbar = QProgressBar()
        self._stack.addWidget(self.progressbar)

        layout.addWidget(self._stack)

        self.searchline = QLineEdit()
        self.searchline.setPlaceholderText("Filter...")
        completer = QCompleter()
        self.searchline.setCompleter(completer)
        layout.addWidget(self.searchline)

        self.gsview = QTreeView()
        self.gsview.setAlternatingRowColors(True)
        self.gsview.setRootIsDecorated(False)
        self.gsview.setSelectionMode(QTreeView.ExtendedSelection)
        self.gsview.setSortingEnabled(True)
        self.gsview.setUniformRowHeights(True)
        self.proxymodel = QSortFilterProxyModel(
            filterKeyColumn=1, sortCaseSensitivity=Qt.CaseInsensitive
        )

        self.gsview.setModel(self.proxymodel)
        self.gsview.selectionModel().selectionChanged.connect(
            self._on_selectionChanged)

        self.searchline.textChanged.connect(
            self.proxymodel.setFilterFixedString)

        layout.addWidget(self.gsview)
        self.setLayout(layout)

        self._executor = ThreadExecutor(self)
        self.initialize()

    def initialize(self):
        self.gs_hierarchy = gs = genesets.list_all()
        taxids = set(taxid for _, taxid, _ in gs)
        self.organisms = [(taxid, taxonomy.name(taxid)) for taxid in taxids]
        for taxid, name in self.organisms:
            self.orgcombo.addItem(name, taxid)

        self.orgcombo.setCurrentIndex(-1)

    def sizeHint(self):
        return QSize(500, 550)

    def setCurrentOrganism(self, taxid):
        taxids = [tid for tid, _ in self.organisms]
        if taxid is not None and taxid not in taxids:
            taxid = None

        if taxid != self._taxid:
            self._taxid = taxid
            if taxid is None:
                self.orgcombo.setCurrentIndex(-1)
            else:
                index = taxids.index(taxid)
                self.orgcombo.setCurrentIndex(index)
            self._updateGeneSetsModel()
            self.selectedOrganismChanged.emit(taxid)

    def currentOrganism(self):
        return self._taxid

    def selectedGeneSets(self):
        selmod = self.gsview.selectionModel()
        model = self.proxymodel.sourceModel()
        rows = [self.proxymodel.mapToSource(row)
                for row in selmod.selectedRows(1)]
        gsets = [model.data(row, Qt.UserRole) for row in rows]
        return map(toPyObject, gsets)

    def _updateGeneSetsModel(self):
        taxid = self._taxid
        if taxid is None:
            self.proxymodel.setSourceModel(None)
        else:
            currentsets = [(hier, tid)
                           for hier, tid, _ in self.gs_hierarchy
                           if tid == taxid]

            gsmissing = [(hier, tid, local)
                         for hier, tid, local in self.gs_hierarchy
                         if tid == taxid and not local]

            self._stack.setCurrentWidget(self.progressbar)

            if gsmissing:
                self.progressbar.setRange(0, 100)
                progress_info = methodinvoke(
                    self.progressbar, "setValue", (int,))
            else:
                self.progressbar.setRange(0, 0)
                progress_info = None

            def load():
                gs_ensure_downloaded(
                    gsmissing,
                    progress_info=progress_info)

                return [((hier, tid), genesets.load(hier, tid))
                        for hier, tid in currentsets]

            self._task = Task(function=load)
            self._task.finished.connect(self._on_loadFinished)
            self._executor.submit(self._task)

    def _on_loadFinished(self):
        assert QThread.currentThread() is self.thread()
        self._stack.setCurrentWidget(self.orgcombo)

        try:
            sets = self._task.result()
        except Exception:
            # Should do something better here.
            sys.excepthook(*sys.exc_info())
            sets = []

        model = sets_to_model(sets)
        self.proxymodel.setSourceModel(model)
        self.gsview.resizeColumnToContents(0)
        self.geneSetsLoaded.emit()

    def _on_organismSelected(self, index):
        if index != -1:
            item = self.orgcombo.model().item(index)
            taxid = toString(item.data(Qt.UserRole))
            self.setCurrentOrganism(taxid)

    def _on_selectionChanged(self, *args):
        self.selectionChanged.emit()


def sets_to_model(gsets):
    model = QStandardItemModel()
    model.setHorizontalHeaderLabels(["Category", "Name", "#"])

    for (hier, tid), sets in gsets:
        for gset in sets:
            ngenes = len(gset.genes)
            names = [escape(name) for name in list(gset.genes)[:30]]
            names = ", ".join(names)
            tooltip = "<p>{0}</p>{1}".format(escape(gset.name), names)
            if ngenes > 30:
                tooltip += ", ... ({0} names not shown)".format(ngenes - 30)

            category = QStandardItem(" ".join(hier))
            category.setData((hier, tid), Qt.UserRole)
            category.setEditable(False)
            category.setToolTip(tooltip)
            name = QStandardItem(gset.name)
            name.setData(gset, Qt.UserRole)
            name.setEditable(False)
            name.setToolTip(tooltip)
            count = QStandardItem()
            count.setData(ngenes, Qt.DisplayRole)
            count.setToolTip(tooltip)
            model.appendRow([category, name, count])

    return model


def gs_ensure_downloaded(gslist, progress_info=None):
    hierlist = [(hier, taxid) for hier, taxid, local in gslist
                if not local]

    files = [(genesets.sfdomain, genesets.filename(hier, taxid))
             for hier, taxid in hierlist]

    download_list(files, progress_info)


def download_list(files_list, progress_callback=None):
    nfiles = len(files_list)
    count = 100 * nfiles
    counter = itertools.count()

    def advance():
        progress_callback(100.0 * next(counter) / count)

    for domain, filename in files_list:
        serverfiles.download(domain, filename,
                             callback=advance if progress_callback else None)


class GeneSetDialog(QDialog):
    selectionChanged = Signal()

    def __init__(self, *args, **kwargs):
        super(GeneSetDialog, self).__init__(*args, **kwargs)
        layout = QVBoxLayout()
        layout.setContentsMargins(4, 4, 4, 4)
        self.gsview = GeneSetView()
        self.gsview.selectionChanged.connect(self.selectionChanged)

        layout.addWidget(self.gsview)

        buttonbox = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal
        )
        buttonbox.accepted.connect(self.accept)
        buttonbox.rejected.connect(self.reject)

        layout.addWidget(buttonbox)

        self.setLayout(layout)

    def selectedGeneSets(self):
        return self.gsview.selectedGeneSets()

    def setCurrentOrganism(self, taxid):
        self.gsview.setCurrentOrganism(taxid)

    def currentOrganism(self):
        return self.gsview.currentOrganism()


def test1():
    app = QApplication([])
    dlg = GeneSetDialog()
    dlg.exec_()
    del dlg
    app.processEvents()


def test():
    app = QApplication([])
    w = OWSelectGenes()
    data = Orange.data.Table("brown-selected")
    w.setData(data)
    w.setGeneSubset(Orange.data.Table(data[:10]))
    w.show()
    app.exec_()
    w.saveSettings()
    w.deleteLater()
    del w
    app.processEvents()

if __name__ == "__main__":
    test()
