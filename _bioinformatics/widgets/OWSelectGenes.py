import re
import unicodedata

from PyQt4.QtGui import (
    QLabel, QWidget, QPlainTextEdit, QSyntaxHighlighter, QTextCharFormat,
    QTextCursor, QCompleter, QStringListModel, QListView
)

from PyQt4.QtCore import Qt, pyqtSignal as Signal

import Orange

from Orange.OrangeWidgets.OWWidget import *
from Orange.OrangeWidgets.OWItemModels import VariableListModel
from Orange.OrangeWidgets import OWGUI


NAME = "Select Genes"
DESCRIPTION = "Select a specified subset of the input genes."
ICON = "icons/SelectGenes.svg"

INPUTS = [("Data", Orange.data.Table, "set_data")]
OUTPUTS = [("Selected Data", Orange.data.Table)]


class OWSelectGenes(OWWidget):

    contextHandlers = {
        "": DomainContextHandler(
            "", ["geneIndex", "selection"]
        )
    }

    settingsList = ["autoCommit"]

    def __init__(self, parent=None, signalManager=None, title=NAME):
        OWWidget.__init__(self, parent, signalManager, title,
                          wantMainArea=False)

        self.selection = []
        self.geneIndex = None
        self.autoCommit = False

        self.loadSettings()

        # Input variables that could contain names
        self.variables = VariableListModel()
        # All gene names from the input (in self.geneIndex column)
        self.geneNames = []
        # Output changed flag
        self._changedFlag = False

        box = OWGUI.widgetBox(self.controlArea, "Gene Attribute")
        self.attrsCombo = OWGUI.comboBox(
            box, self, "geneIndex",
            callback=self._onGeneIndexChanged,
            tooltip="Column with gene names"
        )
        self.attrsCombo.setModel(self.variables)

        box = OWGUI.widgetBox(self.controlArea, "Gene Selection")
        self.entryField = ListTextEdit(box)
        self.entryField.setTabChangesFocus(True)
        self.entryField.setToolTip("Enter selected gene names")
        self.entryField.textChanged.connect(self._textChanged)

        box.layout().addWidget(self.entryField)

        completer = QCompleter(self)
        completer.setCompletionMode(QCompleter.PopupCompletion)
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        completer.popup().setAlternatingRowColors(True)
        completer.setModel(QStringListModel([], self))

        self.entryField.setCompleter(completer)

        self.hightlighter = NameHighlight(self.entryField.document())

        box = OWGUI.widgetBox(self.controlArea, "Output")

        cb = OWGUI.checkBox(box, self, "autoCommit", "Auto commit")
        button = OWGUI.button(box, self, "Commit", callback=self.commit)

        OWGUI.setStopper(self, button, cb, "_changedFlag", self.commit)

    def set_data(self, data):
        """
        Set the input data.
        """
        self.closeContext("")
        self.warning()
        self.data = data
        if data is not None:
            attrs = gene_candidates(data)
            self.variables[:] = attrs
            self.attrsCombo.setCurrentIndex(0)
            self.geneIndex = 0
            self.selection = []
        else:
            self.variables[:] = []
            self.geneIndex = -1
            self.warning(0, "No suitable columns for gene names.")

        self._changedFlag = True
        self._updateCompletionModel()

        self.openContext("", data)

        self.entryField.setPlainText(" ".join(self.selection))

        self.commit()

    @property
    def geneVar(self):
        if self.data is not None and self.geneIndex >= 0:
            return self.variables[self.geneIndex]
        else:
            return None

    def invalidateOutput(self):
        if self.autoCommit:
            self.commit()
        else:
            self._changedFlag = True

    def commit(self):
        gene = self.geneVar

        if gene is not None:
            selection = set(self.selection)

            sel = [inst for inst in self.data
                   if str(inst[gene]) in selection]

            if sel:
                data = Orange.data.Table(self.data.domain, sel)
            else:
                data = Orange.data.Table(self.data.domain)

        else:
            data = None

        self.send("Selected Data", data)
        self._changedFlag = False

    def _updateCompletionModel(self):
        var = self.geneVar
        if var is not None:
            names = [str(inst[var]) for inst in self.data
                     if not inst[var].isSpecial()]
        else:
            names = []

        self.geneNames = names
        self.entryField.completer().model().setStringList(sorted(set(names)))
        self.hightlighter.setNames(names)

    def _onGeneIndexChanged(self):
        self._updateCompletionModel()
        self.invalidateOutput()

    def _textChanged(self):
        names = self.entryField.list()
        selection = set(names).intersection(self.geneNames)
        curr_selection = set(self.selection).intersection(self.geneNames)

        if selection != curr_selection:
            self.selection = names
            self.invalidateOutput()

            names = set(self.geneNames) - set(names)
            self.entryField.completer().model().setStringList(sorted(names))


def is_string(feature):
    return isinstance(feature, Orange.feature.String)


def domain_variables(domain):
    vars = (domain.features +
            domain.class_vars +
            domain.getmetas().values())
    return vars


def gene_candidates(data):
    vars = domain_variables(data.domain)
    vars = filter(is_string, vars)
    return vars


class ListTextEdit(QPlainTextEdit):
    listChanged = Signal()

    def __init__(self, parent=None, **kwargs):
        QPlainTextEdit.__init__(self, parent, **kwargs)

        self._completer = None

    def setCompleter(self, completer):
        """
        Set a completer for list items.
        """
        if self._completer is not None:
            self._completer.setWidget(None)
            if isinstance(self._completer, ListCompleter):
                self._completer.activatedList.disconnect(self._insertCompletion)
            else:
                self._completer.activated.disconnect(self._insertCompletion)

        self._completer = completer

        if self._completer:
            self._completer.setWidget(self)
            if isinstance(self._completer, ListCompleter):
                self._completer.activatedList.connect(self._insertCompletion)
            else:
                self._completer.activated.connect(self._insertCompletion)

    def completer(self):
        """
        Return the completer.
        """
        return self._completer

    def setList(self, list):
        text = " ".join(list)
        self.setPlainText(text)

    def list(self):
        return [name for name in unicode(self.toPlainText()).split()
                if name.strip()]

    def keyPressEvent(self, event):
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
        pos = self.textCursor().position()

        if pos == len(text) or not(text[pos].strip()):
            # At end of text or whitespace
            # TODO: Match all whitespace characters.
            start_sp = text.rfind(" ", 0, pos) + 1
            start_n = text.rfind("\n", 0, pos) + 1
            start = max(start_sp, start_n)

            prefix = text[start:pos]

            if prefix:
                if self._completer.completionPrefix() != prefix:
                    self._completer.setCompletionPrefix(text[start:pos])

                rect = self.cursorRect()
                popup = self._completer.popup()
                if popup.isVisible():
                    rect.setWidth(popup.width())
                else:
                    rect.setWidth(popup.sizeHintForColumn(0) +
                                  popup.verticalScrollBar().sizeHint().width())

                self._completer.complete(rect)

            elif self._completer.popup().isVisible():
                self._completer.popup().hide()

    def _insertCompletion(self, item):
        completion = unicode(item)
        prefix = self._completer.completionPrefix()

        cursor = self.textCursor()
        # Replace the prefix with the full completion (correcting for the
        # case-insensitive search).
        cursor.setPosition(cursor.position() - len(prefix),
                           QTextCursor.KeepAnchor)

        cursor.insertText(completion + " ")


class NameHighlight(QSyntaxHighlighter):
    def __init__(self, parent=None, **kwargs):
        super(NameHighlight, self).__init__(parent)

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


def toString(variant):
    if isinstance(variant, QVariant):
        return unicode(variant.toString())
    else:
        return unicode(variant)


class ListCompleter(QCompleter):
    activatedList = Signal(list)

    def __init__(self, *args, **kwargs):
        QCompleter.__init__(self, *args, **kwargs)

        popup = QListView()
        popup.setSelectionMode(QListView.ExtendedSelection)
        self.setPopup(popup)

    def setPopup(self, popup):
        QCompleter.setPopup(self, popup)

        popup.selectionModel().selectionChanged.connect(
            self._completionSelected)

    def _completionSelected(self, selected, deselected):
        selection = self.popup().selectionModel().selection()
        indexes = selection.indexes()

        items = [toString(index.data(self.completionRole()))
                 for index in indexes]

        self.activatedList.emit(items)


# All control character categories.
_control = set(["Cc", "Cf", "Cs", "Co", "Cn"])


def is_printable(unichar):
    """
    Return True if the unicode character `unichar` is a printable character.
    """
    return unicodedata.category(unichar) not in _control


def test():
    app = QApplication([])
    w = OWSelectGenes()
    data = Orange.data.Table("brown-selected")
    w.set_data(data)
    w.show()
    app.exec_()
    w.deleteLater()
    del w
    app.processEvents()

if __name__ == "__main__":
    test()
