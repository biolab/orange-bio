import re
import unicodedata
from collections import defaultdict, namedtuple

from contextlib import contextmanager

from PyQt4.QtGui import (
    QLabel, QWidget, QPlainTextEdit, QSyntaxHighlighter, QTextCharFormat,
    QTextCursor, QCompleter, QStringListModel, QListView
)

from PyQt4.QtCore import Qt, QEvent, pyqtSignal as Signal

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

    settingsList = ["autoCommit", "preserveOrder"]

    def __init__(self, parent=None, signalManager=None, title=NAME):
        OWWidget.__init__(self, parent, signalManager, title,
                          wantMainArea=False)

        self.selection = []
        self.geneIndex = None
        self.autoCommit = False
        self.preserveOrder = True

        self.loadSettings()

        # Input variables that could contain names
        self.variables = VariableListModel()
        # All gene names from the input (in self.geneIndex column)
        self.geneNames = []
        # Output changed flag
        self._changedFlag = False
        self.data = None

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
        self.entryField.itemsChanged.connect(self._itemsChanged)

        box.layout().addWidget(self.entryField)

        completer = ListCompleter(self)
        completer.setCompletionMode(QCompleter.PopupCompletion)
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        completer.setMaxVisibleItems(10)
        completer.popup().setAlternatingRowColors(True)
        completer.setModel(QStringListModel([], self))

        self.entryField.setCompleter(completer)

        self.hightlighter = NameHighlight(self.entryField.document())

        box = OWGUI.widgetBox(self.controlArea, "Output")
        OWGUI.checkBox(box, self, "preserveOrder", "Preserve input order",
                       tooltip="Preserve the order of the input data "
                               "instances.",
                       callback=self.invalidateOutput)
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
            if attrs:
                self.geneIndex = 0
            else:
                self.geneIndex = -1
                self.warning(0, "No suitable columns for gene names.")

            self.selection = []
        else:
            self.variables[:] = []
            self.geneIndex = -1

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
            if self.preserveOrder:
                selection = set(self.selection)
                sel = [inst for inst in self.data
                       if str(inst[gene]) in selection]
            else:
                by_genes = defaultdict(list)
                for inst in self.data:
                    by_genes[str(inst[gene])].append(inst)

                sel = []
                for name in self.selection:
                    sel.extend(by_genes.get(name, []))

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

    def _itemsChanged(self, names):
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


_CompletionState = namedtuple(
    "_CompletionState",
    ["start",  # completion prefix start position
     "pos",  # cursor position
     "anchor"]  # anchor position (inline completion end)
)


class ListTextEdit(QPlainTextEdit):
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

        if pos == len(text) or not(text[pos].strip()):
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
                    rect.setWidth(popup.sizeHintForColumn(0) +
                                  popup.verticalScrollBar().sizeHint().width())

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
        text = unicode(self.toPlainText())
        if self._completionState[0] != -1:
            # Remove the inline completion text from the text
            _, pos, end = self._completionState
            text = text[:pos] + text[end:]
        return [item for item in text.split() if item.strip()]


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


def toString(variant):
    if isinstance(variant, QVariant):
        return unicode(variant.toString())
    else:
        return unicode(variant)


@contextmanager
def signals_blocked(obj):
    blocked = obj.signalsBlocked()
    obj.blockSignals(True)
    try:
        yield
    finally:
        obj.blockSignals(blocked)


class ListCompleter(QCompleter):
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
        selection = self.popup().selectionModel().selection()
        indexes = selection.indexes()

        items = [toString(index.data(self.completionRole()))
                 for index in indexes]

        if self.popup().isVisible():
            self.popup().hide()

        if items:
            self.activated.emit(items)


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
