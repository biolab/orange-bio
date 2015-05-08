
import sys
import os
import time

from collections import defaultdict

from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QTreeWidget, QTreeWidgetItem

import Orange.data

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.datacaching import data_hints

from ..utils import environ
from .. import dicty


# set buffer file
bufferpath = os.path.join(environ.buffer_dir, "dicty")
try:
    os.makedirs(bufferpath)
except OSError:
    pass

bufferfile = os.path.join(bufferpath, "database.sq3")


class OWDicty(widget.OWWidget):
    name = "dictyExpress"
    description = "Access to data in dictyExpress database."
    icon = "../widgets/icons/DictyExpress.svg"
    priority = 40

    inputs = []
    outputs = [("Data", Orange.data.Table)]

    serverToken = settings.Setting("")
    platform = settings.Setting(None)
    selectedExperiemtns = settings.Setting([])
    server = settings.Setting(dicty.defaddress)
    excludeconstant = settings.Setting(False)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.buffer = dicty.CacheSQLite(bufferfile)

        self.dbc = None
        self.chipsl = []
        self.items = []
        self.searchString = ""

        box = gui.widgetBox(self.controlArea, "Cache")
        gui.button(box, self, "Clear cache", callback=self.clear_buffer)

        gui.checkBox(self.controlArea, self, "excludeconstant",
                     "Exclude labels with constant values")

        gui.button(self.controlArea, self, "&Commit", callback=self.Commit)
        box = gui.widgetBox(self.controlArea, "Server")

        gui.lineEdit(box, self, "serverToken", "Token", callback=self.Connect)
        gui.rubber(self.controlArea)

        gui.lineEdit(self.mainArea, self, "searchString", "Search",
                     callbackOnType=True, callback=self.SearchUpdate)

        self.experimentsWidget = QTreeWidget(
            selectionMode=QTreeWidget.ExtendedSelection,
            rootIsDecorated=True,
            sortingEnabled=True)

        self.experimentsWidget.setHeaderLabels(
            ["Strain", "Treatment", "Growth condition",
             "Platform", "N", "Chips"])

        self.mainArea.layout().addWidget(self.experimentsWidget)

        QtCore.QTimer.singleShot(0, self.UpdateExperiments)

    def sizeHint(self):
        return QtCore.QSize(800, 600)

    def __updateSelectionList(self, oldList, oldSelection, newList):
        oldList = [oldList[i] for i in oldSelection]
        return [i for i, new in enumerate(newList) if new in oldList]

    def Connect(self):
        address = self.server
        if not address.endswith("?"):
            address = address + "?"

        if self.serverToken:
            address += "token=" + self.serverToken + "&"
        try:
            self.dbc = dicty.DatabaseConnection(address, cache=self.buffer)
        except Exception as ex:
            sys.excepthook(*sys.exc_info())
            self.error(0, "Error connecting to server" + str(ex))
            return
        self.error(0)

    def clear_buffer(self):
        self.buffer.clear()
        self.UpdateExperiments()

    def UpdateExperiments(self):
        self.chipsl = []
        self.experimentsWidget.clear()
        self.items = []

        if not self.dbc:
            self.Connect()

        annotations = self.dbc.annotations("norms", None)

        elements = []

        for chip, annot in annotations:
            d = dict(annot)
            elements.append(((d.get("treatment", ""),
                              d.get("growthCond", ""),
                              d.get("platform", ""),
                              d.get("sample", "")),
                             chip))

        def different_chips(li):
            # Returns a map, where keys are different elements in li and
            # values their chip ids
            dc = defaultdict(list)
            for a, chip in li:
                dc[a].append(chip)
            return dc

        typeswchips = different_chips(elements)  # types with counts

        for (treatment, cond, platform, strain), cchips in typeswchips.items():
            self.chipsl.append(cchips)
            num = len(cchips)
            experiment = [strain, treatment, cond, platform,
                          str(num), ','.join(cchips)]
            self.items.append(QTreeWidgetItem(
                self.experimentsWidget, experiment))

        for i in range(5):
            self.experimentsWidget.resizeColumnToContents(i)

        self.progressBarFinished()

    def SearchUpdate(self):

        def matches(item, text):
            return any(text.upper() in str(item.text(i)).upper()
                       for i in range(item.columnCount()))

        for item in self.items:
            item.setHidden(any(not matches(item, s)
                               for s in self.searchString.split()))

    def Commit(self):
        if not self.dbc:
            self.Connect()

        start = time.time()

        pb = gui.ProgressBar(self, iterations=1000)

        ids = []
        for item in self.experimentsWidget.selectedItems():
            ids += str(item.text(5)).split(",")

        table = self.dbc.get_single_data(
            ids=ids, exclude_constant_labels=self.excludeconstant,
            callback=pb.advance)

        end = int(time.time() - start)

        pb.finish()

        data_hints.set_hint(table, "taxid", "352472", 10.0)
        data_hints.set_hint(table, "genesinrows", False, 10.0)

        self.send("Data", table)


def test_main(argv=sys.argv):
    app = QtGui.QApplication(argv)
    w = OWDicty()
    w.show()
    r = app.exec_()
    w.saveSettings()
    return r

if __name__ == "__main__":
    sys.exit(test_main())
