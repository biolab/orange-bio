"""
<name>Array Express</name>
<description>Array Express datasets<description>

"""

import sys
import os
from datetime import date

import Orange

from OWWidget import *
import OWGUI

from Orange.bio import obiArrayExpress


class OWArrayExpress(OWWidget):
    settingsList = ["current_experiement", "search_string"]

    HEADER_LABELS = ["ID", "Title", "Species", "Assays", "Date"]

    def __init__(self, parent=None, signalManager=None, title="Array Express"):
        OWWidget.__init__(self, parent, signalManager, title)

        self.outputs = [("Data Table", Orange.data.Table)]
        self.current_experiment = None
        self.search_string = ""

        self.loadSettings()

        #####
        # GUI
        #####

        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info = OWGUI.widgetLabel(box, "\n")

        OWGUI.rubber(self.controlArea)
        OWGUI.button(self.controlArea, self, "Commit", callback=self.commit)

        self.experiments_view = QTreeView(self)
        self.experiments_view.setSortingEnabled(True)
        self.experiments_view.viewport().setMouseTracking(True)
        self.experiments_view.setItemDelegateForColumn(
            0, OWGUI.LinkStyledItemDelegate(self.experiments_view)
        )

        model = QStandardItemModel()
        model.setHorizontalHeaderLabels(self.HEADER_LABELS)

        self.experiments_view.setModel(model)

        self.mainArea.layout().addWidget(self.experiments_view)

        self.setEnabled(False)
        QTimer.singleShot(5, self.fill_experiments)

    def fill_experiments(self):
        self.connection = obiArrayExpress.ArrayExpressConnection()
        self.all_experiments = []
        res = obiArrayExpress.query_experiments(gxa=True)  # only gxa for now
        res = res["experiments"]
        experiments = res["experiment"]
        model = QStandardItemModel(self)
        model.setHorizontalHeaderLabels(self.HEADER_LABELS)

        for exp in experiments:
            accession = str(exp["accession"])
            species = str(exp.get("species", ""))
            title = exp.get("name", "")
            assays = str(exp["assays"])
            date = exp.get("releasedate", "") or ""

            row = map(PyStandardItem,
                      [accession, title, species, assays, date])
            url = "http://www.ebi.ac.uk/arrayexpress/experiments/" + accession
            row[0].setData(QVariant(url), OWGUI.LinkRole)
            if not exp.get("processeddatafiles", {}).get("available", False):
                continue

            model.appendRow(row)

        self.experiments_view.setModel(model)

        self.info.setText("%i experiments" % model.rowCount())
        self.setEnabled(True)

    def commit(self):
        selected = self.experiments_view.selectionModel().selectedRows()

        if selected:
            i = selected[0].row()
            model = self.experiments_view.model()
            item = model.item(i)
            accession = str(item.data(Qt.DisplayRole).toPyObject())
            experiment = obiArrayExpress.ArrayExpressExperiment(accession)
            table = experiment.fgem_to_table()

            self.send("Data Table", table)


class PyStandardItem(QStandardItem):
    def __lt__(self, other):
        return self.data(Qt.DisplayRole).toPyObject() < \
                other.data(Qt.DisplayRole).toPyObject()


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    w = OWArrayExpress()
    w.show()
    app.exec_()
