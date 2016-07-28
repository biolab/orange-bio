"""dictyExpress widget"""
import sys
import os
import requests_cache


from requests.exceptions import HTTPError, ConnectionError
from PyQt4.QtGui import QTreeWidget, QTreeWidgetItem, QLineEdit, QLabel, QFont, QFrame, QApplication
from PyQt4.QtCore import Qt
from Orange.data import Table
from Orange.widgets.widget import OWWidget
from Orange.widgets import gui, settings
from orangecontrib.bio import resolwe
from .. import dicty
from ..utils import environ

CallBack = dicty.CallBack
Labels = [
    (" ", " "),
    ("var.project", "Project"),
    ("static.name", "Experiment"),
    ("static.cite", "Citation"),
    ("var.growth", "Growth"),
    ("var.treatment", "Treatment"),
    ("var.strain", "Parental strain")]


# Creates line separator
def h_line():
    line = QFrame()
    line.setFrameShape(QFrame.HLine)
    line.setFrameShadow(QFrame.Sunken)
    return line

#  Support cache with requests_cache module
cache_path = os.path.join(environ.buffer_dir, "resolwe")

try:
    os.makedirs(cache_path)
except OSError:
    pass

cache_file = os.path.join(cache_path, 'dictyExpress_cache')
requests_cache.install_cache(cache_name=cache_file, backend='sqlite')


class OWdictyExpress(OWWidget):

    name = "dictyExpress"
    description = "Time-course gene expression data"
    icon = "../widgets/icons/logo_dE.png"
    want_main_area = True
    priority = 37

    inputs = []
    outputs = [("Data", Table)]

    username = settings.Setting('anonymous@genialis.com')
    password = settings.Setting('')
    setTimeVariable = settings.Setting(False)

    def __init__(self):
        super().__init__()

        self.res = None
        self.server = 'https://dictyexpress.research.bcm.edu'
        self.headerLabels = [x[1] for x in Labels]
        self.searchString = ""
        self.items = []
        self.lastSelected = None  # store last selected customTreeItem

        # Login Section

        box = gui.widgetBox(self.controlArea, 'Login')

        self.namefield = gui.lineEdit(box, self, "username", "Username:", labelWidth=100, orientation='horizontal',
                                      callback=self.auth_changed)

        self.passfield = gui.lineEdit(box, self, "password", "Password:", labelWidth=100, orientation='horizontal',
                                      callback=self.auth_changed)

        self.passfield.setEchoMode(QLineEdit.Password)

        self.controlArea.layout().addWidget(h_line())

        self.time_var_checkBox = gui.checkBox(self.controlArea, self, "setTimeVariable",
                                              "Set Time variable")
        self.time_var_checkBox.setToolTip('Create new column where each row represents one time point')

        self.controlArea.layout().addWidget(h_line())

        self.commit_button = gui.button(self.controlArea, self, "Commit", callback=self.commit)
        self.handle_commit_button(False)

        self.refresh_button = gui.button(self.controlArea, self, "Refresh", callback=self.refresh)
        self.handle_cache_button(True)

        gui.rubber(self.controlArea)

        # Experiment Section

        label = QLabel("Available projects:")
        my_font = QFont()
        my_font.setBold(True)
        label.setFont(my_font)
        self.mainArea.layout().addWidget(label)

        self.mainArea.layout().addWidget(h_line())

        self.filter = gui.lineEdit(self.mainArea, self, "searchString", "Filter:", callbackOnType=True,
                                   callback=self.search_update)

        self.experimentsWidget = QTreeWidget(alternatingRowColors=True,
                                             rootIsDecorated=False,
                                             uniformRowHeights=True,
                                             sortingEnabled=True)

        self.experimentsWidget.setItemDelegateForColumn(
            0, gui.IndicatorItemDelegate(self, role=Qt.DisplayRole))

        self.experimentsWidget.selectionModel().selectionChanged.connect(
            self.onSelectionChanged)

        self.experimentsWidget.setHeaderLabels(self.headerLabels)
        self.mainArea.layout().addWidget(self.experimentsWidget)

        if self.username and self.password:
            self.connect()

    def auth_set(self):
        self.passfield.setDisabled(not self.username)

    def auth_changed(self):
        self.auth_set()
        self.connect()

    def refresh(self):
        self.reset()
        self.load_experiments()

    def reset(self):
        self.experimentsWidget.clear()  # clear QTreeWidget
        self.items = []
        self.lastSelected = None
        self.searchString = ""
        self.handle_commit_button(False)

    def search_update(self):
        parts = self.searchString.split()
        for item in self.items:
            item.setHidden(not all(s in item for s in parts))

    def connect(self):
        self.res = None
        self.error(1)
        self.reset()
        self.handle_commit_button(False)
        self.handle_cache_button(False)

        if self.username and self.password:
            try:
                self.res = resolwe.connect(self.username, self.password, self.server, 'genesis')
            except resolwe.ResolweAuthException as e:
                self.error(1, e.args[0])
            else:
                self.load_experiments()
                self.handle_cache_button(True)

    def load_experiments(self):
        with self.progressBar(100) as progress:
            if self.res:
                try:
                    experiments = self.res.fetch_etc_objects(progress.advance)  # returns a list of experiments
                except ConnectionError as e:
                    self.error(1, e.args[0])
                else:
                    self.load_tree_items(experiments)

    def load_tree_items(self, list_of_exp):
        self.items = [CustomTreeItem(self.experimentsWidget, item) for item in list_of_exp]
        for i in range(len(self.headerLabels)):
            self.experimentsWidget.resizeColumnToContents(i)

    def onSelectionChanged(self):
        self.handle_commit_button(True)

    def handle_commit_button(self, handle):
        self.commit_button.setEnabled(handle)

    def handle_cache_button(self, handle):
        self.refresh_button.setEnabled(handle)

    def commit(self):
        self.error(1)
        with self.progressBar(100) as progress:
            selected_item = self.experimentsWidget.currentItem()  # get selected TreeItem
            if self.lastSelected:
                self.lastSelected.setData(0, Qt.DisplayRole, "")
            self.lastSelected = selected_item
            selected_item.setData(0, Qt.DisplayRole, " ")

            try:
                etc_json = self.res.download_etc_data(selected_item.gen_data_id, progress.advance)
            except (HTTPError, ConnectionError, ValueError) as e:
                self.error(1, e.args[0])
            else:
                data = self.res.etc_to_table(etc_json, self.setTimeVariable,  progress.advance)
                self.send("Data", data)


class CustomTreeItem(QTreeWidgetItem):

    def __init__(self, parent, gen_data):
        super(CustomTreeItem, self).__init__(parent)  # Init super class (QtGui.QTreeWidgetItem )

        self._gen_data = gen_data  # GenData object
        self.set_rows(self._gen_data.annotation)  # set rows in QTreeWidget

    def __contains__(self, text):
        return any(text.upper() in str(self.text(i)).upper() for i in range(self.columnCount()))

    @property
    def gen_data_id(self):
        return self._gen_data.id

    def set_rows(self, row):
        for index, label in enumerate(Labels):
            if index > 0:
                if type(row[label[0]]["value"]) == list:
                    self.setText(index, row[label[0]]["value"][0]["name"])
                else:
                    self.setText(index, row[label[0]]["value"])


def test_main():
    app = QApplication(sys.argv)
    dicty.verbose = True
    w = OWdictyExpress()
    w.show()
    r = app.exec_()
    w.saveSettings()
    return r

if __name__ == "__main__":
    sys.exit(test_main())
