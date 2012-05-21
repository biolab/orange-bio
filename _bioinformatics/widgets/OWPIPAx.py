"""
<name>PIPAx</name>
<description>Access data from PIPA RNA-Seq database.</description>
<icon>icons/PIPA.png</icon>
<priority>35</priority>
"""

from __future__ import absolute_import

import sys, os
from collections import defaultdict
import math
from datetime import date

from Orange.orng import orngEnviron
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import obiDicty

from .OWPIPA import (MyTreeWidgetItem, ListItemDelegate,
                    SelectionSetsWidget, SortedListWidget)

try:
    from ast import literal_eval
except ImportError:
    # Compatibility with Python 2.5
    literal_eval = eval


def tfloat(s):
    try:
        return float(s)
    except:
        return None


class MyTreeWidgetItem(QTreeWidgetItem):

    def __init__(self, parent, *args):
        QTreeWidgetItem.__init__(self, parent, *args)
        self.par = parent

    def __contains__(self, text):
        return any(text.upper() in str(self.text(i)).upper() \
                   for i in range(self.columnCount()))

    def __lt__(self, o1):
        col = self.par.sortColumn()
        if col in [8, 9]:  # WARNING: hardcoded column numbers
            return tfloat(self.text(col)) < tfloat(o1.text(col))
        else:
            return QTreeWidgetItem.__lt__(self, o1)


# set buffer file
bufferpath = os.path.join(orngEnviron.directoryNames["bufferDir"], "pipax")

try:
    os.makedirs(bufferpath)
except:
    pass

bufferfile = os.path.join(bufferpath, "database.sq3")


class SelectionByKey(object):
    """An object stores item selection by unique key values
    (works only for row selections in list and table models)
    Example::

        ## Save selection by unique tuple pairs (DisplayRole of column 1 and 2)
        selection = SelectionsByKey(itemView.selectionModel().selection(),
                                    key = (1,2))
        ## restore selection (Possibly omitting rows not present in the model)
        selection.select(itemView.selectionModel())

    """

    def __init__(self, itemSelection, name="", key=(0,)):
        self._key = key
        self.name = name
        self._selected_keys = []
        if itemSelection:
            self.setSelection(itemSelection)

    def _row_key(self, model, row):
        def key(row, col):
            return str(model.data(model.index(row, col),
                                  Qt.DisplayRole).toString())

        return tuple(key(row, col) for col in self._key)

    def setSelection(self, itemSelection):
        self._selected_keys = [self._row_key(ind.model(), ind.row()) \
                               for ind in itemSelection.indexes() \
                               if ind.column() == 0]

    def select(self, selectionModel):
        model = selectionModel.model()
        selectionModel.clear()
        for i in range(model.rowCount()):
            if self._row_key(model, i) in self._selected_keys:
                selectionModel.select(model.index(i, 0),
                    QItemSelectionModel.Select | QItemSelectionModel.Rows)

    def __len__(self):
        return len(self._selected_keys)


# Mapping from PIPAx.results_list annotation keys to Header names.
HEADER = [("_cached", ""),
          ("data_name", "Name"),
          ("species_name", "Species"),
          ("strain", "Strain"),
          ("Experiment", "Experiment"),
          ("genotype", "Genotype"),
          ("treatment", "Treatment"),
          ("growth", "Growth"),
          ("tp", "Timepoint"),
          ("replicate", "Replicate"),
          ("unique_id", "ID"),
          ("date_rnaseq", "Date RNAseq"),
          ("adapter_type", "Adapter"),
          ("experimenter", "Experimenter"),
          ("band", "Band"),
          ("polya", "Polya"),
          ("primer", "Primer"),
          ("shearing", "Shearing")
          ]

# Index of unique_id
ID_INDEX = 10

# Index of 'date_rnaseq'
DATE_INDEX = 11

SORTING_MODEL_LIST = \
    ["Strain", "Experiment", "Genotype",
     "Timepoint", "Growth", "Species",
     "ID", "Name", "Replicate"]


class OWPIPAx(OWWidget):
    settingsList = ["server", "excludeconstant", "username", "password",
                    "joinreplicates", "selectionSetsWidget.selections",
                    "columnsSortingWidget.sortingOrder", "currentSelection",
                    "log2", "experimentsHeaderState", "rtypei"]

    def __init__(self, parent=None, signalManager=None, name="PIPAx"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.outputs = [("Example table", ExampleTable)]

        self.username = ""
        self.password = ""
        self.log2 = False
        self.rtypei = 0

        self.selectedExperiments = []
        self.buffer = obiDicty.BufferSQLite(bufferfile)

        self.searchString = ""
        self.excludeconstant = False
        self.joinreplicates = False
        self.currentSelection = None

        self.experimentsHeaderState = \
                dict(((name, False) for _, name in HEADER[:ID_INDEX + 1]))

        self.result_types = []
        self.mappings = {}

        self.controlArea.setMaximumWidth(250)
        self.controlArea.setMinimumWidth(250)

        OWGUI.button(self.controlArea, self, "Reload",
                     callback=self.Reload)
        OWGUI.button(self.controlArea, self, "Clear cache",
                     callback=self.clear_cache)

        b = OWGUI.widgetBox(self.controlArea, "Experiment Sets")
        self.selectionSetsWidget = SelectionSetsWidget(self)
        self.selectionSetsWidget.setSizePolicy(QSizePolicy.Preferred,
                                               QSizePolicy.Maximum)
        b.layout().addWidget(self.selectionSetsWidget)

        OWGUI.separator(self.controlArea)

        b = OWGUI.widgetBox(self.controlArea, "Sort output columns")
        self.columnsSortingWidget = SortedListWidget(self)
        self.columnsSortingWidget.setSizePolicy(QSizePolicy.Preferred,
                                                QSizePolicy.Maximum)
        b.layout().addWidget(self.columnsSortingWidget)
        sorting_model = QStringListModel(SORTING_MODEL_LIST)
        self.columnsSortingWidget.setModel(sorting_model)

        self.columnsSortingWidget.sortingOrder = \
                ["Strain", "Experiment", "Genotype", "Timepoint"]

        OWGUI.separator(self.controlArea)

        box = OWGUI.widgetBox(self.controlArea, 'Expression Type')
        self.expressionTypesCB = OWGUI.comboBox(box, self, "rtypei",
                items=[],
                callback=self.UpdateResultsList)

        OWGUI.checkBox(self.controlArea, self, "excludeconstant",
                       "Exclude labels with constant values"
                       )

        OWGUI.checkBox(self.controlArea, self, "joinreplicates",
                       "Average replicates (use median)"
                       )

        OWGUI.checkBox(self.controlArea, self, "log2",
                       "Logarithmic (base 2) transformation"
                       )

        self.commit_button = OWGUI.button(self.controlArea, self, "&Commit",
                                          callback=self.Commit)
        self.commit_button.setDisabled(True)

        OWGUI.rubber(self.controlArea)

        box = OWGUI.widgetBox(self.controlArea, "Authentication")

        OWGUI.lineEdit(box, self, "username", "Username:",
                       labelWidth=100,
                       orientation='horizontal',
                       callback=self.AuthChanged)

        self.passf = OWGUI.lineEdit(box, self, "password", "Password:",
                                    labelWidth=100,
                                    orientation='horizontal',
                                    callback=self.AuthChanged)

        self.passf.setEchoMode(QLineEdit.Password)

        OWGUI.lineEdit(self.mainArea, self, "searchString", "Search",
                       callbackOnType=True,
                       callback=self.SearchUpdate)

        self.headerLabels = [t[1] for t in HEADER]

        self.experimentsWidget = QTreeWidget()
        self.experimentsWidget.setHeaderLabels(self.headerLabels)
        self.experimentsWidget.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.experimentsWidget.setRootIsDecorated(False)
        self.experimentsWidget.setSortingEnabled(True)

        contextEventFilter = OWGUI.VisibleHeaderSectionContextEventFilter(
                            self.experimentsWidget, self.experimentsWidget
                            )

        self.experimentsWidget.header().installEventFilter(contextEventFilter)
        self.experimentsWidget.setItemDelegateForColumn(0,
                    OWGUI.IndicatorItemDelegate(self, role=Qt.DisplayRole)
                    )

        self.experimentsWidget.setAlternatingRowColors(True)

        self.connect(self.experimentsWidget.selectionModel(),
                 SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),
                 self.onSelectionChanged)

        self.selectionSetsWidget.setSelectionModel(
                            self.experimentsWidget.selectionModel()
                            )

        self.mainArea.layout().addWidget(self.experimentsWidget)

        self.loadSettings()

        self.restoreHeaderState()

        self.connect(self.experimentsWidget.header(),
                     SIGNAL("geometriesChanged()"),
                     self.saveHeaderState)

        self.dbc = None

        self.AuthSet()

        QTimer.singleShot(100, self.UpdateExperiments)

        self.resize(800, 600)

    def AuthSet(self):
        if len(self.username):
            self.passf.setDisabled(False)
        else:
            self.passf.setDisabled(True)

    def AuthChanged(self):
        self.AuthSet()
        self.ConnectAndUpdate()

    def ConnectAndUpdate(self):
        self.Connect()
        self.UpdateExperiments(reload=True)

    def Connect(self):
        self.error(1)
        self.warning(1)

        def en(x):
            return x if len(x) else None

        self.dbc = obiDicty.PIPAx(buffer=self.buffer,
                                  username=en(self.username),
                                  password=self.password)

        #check password
        if en(self.username) != None:
            try:
                self.dbc.mappings(reload=True)
            except obiDicty.AuthenticationError:
                self.error(1, "Wrong username or password")
                self.dbc = None
            except Exception, ex:
                print "Error when contacting the PIPA database", ex
                import traceback
                print traceback.format_exc()
                try:  # maybe cached?
                    self.dbc.mappings()
                    self.warning(1, "Can not access database - using cached data.")
                except Exception, ex:
                    self.dbc = None
                    self.error(1, "Can not access database.")

    def Reload(self):
        self.UpdateExperiments(reload=True)

    def clear_cache(self):
        self.buffer.clear()
        self.Reload()

    def rtype(self):
        """Return selected result template type """
        if self.result_types:
            return self.result_types[self.rtypei][0]
        else:
            return "-1"

    def UpdateExperimentTypes(self):
        self.expressionTypesCB.clear()
        items = [desc for _, desc  in self.result_types]
        self.expressionTypesCB.addItems(items)
        self.rtypei = max(0, min(self.rtypei, len(self.result_types) - 1))

    def UpdateExperiments(self, reload=False):
        self.experimentsWidget.clear()
        self.items = []

        self.progressBarInit()

        if not self.dbc:
            self.Connect()

        mappings = {}
        result_types = []
        sucind = False  # success indicator for database index

        try:
            mappings = self.dbc.mappings(reload=reload)
            result_types = self.dbc.result_types(reload=reload)
            sucind = True
        except Exception, ex:
            try:
                mappings = self.dbc.mappings()
                result_types = self.dbc.result_types()
                self.warning(0, "Can not access database - using cached data.")
                sucind = True
            except Exception, ex:
                self.error(0, "Can not access database.")

        if sucind:
            self.warning(0)
            self.error(0)

        self.mappings = mappings
        self.result_types = result_types

        self.UpdateExperimentTypes()

        results_list = {}
        try:
            results_list = self.dbc.results_list(self.rtype(), reload=reload)
        except Exception, ex:
            try:
                results_list = self.dbc.results_list(self.rtype())
            except Exception, ex:
                self.error(0, "Can not access database.")

        self.results_list = results_list
        mappings_key_dict = dict(((m["data_id"], m["id"]), key) \
                                 for key, m in mappings.items())

        def mapping_unique_id(annot):
            """Map annotations dict from results_list to unique
            `mappings` ids.
            """
            data_id, mappings_id = annot["data_id"], annot["mappings_id"]
            return mappings_key_dict[data_id, mappings_id]

        elements = []
        pos = 0

        for r_id, r_annot in self.results_list.items():
            pos += 1
            d = defaultdict(lambda: "?", r_annot)
            row_items = [""] + [d.get(key, "?") for key, _ in HEADER[1:]]
            date_string = row_items[DATE_INDEX]
            try:
                time_dict = literal_eval(date_string)
            except Exception:
                time_dict = {}

            if time_dict and "dateUTC" in time_dict and \
                    "monthUTC" in time_dict and "fullYearUTC" in time_dict:
                date_rna = date(time_dict["fullYearUTC"],
                                time_dict["monthUTC"] + 1,  # Why is month 0 based?
                                time_dict["dateUTC"])

                row_items[DATE_INDEX] = date_rna.strftime("%x")

            row_items[ID_INDEX] = mapping_unique_id(r_annot)
            elements.append(row_items)

            ci = MyTreeWidgetItem(self.experimentsWidget, row_items)

            self.items.append(ci)

        for i in range(len(self.headerLabels)):
            self.experimentsWidget.resizeColumnToContents(i)

        # which is the ok buffer version
        # FIXME: what attribute to use for version?
        self.wantbufver = \
            lambda x, ad=self.results_list: \
                defaultdict(lambda: "?", ad[x])["date"]

        self.wantbufver = lambda x: "0"

        self.UpdateCached()

        self.progressBarFinished()

        if self.currentSelection:
            self.currentSelection.select(self.experimentsWidget.selectionModel())

        self.handle_commit_button()

    def UpdateResultsList(self, reload=False):
        results_list = {}
        try:
            results_list = self.dbc.results_list(self.rtype(), reload=reload)
        except Exception:
            results_list = self.dbc.results_list(self.rtype())
        self.results_list = results_list
        self.UpdateCached()

    def UpdateCached(self):
        if self.wantbufver and self.dbc:
            fn = self.dbc.download_key_function()
            result_id_key = dict(((m["data_id"], m["mappings_id"]), key) \
                                 for key, m in self.results_list.items())

            for item in self.items:
                c = str(item.text(10))
                mapping = self.mappings[c]
                data_id, mappings_id = mapping["data_id"], mapping["id"]
                r_id = result_id_key[data_id, mappings_id]
                # Get the buffered version
                buffered = self.dbc.inBuffer(fn(r_id))
                value = " " if buffered == self.wantbufver(r_id) else ""
                item.setData(0, Qt.DisplayRole, QVariant(value))

    def SearchUpdate(self, string=""):
        for item in self.items:
            item.setHidden(not all(s in item \
                                   for s in self.searchString.split())
                           )

    def Commit(self):
        if not self.dbc:
            self.Connect()

        pb = OWGUI.ProgressBar(self, iterations=100)

        table = None

        ids = []
        for item in self.experimentsWidget.selectedItems():
            unique_id = str(item.text(10))
            annots = self.mappings[unique_id]
            ids.append((annots["data_id"], annots["id"]))

        transfn = None
        if self.log2:
            transfn = lambda x: math.log(x + 1.0, 2)

        reverse_header_dict = dict((name, key) for key, name in HEADER)

        hview = self.experimentsWidget.header()
        shownHeaders = [label for i, label in \
                        list(enumerate(self.headerLabels))[1:] \
                        if not hview.isSectionHidden(i)
                        ]

        allowed_labels = [reverse_header_dict.get(label, label) \
                          for label in shownHeaders]

        if self.joinreplicates and "id" not in allowed_labels:
            # need 'id' labels in join_replicates for attribute names
            allowed_labels.append("id")

        if len(ids):
            table = self.dbc.get_data(ids=ids, result_type=self.rtype(),
                          callback=pb.advance,
                          exclude_constant_labels=self.excludeconstant,
#                          bufver=self.wantbufver,
                          transform=transfn,
                          allowed_labels=allowed_labels)

            if self.joinreplicates:
                table = obiDicty.join_replicates(table,
                    ignorenames=["replicate", "data_id", "mappings_id",
                                 "data_name", "id", "unique_id"],
                    namefn=None,
                    avg=obiDicty.median
                    )

            # Sort attributes
            sortOrder = self.columnsSortingWidget.sortingOrder

            def sorting_key(attr):
                atts = attr.attributes
                return tuple([atts.get(reverse_header_dict[name], "") \
                              for name in sortOrder])

            attributes = sorted(table.domain.attributes,
                                key=sorting_key)

            domain = orange.Domain(attributes, table.domain.classVar)
            domain.addmetas(table.domain.getmetas())
            table = orange.ExampleTable(domain, table)

            from Orange.orng.orngDataCaching import data_hints
            data_hints.set_hint(table, "taxid", "352472")
            data_hints.set_hint(table, "genesinrows", False)

            self.send("Example table", table)

            self.UpdateCached()

        pb.finish()

    def onSelectionChanged(self, selected, deselected):
        self.handle_commit_button()

    def handle_commit_button(self):
        self.currentSelection = \
            SelectionByKey(self.experimentsWidget.selectionModel().selection(),
                           key=(1, 2, 3, 10))
        self.commit_button.setDisabled(not len(self.currentSelection))

    def saveHeaderState(self):
        hview = self.experimentsWidget.header()
        for i, label in enumerate(self.headerLabels):
            self.experimentsHeaderState[label] = hview.isSectionHidden(i)

    def restoreHeaderState(self):
        hview = self.experimentsWidget.header()
        state = self.experimentsHeaderState
        for i, label in enumerate(self.headerLabels):
            hview.setSectionHidden(i, state.get(label, True))
            self.experimentsWidget.resizeColumnToContents(i)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    obiDicty.verbose = True
    w = OWPIPAx()
    w.show()
    app.exec_()
    w.saveSettings()
