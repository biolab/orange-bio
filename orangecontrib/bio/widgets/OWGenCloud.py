"""
<name>GenCloud</name>
<description>Access expression data from the GenCloud platform.</description>
<icon>icons/GenCloud.svg</icon>
<priority>36</priority>
"""

from __future__ import absolute_import

import Orange
import genapi

import sys, os
from collections import defaultdict
import math
from datetime import date

from Orange.orng import orngEnviron
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

import orangecontrib.bio.widgets.OWPIPAx as OWPIPAx

import requests
from orangecontrib.bio import obiDicty

NAME = "GenCloud"
DESCRIPTION = "Access expression data from the GenCloud platform."
ICON = "icons/GenCloud.svg"
PRIORITY = 36

INPUTS = []
OUTPUTS = [("Example table", Orange.data.Table)]

median = obiDicty.median
transformValues = obiDicty.transformValues
averageAttributes = obiDicty.averageAttributes
example_tables = obiDicty.example_tables
CallBack = obiDicty.CallBack

def to_text(x):
    if isinstance(x, dict):
        return str(x["value"])
    else:
        return x

class GenCountConnectionException(Exception):
    pass

class Genesis(object):

    #ignore when joining replicates
    IGNORE_REPLICATE = ["Replicate", "id", "ID", "Name" ]

    #FIXME get this from the server
    LABELS = [("static.name", "Name"),
          ("input.alignment.input.genome.static.name", "Species"),
          ("input.alignment.input.reads.var.sample.strain", "Strain"),
          ("Experiment", "Experiment"),
          ('input.alignment.input.reads.var.sample.genotype', "Genotype"),
          ("input.alignment.input.reads.var.sample.treatment", "Treatment"),
          ("input.alignment.input.reads.var.sample.growth", "Growth"),
          ("input.alignment.input.reads.var.sample.timepoint", "Timepoint"),
          ("input.alignment.input.reads.var.replicates.replicate", "Replicate"),
          ("unique_id", "ID"),
          ("adapter_type", "Adapter"),
          ("experimenter", "Experimenter"),
          ("band", "Band"),
          ("polya", "Polya"),
          ("primer", "Primer"),
          ("shearing", "Shearing")]

    def __init__(self, address, cache=None,
                 username="", password="", connect=True):

        """
        :param str address: The address of the API. 
        :param str username:
        :param str password: Login info; None for public access.
        :param CacheSQLite cache: A cache that stores results locally (an
            :obj:`CacheSQLite`).
        """
        self.address = address
        self.buffer = cache
        self.username = username
        self.password = password
        self._gen = None
        if connect:
            self._gen = genapi.GenCloud(username, password, address)
        self._project = None
        self.projectid = None

    @property
    def gen(self):
        if not self._gen:
            try:
                self._gen = genapi.GenCloud(self.username, self.password, self.address)
            except:
                raise GenCountConnectionException("Connection needed")
        return self._gen

    @property
    def project(self):
        if not self._project or self._project.id != self.projectid:
            self._project = self.gen.projects()[self.projectid]
        return self._project

    def projects(self, reload=False, bufver="0"):
        def a(self):
            return { k:str(p) for k,p in self.gen.projects().items() }
        return self._buffer_fn("projects" + "|||" + self.username, bufver, reload, a, self)

    def result_types(self, reload=False, bufver="0"):
        """Return a list of available result types.
        """
        def a(self):
            objects = self.project.objects(type__startswith='data:expression').values()
            types = set()
            for o in objects:
                an = o.annotation
                for path, a in an.iteritems():
                    if path.startswith('output') and a['type'] == 'basic:file:' \
                        and not path.startswith('output.proc.'):
                            types.add(a["name"])
            return sorted(types)
        return self._buffer_fn(self.projectid + "|||" + self.username + "|||" + "result_types", bufver, reload, a, self)

    def results_list(self, rtype, reload=False, bufver="0"):
        """Return a list of available gene expressions for a specific
        result type. Returns a dictionary, where the keys are ID
        and values are dictionaries of sample annotations.

        :param str rtype: Result type to use (see :obj:`result_types`).
        """
        def a(self):
            objects = self.project.objects(type__startswith='data:expression').values()
            rdict = {}
            for o in objects:
                an = o.annotation
                ok = False
                for path, a in an.iteritems():
                    if path.startswith('output') and a['type'] == 'basic:file:' \
                        and a["name"] == rtype:
                            ok = True  
                if ok:
                    rdict[o.id] = an
                    rdict[o.id]["date_modified"] = o.date_modified
            return rdict
        return self._buffer_fn(self.projectid + "|||" + self.username + "|||" + "results_list"  + "|||" + str(rtype), bufver, reload, a, self)

    def _from_buffer(self, addr):
        return self.buffer.get(self.address + "|v1||" + addr)

    def _to_buffer(self, addr, cont, version="0", autocommit=True):
        if self.buffer:
            return self.buffer.add(self.address + "|v1||" + addr, cont, version=version, autocommit=autocommit)

    def _buffer_commit(self):
        if self.buffer:
            self.buffer.commit()

    def _buffer_fn(self, bufkey, bufver, reload, fn, *args, **kwargs):
        """
        If bufkey is already present in buffer, return its contents.
        If not, run function with arguments and save its result
        into the buffer.
        """
        if self._in_buffer(bufkey) == bufver and reload == False:
            #print "IN", bufkey
            res = self._from_buffer(bufkey)
        else:
            #print "NOT IN", bufkey
            res = fn(*args, **kwargs)
            self._to_buffer(bufkey, res, bufver)
        return res

    def _in_buffer(self, addr):
        if self.buffer:
            return self.buffer.contains(self.address + "|v1||" + addr)
        else:
            return False

    def download(self, ids, rtype, reload=False, bufver="0"):
        objdic = self.results_list(rtype)

        downloads = [] #what to download
        for id in ids:
            o = objdic[id]
            field = None
            for path, a in o.iteritems():
                if path.startswith('output') and a['type'] == 'basic:file:' \
                    and not path.startswith('output.proc.'):
                        if a["name"] == rtype:
                            field = a["value"]["file"]
            downloads.append((id, field))

        bufverfn = (lambda x: bufver) if isinstance(bufver, basestring) else bufver

        unbuffered = [] #what is missing
        for id,field in downloads:
            if not self._in_buffer(id + "|||" + rtype) == bufverfn(id) or reload:
                unbuffered.append((id, field))
        unbufferedset = set(unbuffered)

        newgen = [].__iter__()
        if unbuffered:
            newgen = self.gen.download(unbuffered)
        for id,field in downloads:
            if (id, field) in unbufferedset:
                response = newgen.next()
                out = []
                for l in response.text.split('\n')[1:]:
                    if l:
                        gene, val = l.split('\t')
                        out.append((str(gene), str(val)))
                self._to_buffer(id + "|||" + rtype, out, version=bufverfn(id), autocommit=True)
                yield out
            else:
                yield self._from_buffer(id + "|||" + rtype)

    def get_data(self, ids=None, result_type=None,
                 exclude_constant_labels=False, average=median,
                 callback=None, bufver="0", transform=None,
                 allowed_labels=None, reload=False):
        """
        Return data in a :obj:`Orange.data.Table`. Each feature represents 
        a sample and each row is a gene. The feature's ``.attributes`` 
        contain annotations.

        :param list ids: List of ids as returned by :obj:`results_list`
            if `result_type` is None; list of ids as returned by :obj:`mappings` 
            if `result_type` is set.

        :param str result_type: Result template type id as returned by
             :obj:`result_types`.

        :param bool exclude_constant_labels: If a label has the same value
            in whole example table, remove it.

        :param function average: Function that combines multiple reading of
            the same gene on a chip. If None, no averaging is done.
            Function should take a list of floats and return an "averaged"
            float (the default functions returns the median).

        :param function transform: A function that transforms individual values.
            It should take and return a float. Example use: logarithmic 
            transformation. Default: None.

        """

        def optcb():
            if callback:
                callback()

        cbc = CallBack(len(ids), optcb, callbacks=10)

        res_list = self.results_list(result_type, reload=reload)

        #annotations
        read = {}
        for a in ids:
            an = res_list[a]
            read[a] = []
            for key, pretty in self.LABELS:
                if key in an:
                    read[a].append((pretty, to_text(an[key])))
        cbc.end()
    
        download_func = lambda x: self.download(x, result_type, reload=reload,
                                                bufver=bufver)

        cbc = CallBack(len(ids) + 3, optcb,
                       callbacks=99 - 20)
        et = example_tables(ids, spotmap={}, callback=cbc,
                            annots=read,
                            exclude_constant_labels=exclude_constant_labels,
                            chipfn=download_func,
                            allowed_labels=allowed_labels)
        cbc.end()

        cbc = CallBack(2, optcb, callbacks=10)

        #transformation is performed prior to averaging
        if transform != None:
            transformValues(et, fn=transform)  # in place transform
            cbc()

        #if average function is given, use it to join same spotids
        if average != None:
            et = averageAttributes(et, fn=average)
            cbc()

        cbc.end()

        return et


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
bufferpath = os.path.join(orngEnviron.directoryNames["bufferDir"], "genesis")

try:
    os.makedirs(bufferpath)
except:
    pass

bufferfile = os.path.join(bufferpath, "database.sq3")

SelectionByKey = OWPIPAx.SelectionByKey
ListItemDelegate = OWPIPAx.ListItemDelegate
SelectionSetsWidget = OWPIPAx.SelectionSetsWidget
SortedListWidget = OWPIPAx.SortedListWidget
               
# Mapping from PIPAx.results_list annotation keys to Header names.
HEADER = [("_cached", "")] + Genesis.LABELS

# Index of unique_id
ID_INDEX = 10

SORTING_MODEL_LIST = \
    ["Strain", "Experiment", "Genotype",
     "Timepoint", "Growth", "Species",
     "ID", "Name", "Replicate"]


class OWGenCloud(OWWidget):
    settingsList = ["server", "excludeconstant", "username", "password",
                    "joinreplicates", "selectionSetsWidget.selections",
                    "columnsSortingWidget.sortingOrder", "currentSelection",
                    "log2", "experimentsHeaderState", "rtypei", "projecti" ]

    def __init__(self, parent=None, signalManager=None, name="Genesis"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.outputs = [("Example table", ExampleTable)]

        self.username = ""
        self.password = ""
        self.log2 = False
        self.rtypei = 0
        self.projecti = 0

        self.selectedExperiments = []
        self.buffer = obiDicty.CacheSQLite(bufferfile)

        self.searchString = ""
        self.excludeconstant = False
        self.joinreplicates = False
        self.currentSelection = None
    
        self.items = []

        self.experimentsHeaderState = \
                dict(((name, False) for _, name in HEADER[:ID_INDEX + 1]))

        self.result_types = []

        self.controlArea.setMaximumWidth(250)
        self.controlArea.setMinimumWidth(250)
    
        """
        OWGUI.button(self.controlArea, self, "Reload",
                     callback=self.Reload)
        """
        OWGUI.button(self.controlArea, self, "Clear cache",
                     callback=self.clear_cache)

        box = OWGUI.widgetBox(self.controlArea, 'Project')
        self.projectCB = OWGUI.comboBox(box, self, "projecti",
                items=[],
                callback=self.ProjectChosen)

        self.projects = []

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

        QTimer.singleShot(100, self.ConnectAndUpdate)

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
        if self.dbc:
            self.projects = sorted(self.dbc.projects().items(), key=lambda x: x[1])
            self.UpdateProjects()
            self.ProjectChosen()
            self.UpdateExperimentTypes()

    def Connect(self):
        self.error(1)
        self.warning(1)

        username = "anonymous@genialis.com"
        password = "anonymous"
        if self.username:
            username = self.username
            password = self.password

        self.dbc = None
        self.projects = []
        self.result_types = []

        try:
            self.dbc = Genesis(address="http://cloud.genialis.com/",
                                  username=username,
                                  password=password,
                                  cache=self.buffer)
        except requests.exceptions.ConnectionError:
            self.dbc = Genesis(address="http://cloud.genialis.com/", 
                                  username=username,
                                  password=password,
                                  connect=False,
                                  cache=self.buffer)
            self.warning(1, "Could not connect to server, working from cache.")
        except Exception, ex:
            self.error(1, "Wrong username or password.")

        self.UpdateProjects()
        self.UpdateExperimentTypes() #clear lists

    def Reload(self):
        self.UpdateExperiments(reload=True)

    def clear_cache(self):
        self.buffer.clear()
        self.Reload()

    def rtype(self):
        """Return selected result template type """
        if self.result_types:
            return self.result_types[self.rtypei]
        else:
            return None

    def UpdateExperimentTypes(self):
        self.expressionTypesCB.clear()
        items = [desc for desc  in self.result_types]
        self.expressionTypesCB.addItems(items)
        #do not update anything if the list is empty
        if len(self.result_types):
            self.rtypei = max(0, min(self.rtypei, len(self.result_types) - 1))

    def UpdateProjects(self):
        self.projectCB.clear()
        items = [desc for pid,desc in self.projects]
        self.projectCB.addItems(items)
        #do not update anathing if the list if empty
        if len(self.projects) > 0:
            self.projecti = max(0, min(self.projecti, len(self.projects) - 1))

    def UpdateExperiments(self, reload=False):

        self.experimentsWidget.clear()

        if not self.dbc or not self.dbc.projectid: #the connection did not succeed
            return 

        self.items = []

        self.progressBarInit()

        result_types = []
        sucind = False  # success indicator for database index

        try:
            result_types = self.dbc.result_types(reload=reload)
            sucind = True
        except Exception, ex:
            try:
                result_types = self.dbc.result_types()
                self.warning(0, "Can not access database - using cached data.")
                sucind = True
            except Exception, ex:
                self.error(0, "Can not access database.")

        if sucind:
            self.warning(0)
            self.error(0)

        self.result_types = result_types

        self.UpdateExperimentTypes()
        self.UpdateResultsList(reload=reload)

        self.progressBarFinished()
        
        if self.currentSelection:
            self.currentSelection.select(self.experimentsWidget.selectionModel())

        self.handle_commit_button()

    def ProjectChosen(self, reload=False):
        if self.projects:
            self.dbc.projectid = self.projects[self.projecti][0]
            self.UpdateExperiments(reload=reload)

    def UpdateResultsList(self, reload=False):

        results_list = {}
        results_list = self.dbc.results_list(self.rtype(), reload=reload)
        try:
            results_list = self.dbc.results_list(self.rtype(), reload=reload)
        except Exception, ex:
            try:
                results_list = self.dbc.results_list(self.rtype())
            except Exception, ex:
                self.error(0, "Can not access database.")

        self.results_list = results_list

        #softly change the view so that the selection stays the same

        items_shown = {}
        for i,item in enumerate(self.items):
            c = str(item.text(ID_INDEX))
            items_shown[c] = i

        items_to_show = set(id_ for id_ in self.results_list)

        add_items = set(items_to_show) - set(items_shown)
        delete_items = set(items_shown) - set(items_to_show)

        i = 0
        while i < self.experimentsWidget.topLevelItemCount():
            it = self.experimentsWidget.topLevelItem(i)
            if str(it.text(ID_INDEX)) in delete_items:
                self.experimentsWidget.takeTopLevelItem(i)
            else:
                i += 1

        delete_ind = set([ items_shown[i] for i in delete_items ])
        self.items = [ it for i, it in enumerate(self.items) if i not in delete_ind ]

        for r_annot in add_items:
            d = defaultdict(lambda: "?", self.results_list[r_annot])
            row_items = [""] + [to_text(d.get(key, "?")) for key, _ in HEADER[1:]]
            row_items[ID_INDEX] = r_annot

            ci = MyTreeWidgetItem(self.experimentsWidget, row_items)
            self.items.append(ci)

        for i in range(len(self.headerLabels)):
            self.experimentsWidget.resizeColumnToContents(i)

        self.wantbufver = lambda x: self.results_list[x]["date_modified"]

        self.UpdateCached()


    def UpdateCached(self):

        if self.wantbufver and self.dbc:

            for item in self.items:
                id = str(item.text(ID_INDEX))
                version = self.dbc._in_buffer(id + "|||" + self.rtype())
                value = " " if version == self.wantbufver(id) else ""
                item.setData(0, Qt.DisplayRole, QVariant(value))

    def SearchUpdate(self, string=""):
        for item in self.items:
            item.setHidden(not all(s in item \
                                   for s in self.searchString.split())
                           )

    def Commit(self):

        pb = OWGUI.ProgressBar(self, iterations=100)

        table = None

        ids = []
        for item in self.experimentsWidget.selectedItems():
            unique_id = str(item.text(ID_INDEX))
            ids.append(unique_id)

        transfn = None
        if self.log2:
            transfn = lambda x: math.log(x + 1.0, 2)

        reverse_header_dict = dict((name, key) for key, name in HEADER)

        hview = self.experimentsWidget.header()
        shownHeaders = [label for i, label in \
                        list(enumerate(self.headerLabels))[1:] \
                        if not hview.isSectionHidden(i)
                        ]

        #allowed_labels = [reverse_header_dict.get(label, label) \
        #for label in shownHeaders]
        allowed_labels = None

        #if self.joinreplicates and "id" not in allowed_labels:
        #    # need 'id' labels in join_replicates for attribute names
        #    allowed_labels.append("id")

        if len(ids):
            table = self.dbc.get_data(ids=ids, result_type=self.rtype(),
                          callback=pb.advance,
                          exclude_constant_labels=self.excludeconstant,
                          bufver=self.wantbufver,
                          transform=transfn,
                          allowed_labels=allowed_labels)

            if self.joinreplicates:
                table = obiDicty.join_replicates(table,
                    ignorenames=self.dbc.IGNORE_REPLICATE,
                    namefn=None,
                    avg=obiDicty.median
                    )

            # Sort attributes
            sortOrder = self.columnsSortingWidget.sortingOrder

            all_values = defaultdict(set)
            for at in table.domain.attributes:
                atts = at.attributes
                for name in sortOrder:
                    all_values[name].add(atts.get(reverse_header_dict[name], ""))

            isnum = {}
            for at, vals in all_values.items():
                vals = filter(None, vals)
                try:
                    for a in vals:
                        float(a)
                    isnum[at] = True
                except:
                    isnum[at] = False

            def optfloat(x, at):
                if x == "":
                    return ""
                else:
                    return float(x) if isnum[at] else x

            def sorting_key(attr):
                atts = attr.attributes
                return tuple([optfloat(atts.get(reverse_header_dict[name], ""), name) \
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
                           key=(ID_INDEX,))
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
    w = OWGenCloud()
    w.show()
    app.exec_()
    w.saveSettings()
