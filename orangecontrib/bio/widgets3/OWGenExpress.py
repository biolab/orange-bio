import sys
import os
import math
import io
import gzip

from collections import defaultdict

from AnyQt.QtWidgets import (
    QSizePolicy, QTreeWidget, QTreeWidgetItem, QLineEdit,
)
from AnyQt.QtCore import Qt, QSize, QTimer, QStringListModel

import requests
import genesis

import Orange.data

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.datacaching import data_hints

from .. import dicty
from ..utils import environ, compat
from .OWPIPAx import SelectionByKey, SelectionSetsWidget, SortedListWidget

median = dicty.median
transformValues = dicty.transformValues
averageAttributes = dicty.averageAttributes
example_tables = dicty.example_tables
CallBack = dicty.CallBack


def to_text(x):
    if isinstance(x, dict):
        return str(x["value"])
    else:
        return x


class GenCountConnectionException(Exception):
    pass


class Genesis(object):

    #ignore when joining replicates
    IGNORE_REPLICATE = ["Replicate", "id", "ID", "Name"]

    #FIXME get this from the server
    LABELS = [
        ("input.alignment.input.reads.var.experiment", "Experiment"),
        ("input.alignment.input.reads.var.sample.time", "Timepoint"),
        ("input.alignment.input.reads.var.replicates.replicate", "Replicate"),
        ("input.alignment.input.reads.var.sample.strain", "Strain"),
        ('input.alignment.input.reads.var.sample.genotype', "Genotype"),
        ("input.alignment.input.reads.var.sample.treatment", "Treatment"),
        ("input.alignment.input.reads.var.sample.growth", "Growth"),
        ("input.alignment.input.genome.static.name", "Genome"),
        ("static.name", "Name"),
        ("unique_id", "ID"),
        ("adapter_type", "Adapter"),
        ("experimenter", "Experimenter"),
        ("band", "Band"),
        ("polya", "Polya"),
        ("primer", "Primer"),
        ("shearing", "Shearing"),
    ]

    def __init__(self, address, cache=None,
                 username='anonymous@genialis.com', password='anonymous',
                 connect=True):

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
            self._gen = genesis.Genesis(username, password, address)
        self._project = None
        self.projectid = None

    @property
    def gen(self):
        if not self._gen:
            self._gen = genesis.Genesis(self.username, self.password,
                                        self.address)
        return self._gen

    @property
    def project(self):
        if not self._project or self._project.id != self.projectid:
            self._project = self.gen.projects()[self.projectid]
        return self._project

    def projects(self, reload=False, bufver="0"):
        def a(self):
            return {k: str(p) for k, p in self.gen.projects().items()}
        return self._buffer_fn("projects" + "|||" + self.username, bufver,
                               reload, a, self)

    def result_types(self, reload=False, bufver="0"):
        """Return a list of available result types."""
        def a(self):
            objects = self.project.data(type__startswith='data:expression:')
            types = {}
            for o in objects:
                an = o.annotation
                for path, a in an.items():
                    if path.startswith('output') and \
                            a['type'] == 'basic:file:' and \
                            not path.startswith('output.proc.'):
                        types[a['name']] = a['label']

            # Sort keys by labels
            keys = [typ for typ, _ in sorted(types.items(), key=lambda x: x[1])]
            try:
                # Push field exp to the top (exp is considered default)
                keys.remove('exp')
                keys = ['exp'] + keys
            except ValueError:
                # Ignore if no 'exp' field
                pass

            return keys, types
        return self._buffer_fn(self.projectid + "|||" + self.username +
                               "|||" + "result_types", bufver, reload,
                               a, self)

    def results_list(self, rtype, reload=False, bufver="0"):
        """Return a list of available gene expressions for a specific
        result type. Returns a dictionary, where the keys are ID
        and values are dictionaries of sample annotations.

        :param str rtype: Result type to use (see :obj:`result_types`).
        """
        def a(self):
            objects = self.project.data(type__startswith='data:expression:')
            rdict = {}
            for o in objects:
                an = o.annotation
                ok = False
                for path, a in an.items():
                    if path.startswith('output') and \
                            a['type'] == 'basic:file:' and a["name"] == rtype:
                        ok = True
                if ok:
                    rdict[o.id] = an
                    rdict[o.id]["date_modified"] = o.date_modified
            return rdict
        return self._buffer_fn(self.projectid + "|||" + self.username +
                               "|||" + "results_list" + "|||" + str(rtype),
                               bufver, reload, a, self)

    def _from_buffer(self, addr):
        return self.buffer.get(self.address + "|v1||" + addr)

    def _to_buffer(self, addr, cont, version="0", autocommit=True):
        if self.buffer:
            return self.buffer.add(self.address + "|v1||" + addr, cont,
                                   version=version, autocommit=autocommit)

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

        downloads = []  # what to download
        for id in ids:
            o = objdic[id]
            field = None
            for path, a in o.items():
                if path.startswith('output') and a['type'] == 'basic:file:' \
                    and not path.startswith('output.proc.'):
                        if a["name"] == rtype:
                            field = a["value"]["file"]
            downloads.append((id, field))

        bufverfn = (lambda x: bufver) if isinstance(bufver, str) else bufver

        unbuffered = []  # what is missing
        for id, field in downloads:
            if not self._in_buffer(id + "|||" + rtype) == bufverfn(id) or \
                    reload:
                unbuffered.append((id, field))
        unbufferedset = set(unbuffered)

        # newgen = [].__iter__()
        # if unbuffered:
        #     newgen = self.gen.download(unbuffered)
        for id, field in downloads:
            if (id, field) in unbufferedset:
                response = next(self.gen.download([id], 'output.' + rtype))
                response_gzipped = io.BytesIO(response.content)
                response_content = io.TextIOWrapper(
                    gzip.GzipFile(fileobj=response_gzipped),
                    encoding="utf-8")

                out = []
                for l in response_content.read().split('\n')[1:]:
                    if l:
                        gene, val = l.split('\t')
                        out.append((str(gene), str(val)))
                self._to_buffer(id + "|||" + rtype, out, version=bufverfn(id),
                                autocommit=True)
                yield out
            else:
                yield self._from_buffer(id + "|||" + rtype)

    def get_data(self, ids=None, result_type=None,
                 exclude_constant_labels=False, average=median,
                 callback=None, bufver="0", transform=None,
                 allowed_labels=None, reload=False, namefn=None):
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
                            allowed_labels=allowed_labels,
                            namefn=namefn)
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
    except ValueError:
        return None


class MyTreeWidgetItem(QTreeWidgetItem):

    def __init__(self, parent, *args):
        QTreeWidgetItem.__init__(self, parent, *args)
        self.par = parent

    def __contains__(self, text):
        return any(text.upper() in str(self.text(i)).upper() \
                   for i in range(self.columnCount()))

    def __lt__(self, other):
        col = self.par.sortColumn()
        if col in [TIMEPOINT_COLUMN, REPLICATE_COLUMN]:
            left = tfloat(self.text(col))
            right = tfloat(other.text(col))
            if isinstance(left, float) and isinstance(right, float):
                return left < right

        return QTreeWidgetItem.__lt__(self, other)


# set buffer file
bufferpath = os.path.join(environ.buffer_dir, "genesis")

try:
    os.makedirs(bufferpath)
except OSError:
    pass

bufferfile = os.path.join(bufferpath, "database.sq3")

# Mapping from PIPAx.results_list annotation keys to Header names.
HEADER = [("_cached", "")] + Genesis.LABELS

# Index of unique_id
ID_INDEX = 10

TIMEPOINT_COLUMN = 2
REPLICATE_COLUMN = 3

SORTING_MODEL_LIST = [
    "Strain", "Experiment", "Genotype", "Timepoint", "Growth", "Genome",
    "ID", "Name", "Replicate"]


class OWGenExpress(widget.OWWidget):
    name = "GenExpress"
    description = "Expression data from GenExpress."
    icon = "../widgets/icons/GenCloud.svg"
    priority = 36

    inputs = []
    outputs = [("Data", Orange.data.Table)]

    username = settings.Setting("anonymous")
    password = settings.Setting("")
    log2 = settings.Setting(False)
    transpose = settings.Setting(False)
    rtypei = settings.Setting(0)
    projecti = settings.Setting(0)
    serveri = settings.Setting(0)
    exnamei = settings.Setting(6)

    excludeconstant = settings.Setting(False)
    joinreplicates = settings.Setting(False)
    currentSelection = settings.Setting(None)

    experimentsHeaderState = settings.Setting({
        name: False for _, name in HEADER[:ID_INDEX + 1]}
    )

    storedSortOrder = settings.Setting([])
    storedSelections = settings.Setting([])

    def __init__(self, parent=None):
        super().__init__(parent)

        self.servers = [
            ('https://dictyexpress.research.bcm.edu/', 'dictyExpress'),
            ('https://cloud.genialis.com/', 'Genialis'),
        ]

        self.selectedExperiments = []
        self.buffer = dicty.CacheSQLite(bufferfile)

        self.searchString = ""

        self.items = []

        self.result_types = []

        self.controlArea.setMaximumWidth(250)
        self.controlArea.setMinimumWidth(250)

        box = gui.widgetBox(self.controlArea, 'Project')
        self.projectCB = gui.comboBox(
            box, self, "projecti", items=[], callback=self.ProjectChosen)

        self.projects = []

        b = gui.widgetBox(self.controlArea, "Selection bookmarks")
        self.selectionSetsWidget = SelectionSetsWidget(self)
        self.selectionSetsWidget.setSizePolicy(
            QSizePolicy.Preferred, QSizePolicy.Maximum)

        def store_selections(modified):
            if not modified:
                self.storedSelections = self.selectionSetsWidget.selections
        self.selectionSetsWidget.selectionModified.connect(store_selections)

        b.layout().addWidget(self.selectionSetsWidget)

        gui.separator(self.controlArea)

        b = gui.widgetBox(self.controlArea, "Sort output columns")
        self.columnsSortingWidget = SortedListWidget(self)
        self.columnsSortingWidget.setSizePolicy(
            QSizePolicy.Preferred, QSizePolicy.Maximum)

        box = gui.widgetBox(self.controlArea, 'Experiment name')
        self.experimentNameCB = gui.comboBox(
            box, self, "exnamei", items=SORTING_MODEL_LIST)

        b.layout().addWidget(self.columnsSortingWidget)
        sorting_model = QStringListModel(SORTING_MODEL_LIST)
        self.columnsSortingWidget.setModel(sorting_model)
        self.columnsSortingWidget.sortingOrder = self.storedSortOrder

        def store_sort_order():
            self.storedSortOrder = self.columnsSortingWidget.sortingOrder
        self.columnsSortingWidget.sortingOrderChanged.connect(store_sort_order)

        gui.separator(self.controlArea)


        box = gui.widgetBox(self.controlArea, 'Expression Type')
        self.expressionTypesCB = gui.comboBox(
            box, self, "rtypei", items=[], callback=self.UpdateResultsList)

        gui.checkBox(self.controlArea, self, "excludeconstant",
                     "Exclude labels with constant values")

        gui.checkBox(self.controlArea, self, "joinreplicates",
                     "Average replicates (use median)")

        gui.checkBox(self.controlArea, self, "log2",
                     "Logarithmic (base 2) transformation")

        gui.checkBox(self.controlArea, self, "transpose",
                     "Genes as columns")

        self.commit_button = gui.button(self.controlArea, self, "&Commit",
                                        callback=self.Commit)
        self.commit_button.setDisabled(True)

        gui.rubber(self.controlArea)

        box = gui.widgetBox(self.controlArea, 'Server')
        gui.comboBox(box, self, "serveri",
                     items=[title for url, title in self.servers],
                     callback=self.ServerChosen)

        gui.lineEdit(box, self, "username", "Username:",
                     labelWidth=100,
                     orientation='horizontal',
                     callback=self.AuthChanged)

        self.passf = gui.lineEdit(box, self, "password", "Password:",
                                  labelWidth=100,
                                  orientation='horizontal',
                                  callback=self.AuthChanged)

        self.passf.setEchoMode(QLineEdit.Password)

        gui.button(self.controlArea, self, "Clear cache",
                   callback=self.clear_cache)

        gui.lineEdit(self.mainArea, self, "searchString", "Search",
                     callbackOnType=True,
                     callback=self.SearchUpdate)

        self.headerLabels = [t[1] for t in HEADER]

        self.experimentsWidget = QTreeWidget()
        self.experimentsWidget.setHeaderLabels(self.headerLabels)
        self.experimentsWidget.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.experimentsWidget.setRootIsDecorated(False)
        self.experimentsWidget.setSortingEnabled(True)

        contextEventFilter = gui.VisibleHeaderSectionContextEventFilter(
            self.experimentsWidget, self.experimentsWidget)

        self.experimentsWidget.header().installEventFilter(contextEventFilter)
        self.experimentsWidget.setItemDelegateForColumn(
            0, gui.IndicatorItemDelegate(self, role=Qt.DisplayRole))

        self.experimentsWidget.setAlternatingRowColors(True)

        self.experimentsWidget.selectionModel().selectionChanged.connect(
            self.onSelectionChanged)

        self.selectionSetsWidget.setSelectionModel(
            self.experimentsWidget.selectionModel())
        self.selectionSetsWidget.setSelections(self.storedSelections)

        self.mainArea.layout().addWidget(self.experimentsWidget)

        self.restoreHeaderState()

        self.experimentsWidget.header().geometriesChanged.connect(
            self.saveHeaderState)

        self.dbc = None

        self.AuthSet()

        QTimer.singleShot(100, self.ConnectAndUpdate)

    def sizeHint(self):
        return QSize(800, 600)

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
            def get_data_count(project_id):
                # XXX: is there a better way?
                # Note: limit 0 would return all objects
                return self.dbc.gen.api.data.get(case_ids__contains=project_id,
                                                 type__startswith='data:expression:',
                                                 limit=1)['meta']['total_count']

            self.projects = sorted([p for p in self.dbc.projects().items() if
                                    get_data_count(p[0]) > 0], key=lambda x: x[1])
            self.UpdateProjects()
            self.ProjectChosen()
            self.UpdateExperimentTypes()

    def Connect(self):
        self.error(1)
        self.warning(1)

        username = 'anonymous@genialis.com'
        password = 'anonymous'
        url = self.servers[self.serveri][0]

        if self.username:
            username = self.username
            password = self.password

        if username.lower() in ['anonymous@genialis.com', 'anonymous']:
            username = 'anonymous@genialis.com'
            password = 'anonymous'

        self.dbc = None
        self.projects = []
        self.result_types = []

        try:
            self.dbc = Genesis(
                address=url, username=username, password=password,
                cache=self.buffer)
        except requests.exceptions.ConnectionError:
            self.dbc = Genesis(
                address=url, username=username, password=password,
                connect=False, cache=self.buffer)
            self.warning(1, "Could not connect to server, working from cache.")
        except Exception:
            self.error(1, "Wrong username or password.")

        self.UpdateProjects()
        self.UpdateExperimentTypes()  # clear lists

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
        items = [self.result_types_labels[desc] for desc in self.result_types]
        self.expressionTypesCB.addItems(items)
        #do not update anything if the list is empty
        if len(self.result_types):
            self.rtypei = max(0, min(self.rtypei, len(self.result_types) - 1))

    def UpdateProjects(self):
        self.projectCB.clear()
        items = [desc for pid, desc in self.projects]
        self.projectCB.addItems(items)
        #do not update anything if the list if empty
        if len(self.projects) > 0:
            self.projecti = max(0, min(self.projecti, len(self.projects) - 1))

    def UpdateExperiments(self, reload=False):

        self.experimentsWidget.clear()

        if not self.dbc or not self.dbc.projectid:  # the connection did not succeed
            return

        self.items = []

        self.progressBarInit()

        result_types = []
        result_types_labels = []

        sucind = False  # success indicator for database index

        try:
            result_types, result_types_labels = self.dbc.result_types(reload=reload)
            sucind = True
        except Exception:
            try:
                result_types, result_types_labels = self.dbc.result_types()
                self.warning(0, "Can not access database - using cached data.")
                sucind = True
            except Exception:
                self.error(0, "Can not access database.")

        if sucind:
            self.warning(0)
            self.error(0)

        self.result_types = result_types
        self.result_types_labels = result_types_labels

        self.UpdateExperimentTypes()
        self.UpdateResultsList(reload=reload)

        self.progressBarFinished()

        if self.currentSelection:
            self.currentSelection.select(self.experimentsWidget.selectionModel())

        self.handle_commit_button()

    def ProjectChosen(self, reload=False):
        if self.projects:
            self.dbc.projectid = self.projects[self.projecti][0]
        else:
            self.dbc.projectid = None

        self.UpdateExperiments(reload=reload)

    def ServerChosen(self):
        self.ConnectAndUpdate()

    def UpdateResultsList(self, reload=False):
        results_list = self.dbc.results_list(self.rtype(), reload=reload)
        try:
            results_list = self.dbc.results_list(self.rtype(), reload=reload)
        except Exception:
            try:
                results_list = self.dbc.results_list(self.rtype())
            except Exception:
                self.error(0, "Can not access database.")

        self.results_list = results_list

        #softly change the view so that the selection stays the same

        items_shown = {}
        for i, item in enumerate(self.items):
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

        delete_ind = set([items_shown[i] for i in delete_items])
        self.items = [it for i, it in enumerate(self.items)
                      if i not in delete_ind]

        for r_annot in add_items:
            d = defaultdict(lambda: "?", self.results_list[r_annot])
            row_items = [""] + [to_text(d.get(key, "?"))
                                for key, _ in HEADER[1:]]
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
                item.setData(0, Qt.DisplayRole, value)

    def SearchUpdate(self, string=""):
        for item in self.items:
            item.setHidden(
                not all(s in item for s in self.searchString.split()))

    def Commit(self):

        pb = gui.ProgressBar(self, iterations=100)

        table = None

        ids = []
        for item in self.experimentsWidget.selectedItems():
            unique_id = str(item.text(ID_INDEX))
            ids.append(unique_id)

        transfn = None
        if self.log2:
            transfn = lambda x: math.log(x + 1.0, 2)

        reverse_header_dict = {name: name for key, name in HEADER}
        reverse_header_dict["ID"] = "id"

        allowed_labels = None

        def namefn(a):
            name = SORTING_MODEL_LIST[self.exnamei]
            name = reverse_header_dict.get(name, "id")
            return dict(a)[name]

        if len(ids):
            table = self.dbc.get_data(
                ids=ids, result_type=self.rtype(),
                callback=pb.advance,
                exclude_constant_labels=self.excludeconstant,
                bufver=self.wantbufver,
                transform=transfn,
                allowed_labels=allowed_labels,
                namefn=namefn)

            if self.joinreplicates:
                table = dicty.join_replicates(table,
                    ignorenames=self.dbc.IGNORE_REPLICATE,
                    namefn="name",
                    avg=dicty.median,
                    fnshow=lambda x: " | ".join(map(str, x)))

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
                except ValueError:
                    isnum[at] = False

            def optfloat(x, at):
                if x == "":
                    return ""
                else:
                    return float(x) if isnum[at] else x

            def sorting_key(attr):
                atts = attr.attributes
                return tuple([optfloat(atts.get(reverse_header_dict[name], ""), name)
                              for name in sortOrder])

            attributes = sorted(table.domain.attributes, key=sorting_key)

            domain = Orange.data.Domain(
                attributes, table.domain.class_vars, table.domain.metas)

            table = Orange.data.Table.from_table(domain, table)
            table = Orange.data.Table(domain, table)

            if self.transpose:
                experiments = [at for at in table.domain.variables]
                attr = [compat.ContinuousVariable.make(ex['DDB'].value) for ex in table]
                metavars = sorted(table.domain.variables[0].attributes.keys())
                metavars = [compat.StringVariable.make(name) for name in metavars]
                domain = compat.create_domain(attr, None, metavars)
                metas = [[exp.attributes[var.name] for var in metavars] for exp in experiments]
                table = compat.create_table(domain, table.X.transpose(), None, metas)

            data_hints.set_hint(table, "taxid", "352472")
            data_hints.set_hint(table, "genesinrows", False)

            self.send("Data", table)

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


def test_main():
    from AnyQt.QtWidgets import QApplication
    app = QApplication(sys.argv)
    dicty.verbose = True
    w = OWGenExpress()
    w.show()
    r = app.exec_()
    w.saveSettings()
    return r

if __name__ == "__main__":
    sys.exit(test_main())
