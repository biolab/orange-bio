import sys
import warnings
import io
from collections import defaultdict
from functools import partial
import itertools
import xml.sax

import six

from AnyQt.QtWidgets import (
    QWidget, QComboBox, QListView, QLineEdit, QPlainTextEdit, QRadioButton,
    QButtonGroup, QGroupBox, QCheckBox, QLabel, QFrame, QTabWidget,
    QSizePolicy, QGridLayout, QVBoxLayout, QHBoxLayout, QStackedLayout,
)
from AnyQt.QtGui import QRegExpValidator

from AnyQt.QtCore import (
    Qt, QObject, QRegExp, QModelIndex, QThread, QThreadPool, QSize,
    QStringListModel, QItemSelectionModel
)
from AnyQt.QtCore import Slot

import Orange

from Orange.widgets.utils import concurrent
from Orange.widgets import widget, gui, settings

from .. import biomart


def is_hidden(tree):
    return (getattr(tree, "hidden", "false") != "false" or
            getattr(tree, "hideDisplay", "false") != "false")


class Control(object):

    """ Base mixin class for query GUI widgets
    """

    def __init__(self, tree=None, dataset=None, master=None, **kwargs):
        """
            :param tree: configuration element
            :type tree: biomart.ConfigurationNode
            :param dataset: dataset
            :type dataset: biomart.BioMartDataset
            :param master: main widget
            :type master: OWBioMart

        """
        super(Control, self).__init__(**kwargs)
        self.tree = tree
        self.dataset = dataset
        self.master = master
        self.subControls = []

        if tree is not None and isinstance(self, QObject):
            self.setObjectName(tree.internalName)
            if hasattr(tree, "description"):
                self.setToolTip(tree.description)

    def addSubControl(self, tree, control):
        self.subControls.append((tree, control))

    def registerDelayedCall(self, call):
        self.master.registerDelayedCall(call)

    def pushAction(self, action):
        self.master.pushAction(action)

    def query(self):
        return itertools.chain(*[ctrl.query() for _, ctrl in self.subControls])

    def setControlValue(self, name, value):
        if "." in name:
            name, rest = name.split(".", 1)
            controls = {tree.internalName: control
                        for tree, control in self.subControls}
            ctrl = controls.get("name", None)
            if ctrl:
                ctrl.setControlValue(rest, value)


class UnknownFilter(QWidget, Control):

    def __init__(self, tree, dataset, master, parent=None):
        QWidget.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)


class TextFieldFilter(QLineEdit, Control):

    """ A single edit line filter
    """

    def __init__(self, tree, dataset, master, parent=None):
        QLineEdit.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)

        if hasattr(tree, "regexp"):
            self.setValidator(QRegExpValidator(QRegExp(tree.regexp), self))

        if hasattr(tree, "defaultValue"):
            self.setText(tree.defaultValue)

    def get_filter(self):
        return self.tree.internalName, str(self.text())

    def query(self):
        return [("Filter", self.tree, str(self.text()))]

    def setControlValue(self, name, value):
        self.setText(value)


class IdListFilter(QWidget, Control):

    """ Multiple ids filter
    """

    def __init__(self, tree, dataset, master, parent=None):
        QWidget.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)
        self.tree = tree
        self.dataset = dataset
        self.master = master
        self.setObjectName(tree.internalName)

        self.setLayout(QGridLayout())
        self.textWidget = QPlainTextEdit()  # TODO: Subclass to receive drop events from item model views
        self.layout().addWidget(self.textWidget, 0, 0, 1, 1)

    def value(self):
        """ Return filter value for use in a query
        """
        return str(self.textWidget.toPlainText()).split()

    def get_filter(self):
        return self.tree.internalName, self.value()

    def query(self):
        return [("Filter", self.tree, self.value())]

    def setControlValue(self, name, value):
        if type(value) == list:
            value = "\n".join(value)
        self.textWidget.setPlainText(value)


class RadioBooleanFilter(QWidget, Control):

    """ Boolean filter (Only/Exclude)
    """

    def __init__(self, tree, dataset, master, parent=None):
        QWidget.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)

        self.setLayout(QVBoxLayout())
        self.buttonGroup = QButtonGroup(self)
        self.values = []
        for i, option in enumerate(tree.subelements_top("Option")):
            rb = QRadioButton(option.displayName, self)
            self.buttonGroup.addButton(rb)
            self.buttonGroup.setId(rb, i)
            self.layout().addWidget(rb)
            self.values.append(option.value)
        self.buttonGroup.button(0).setChecked(True)

    def value(self):
        return {"excluded": "%i" % self.buttonGroup.checkedId()}

    def get_filter(self):
        return self.tree.internalName, self.value()

    def query(self):
        return [("Filter", self.tree, self.value())]

    def setControlValue(self, name, value):
        for i, v in enumerate(self.values):
            if v == value:
                button = self.buttonGroup.button(i)
                button.setChecked(True)
                break


class DropDownFilter(QComboBox, Control):

    """ List menu filter
    """

    def __init__(self, tree, dataset, master, parent=None):
        QComboBox.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)

        self.options = []
        self.selectedIndex = 0

        if getattr(tree, "graph", "0") == "1":
            self.setOptions(tree.subelements("Option"))
        else:
            self.setOptions(tree.subelements_top("Option"))

        self.currentIndexChanged[int].connect(self.onIndexChange)

    def setOptions(self, options):
        self.options = []
        self.blockSignals(True)
        self.clear()
        for option in options:
            self.addItem(option.displayName)
            self.options.append(option)
        self.selectedIndex = 0
        self.blockSignals(False)
        self.registerDelayedCall(lambda: self.onIndexChange(0))

    def onIndexChange(self, index):
        if self.options:
            option = self.options[index]
            self.selectedIndex = index
            pushActions = option.subelements_top("PushAction")
            for action in pushActions:
                self.master.pushAction(action)

    def value(self):
        option = self.options[self.selectedIndex]
        return option.value

    def query(self):
        return [("Filter", self.tree, self.value())]

    def setControlValue(self, name, value):
        for i, option in enumerate(self.options):
            if option.value == value:
                self.setCurrentIndex(i)


def rows(index_list):
    return map(QModelIndex.row, index_list)


class MultiSelectListFilter(QListView, Control):

    """ Menu list filter with multiple selection
    """

    def __init__(self, tree, dataset, master, parent=None):
        QComboBox.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)

        self.setSelectionMode(QListView.ExtendedSelection)
        model = QStringListModel(self)
        self.setModel(model)
        self.setOptions(tree.subelements_top("Option"))

    def setOptions(self, options):
        self.options = []
        for option in options:
            self.options.append(option)
        self.model().setStringList(
            [option.displayName for option in self.options])

    def value(self):
        value = []
        for index in rows(self.selectedIndexes()):
            value.append(self.options[index].value)
        return ",".join(value)

    def query(self):
        return [("Filter", self.tree, self.value())]

    def setControlValue(self, name, value):
        if isinstance(value, str):
            values = value.split(",")
        else:
            values = value
        selection = self.selectionModel()
        for i, option in enumerate(self.options):
            if option.value in values:
                selection.select(self.model().index(i),
                                 QItemSelectionModel.Select)


class DropDownRadioBooleanFilter(QWidget, Control):
    """Container for multiple boolean filters
    """

    def __init__(self, tree, dataset, master, parent=None):
        QWidget.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)

        self.setLayout(QHBoxLayout())
        self.cb = QComboBox(self)

        self.layout().addWidget(self.cb)

        rblayout = QVBoxLayout()
        self.radioButtons = [QRadioButton("Only", self),
                             QRadioButton("Excluded", self)
                             ]

        for b in self.radioButtons:
            rblayout.addWidget(b)

        self.radioButtons[0].setChecked(True)

        self.layout().addLayout(rblayout)

        self.options = []

        self.setOptions(tree.subelements_top("Option"))

    def setOptions(self, options):
        self.cb.clear()
        self.options = []
        for option in options:
            self.cb.addItem(option.displayName)
            self.options.append(option)

        for op, rb in zip(self.options[0].subelements_top("Option"),
                          self.radioButtons):
            rb.setText(op.displayName)
            rb.setChecked(getattr(op, "default", "false") == "true")

    def value(self):
        return {"excluded": "0" if self.radioButtons[0].isChecked() else "1"}

    def query(self):
        filter = self.options[self.cb.currentIndex()]
        filter = biomart.FilterDescription(
            self.tree.registry, "FilterDescription",
            filter.attributes, filter.children)
        return [("Filter", filter, self.value())]

    def setControlValue(self, name, value):
        for i, option in enumerate(self.options):
            if option.internalName == name:
                self.cb.setCurrentIndex(i)
                if value == "Only":
                    self.radioButtons[0].setChecked(True)


class DropDownIdListFilter(QWidget, Control):

    """Container for multiple id list filters
    """

    def __init__(self, tree, dataset, master, parent=None):
        QWidget.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)

        self.setLayout(QVBoxLayout())
        self.setContentsMargins(0, 0, 0, 0)
        self.cb = QComboBox()
        self.idsEdit = QPlainTextEdit()

        self.layout().addWidget(self.cb)
        self.layout().addWidget(self.idsEdit)

        self.options = []
        self.setOptions(tree.subelements_top("Option"))

    def setOptions(self, options):
        self.cb.clear()
        self.options = []
        for option in options:
            self.cb.addItem(option.displayName)
            self.options.append(option)

    def value(self):
        return str(self.idsEdit.toPlainText()).split()

    def query(self):
        filter = self.options[self.cb.currentIndex()]
        filter = biomart.FilterDescription(
            self.tree.registry, "FilterDescription",
            filter.attributes, filter.children)
        return [("Filter", filter, self.value())]

    def setControlValue(self, name, value):
        if isinstance(value, list):
            value = "\n".join(value)

        for i, op in enumerate(self.options):
            if name == op.internalName:
                self.cb.setCurrentIndex(i)
                self.idsEdit.setPlainText(value)


class CollectionWidget(QGroupBox, Control):
    NO_FLAGS = 0
    FILTER_FLAG = 1
    ATTRIBUTE_FLAG = 2
    SINGLE_FILTER_FLAG = 4

    def __init__(self, tree, dataset, master, parent=None):
        QGroupBox.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)
        self.setFlat(True)
        self.attributes = []
        self.filters = []

    def setCollection(self, tree, dataset, flags=0):
        self.tree = tree
        self.dataset = dataset
        self.collectionFlags = flags


def isFilterCollection(tree):
    """ In case all AttributeDescription nodes contain pointers to filters
    (case of downstream/upstream_flank)
    """
    if tree.tag == "FilterCollection":
        return True
    return all(["pointerFilter" in desc.attributes
                for desc in tree.elements("AttributeDescription")])


class AttributeWidget(QCheckBox, Control):

    """ Single attribute selection widget
    """

    def __init__(self, tree, dataset, master, parent=None):
        QCheckBox.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)

        self.setText(getattr(tree, "displayName", ""))
        self.setChecked(getattr(tree, "default", "false") == "true")

        if hasattr(tree, "description"):
            self.setToolTip(tree.description)

    def query(self):
        if self.isChecked():
            return [("Attribute", self.tree, self.isChecked())]
        else:
            return []

    def setControlValue(self, name, value):
        self.setChecked(value == "true")


class AttributeCollectionWidget(CollectionWidget):

    def __init__(self, tree, dataset, master, parent=None):
        CollectionWidget.__init__(self, tree, dataset, master, parent)
        self.setLayout(QGridLayout(self))
        self.setCollection(tree, dataset)

    def setCollection(self, collection, dataset, flags=0):
        CollectionWidget.setCollection(self, collection, dataset, flags)
        self.setTitle(getattr(collection, "displayName", ""))

        if hasattr(collection, "description"):
            self.setToolTip(collection.description)

        attributes = [attr for attr in collection.elements("AttributeDescription")
                      if not is_hidden(attr)]

        numAttrs = max(len(attributes) / 2 + 1, 1)
        i = 0
        for attr in attributes:
            if attr.is_pointer():
                try:
                    dataset, newattr = attr.get_pointed()
                except ValueError as ex:
                    newattr = None

                if not newattr:
                    continue
                attr = newattr

            attr_widget = AttributeWidget(attr, self.dataset, self.master, self)
            self.layout().addWidget(attr_widget, i % numAttrs, i / numAttrs)
            self.addSubControl(attr, attr_widget)
            i += 1


def filterType(filter):
    return (getattr(filter, "displayType", ""),
            getattr(filter, "style", ""),
            getattr(filter, "multipleValues", ""))


class FilterCollectionWidget(CollectionWidget):

    def __init__(self, tree, dataset, master, parent=None):
        CollectionWidget.__init__(self, tree, dataset, master, parent)
        self.setLayout(QGridLayout())

        self.layout().setContentsMargins(5, 5, 5, 5)
        self.layout().setColumnStretch(0, 10)
        self.layout().setColumnStretch(1, 10)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self.setCollection(tree, dataset)

    def setCollection(self, collection, dataset, flags=0):
        CollectionWidget.setCollection(self, collection, dataset, flags)
        filters = [f for f in collection.elements("FilterDescription") if not is_hidden(f)] + \
                  [f for f in collection.elements("AttributeDescription") if not is_hidden(f)]  # in case of pointers to filters (upstream/downstream flank)
        self.enableCB = QCheckBox(getattr(collection, "displayName", ""), self)

        if hasattr(collection, "description"):
            self.setToolTip(collection.description)

        if len(filters) == 1 and (not hasattr(filters[0], "displayName") or \
                                  getattr(filters[0], "displayName", "") == getattr(collection, "displayName", "")):
            flags = flags | self.SINGLE_FILTER_FLAG

        i = 0
        if not flags & self.SINGLE_FILTER_FLAG:
            self.layout().addWidget(self.enableCB, 0, 0)
            i += 1

        for filter in filters:
            fType = getattr(filter, "type", None)
            if filter.is_pointer():
                try:
                    dataset, newfilter = filter.get_pointed()
                except ValueError as ex:
                    newfilter = None

                if not newfilter:
                    continue
                filter = newfilter
                fType = getattr(filter, "type", None) if fType is None else fType
            label, filter_widget = self.buildFilter(filter, flags, fTypeHint=fType)
            if isinstance(label, six.string_types):
                label = QLabel(label)

            self.layout().addWidget(label, i, 0)
            self.layout().addWidget(filter_widget, i, 1)
            i += 1

            self.addSubControl(filter, filter_widget)
        if self.layout().count() == 0:
            self.layout().addWidget(self.enableCB, 0, 0)

    def buildFilter(self, filter, flags, fTypeHint=None):
        if flags & self.SINGLE_FILTER_FLAG:
            label = self.enableCB
        else:
            label = getattr(filter, "displayName", "")

        fType = filterType(filter)
        if fTypeHint == "drop_down_basic_filter":
            fType = ("list", "menu", "")
        filter_widget = None

        if fType == ("text", "", ""):
            filter_widget = TextFieldFilter(filter, self.dataset, self.master, self)

        elif fType == ("text", "", "1") or fType == ("text", "", "true"):
            filter_widget = IdListFilter(filter, self.dataset, self.master, self)

        elif fType == ("list", "radio", ""):
            filter_widget = RadioBooleanFilter(filter, self.dataset, self.master, self)

        elif fType == ("list", "menu", ""):
            filter_widget = DropDownFilter(filter, self.dataset, self.master, self)

        elif fType == ("list", "menu", "1") or fType == ("list", "menu", "true"):
            filter_widget = MultiSelectListFilter(filter, self.dataset, self.master, self)

        elif fType == ("container", "", ""):
            fType = set(map(filterType, filter.elements_top("Option")))
            if len(fType) != 1:
                warnings.warn("Multiple filter types in a single container!" + str(fType))
            fType = fType.pop()
            if fType[0] == "text":  # ("text", "", "1"):
                filter_widget = DropDownIdListFilter(filter, self.dataset, self.master, self)
            elif fType == ("list", "radio", ""):
                filter_widget = DropDownRadioBooleanFilter(filter, self.dataset, self.master, self)

        if filter_widget is None:
            warnings.warn("Unknown filter type '%s' %s'" % (repr(fType), repr(filter)))
            filter_widget = UnknownFilter(filter, self.dataset, self.master, self)

        filter_widget.setMaximumWidth(400)
        return label, filter_widget

    def query(self):
        if self.enableCB.isChecked():
            return CollectionWidget.query(self)
        else:
            return []


class GroupWidget(QGroupBox, Control):

    def __init__(self, tree, dataset, master, parent=None):
        QGroupBox.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)
        self.setLayout(QVBoxLayout())

        if tree:
            self.setGroup(tree)

    def setGroup(self, group):
        self.setTitle(getattr(group, "displayName", ""))
        if hasattr(group, "description"):
            self.setToolTip(group.description)
        self.setCheckable(True)
        self.setChecked(True)

        if group.tag == "FilterGroup":
            collections = group.elements("FilterCollection")
        else:
            collections = group.elements("AttributeCollection")
        for collection in collections:
            if isFilterCollection(collection):
                collection_widget = FilterCollectionWidget(collection, self.dataset, self.master, self)
            else:
                collection_widget = AttributeCollectionWidget(collection, self.dataset, self.master, self)

            self.layout().addWidget(collection_widget)
            self.addSubControl(collection, collection_widget)


class PageWidget(QFrame, Control):

    def __init__(self, tree, dataset, master, parent=None):
        QFrame.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)
        self.setLayout(QVBoxLayout())
        self.outFormats = getattr(tree, "outFormats", "tsv")

        if tree:
            self.setPage(tree)

    def setPage(self, tree):
        if tree.tag == "FilterPage":
            groups = tree.elements("FilterGroup")
        else:
            groups = tree.elements("AttributeGroup")

        for group in groups:
            group_widget = GroupWidget(group, self.dataset, self.master, self)
            self.layout().addWidget(group_widget)
            self.addSubControl(group, group_widget)

        self.layout().addStretch(10)

        if hasattr(tree, "description"):
            self.setToolTip(tree.description)

# A list of preset mart services available for selection from the pull
# down menu

MartServices = [
    ("Ensembl", "http://www.ensembl.org/biomart/martservice"),
    ("Ensembl Fungi", "http://fungi.ensembl.org/biomart/martservice"),
    ("Ensembl Plants", "http://plants.ensembl.org/biomart/martservice"),
    ("Pancreatic Expression Database",
     "http://www.pancreasexpression.org/biomart/martservice"),
    ("VectorBase", "http://biomart.vectorbase.org/biomart/martservice"),
    ("Phytozome", "https://phytozome.jgi.doe.gov/biomart/martservice"),
]


class OWBioMart(widget.OWWidget):
    name = "BioMart"
    description = "Query BioMart service"
    icon = "../widgets/icons/BioMart.svg"
    priority = 2010

    outputs = [("Data", Orange.data.Table)]

    SHOW_FILTERS = True

    selectedService = settings.Setting(MartServices[0][1])
    selectedDataset = settings.Setting(0)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.selectedDatabase = 0
        self.uniqueRows = True

        gui.button(gui.widgetBox(self.controlArea, "Cache", addSpace=True),
                   self, "Clear cache",
                   tooltip="Clear saved query results",
                   callback=self.clearCache)
        self.serviceindex = 0
        self.serviceCombo = gui.comboBox(
            self.controlArea, self, "serviceindex", "Mart Service",
            callback=self._setServiceUrl
        )
        for name, url in MartServices:
            self.serviceCombo.addItem(name, userData=url)
        idx = self.serviceCombo.findData(self.selectedService, Qt.UserRole)
        self.serviceCombo.setCurrentIndex(idx)
        # self.selectedService = self.serviceCombo.itemData(self.serviceCombo.currentItem())

        self.martsCombo = gui.comboBox(
            self.controlArea, self, "selectedDatabase", "Database",
            callback=self.setSelectedMart,
            addSpace=True)
        self.martsCombo.setMaximumWidth(250)

        self.datasetsCombo = gui.comboBox(
            self.controlArea, self, "selectedDataset", "Dataset",
            callback=self.setSelectedDataset,
            addSpace=True)

        self.datasetsCombo.setMaximumWidth(250)

        gui.rubber(self.controlArea)

        box = gui.widgetBox(self.controlArea, "Results")
        gui.checkBox(
            box, self, "uniqueRows", "Unique results only",
            tooltip="Return unique results only.",)

        self.commitButton = gui.button(
            box, self, "Get Results", callback=self.commit,
            tooltip="Query the BioMart server and output the results",
            autoDefault=True)

        self.commitButton.setEnabled(False)

        self.mainWidget = gui.widgetBox(
            self.mainArea, orientation=QStackedLayout())

        self.mainTab = QTabWidget()

        self.mainWidget.layout().addWidget(self.mainTab)

        self.attributesConfigurationBox = gui.createTabPage(self.mainTab, "Attributes")

        if self.SHOW_FILTERS:  # ??
            self.filtersConfigurationBox = gui.createTabPage(self.mainTab, "Filters")

        self.error(0)
        self.setEnabled(False)
        self._task = None
        self._executor = concurrent.ThreadExecutor(
            threadPool=QThreadPool(maxThreadCount=2)
        )
        service = self.selectedService
        self._task = task = concurrent.Task(
            function=partial(self._get_registry, url=service))
        task.resultReady.connect(self.setBioMartRegistry)
        task.exceptionReady.connect(self._handleException)
        self._executor.submit(task)
        self._setServiceUrl()
        self._afterInitQueue = []

        try:
            from Bio import SeqIO
            self.hasBiopython = True
        except ImportError:
            self.warning(100, "Biopython package not found.\nTo retrieve FASTA sequence data from BioMart install Biopython.")
            self.hasBiopython = False

    def sizeHint(self):
        return QSize(800, 600)

    def _setServiceUrl(self):
        service = self.serviceCombo.itemData(self.serviceCombo.currentIndex())
        if service is not None:
            self.selectedService = service
            self._task = task = concurrent.Task(
                function=partial(self._get_registry, url=service))
            task.resultReady.connect(self.setBioMartRegistry)
            task.exceptionReady.connect(self._handleException)
            self._executor.submit(task)

    @staticmethod
    def _get_registry(url=None, precache=True):
        if url is None:
            url = MartServices[0][1]
        con = biomart.BioMartConnection(address=url, timeout=30)
        reg = biomart.BioMartRegistry(con)
        if precache:
            _ = reg.marts()
        return reg

    @Slot(Exception)
    def _handleException(self, exception):
        assert(QThread.currentThread() is self.thread())
        print("Task failed with:", exception, file=sys.stderr)
        import logging
        log = logging.getLogger(__name__)
        log.exception("Error:", exc_info=exception)
        self.error(0, str(exception))
        self.setEnabled(True)

    @Slot(object)
    def setBioMartRegistry(self, registry):
        assert(QThread.currentThread() is self.thread())
        self.setEnabled(True)
        self.registry = registry
        self.marts = [mart for mart in self.registry.marts()
                      if getattr(mart, "visible", "0") != "0"]

        self.martsCombo.clear()
        for mart in self.marts:
            self.martsCombo.addItem(mart.displayName)

    def setSelectedMart(self):
        self.mart = self.marts[self.selectedDatabase]
        self.error(0)
        self.setEnabled(False)

        self._task = task = concurrent.Task(function=self.mart.datasets)
        task.resultReady.connect(self.setBioMartDatasets)
        task.exceptionReady.connect(self._handleException)
        self._executor.submit(task)

    @Slot(object)
    def setBioMartDatasets(self, datasets):
        assert(QThread.currentThread() is self.thread())
        self.setEnabled(True)
        self.datasets = [data for data in datasets if
                         getattr(data, "visible", "0") != "0"]
        self.datasetsCombo.clear()
        self.datasetsCombo.addItems([data.displayName for data in self.datasets])

    def setSelectedDataset(self):
        self.dataset = self.datasets[self.selectedDataset]
        self.error(0)
        self.setEnabled(False)

        def get_configuration(dataset):
            connection = dataset.connection
            stream = connection.configuration(
                dataset=dataset.internalName,
                virtualSchema=dataset.serverVirtualSchema)
            response = stream.read()
            return response

        self._task = task = concurrent.Task(
            function=partial(get_configuration, self.dataset))

        task.resultReady.connect(self.setBioMartConfiguration)
        task.exceptionReady.connect(self._handleException)

        self._executor.submit(task)

    @Slot(object)
    def setBioMartConfiguration(self, configuration):
        assert(QThread.currentThread() is self.thread())
        self.setEnabled(True)
        # parse the xml in the main thread (a long time ago this step was
        # done in a thread but would frequently cause `expat` to segfault.
        doc = biomart.parseXML(io.BytesIO(configuration))
        config = list(doc.elements("DatasetConfig"))[0]
        configuration = biomart.DatasetConfig(
            self.registry, config.tag, config.attributes, config.children)

        self.clearConfiguration()

        self.configuration = configuration

        def hidden(tree):
            return getattr(tree, "hidden", "false") != "false" or \
                   getattr(tree, "hideDisplay", "false") != "false"

        self.attributePagesTabWidget = tabs = gui.tabWidget(self.attributesConfigurationBox)

        for page in configuration.elements("AttributePage"):
            if not hidden(page):
                page_widget = PageWidget(page, self.dataset, self)
                gui.createTabPage(tabs, getattr(page, "displayName", ""),
                                  widgetToAdd=page_widget, canScroll=True)

        if self.SHOW_FILTERS:
            self.filterPagesTabWidget = tabs = gui.tabWidget(self.filtersConfigurationBox)
            for page in configuration.elements("FilterPage"):
                if not hidden(page):
                    page_widget = PageWidget(page, self.dataset, self)
                    gui.createTabPage(tabs, getattr(page, "displayName", ""),
                                      widgetToAdd=page_widget, canScroll=True)

        self.afterInit()

        self.commitButton.setEnabled(True)

    def clearConfiguration(self):
        self.mainTab.deleteLater()

        self.mainTab = QTabWidget()
        self.mainWidget.layout().addWidget(self.mainTab)
        self.mainWidget.layout().setCurrentWidget(self.mainTab)

        self.attributesConfigurationBox = gui.createTabPage(self.mainTab, "Attributes")
        if self.SHOW_FILTERS:
            self.filtersConfigurationBox = gui.createTabPage(self.mainTab, "Filters")

    def commit(self):
        pageconf = self.attributePagesTabWidget.currentWidget().widget()
        format = pageconf.outFormats

        self.error(100)
        if not self.hasBiopython and format.lower() == "fasta":
            self.error(100, "Cannot parse FASTA format")
            return

        query = pageconf.query()
        bydatasets = defaultdict(lambda: ([], []))

        for conftype, tree, val in query:
            dataset = self.dataset

            if conftype == "Attribute":
                bydatasets[dataset][0].append(tree.internalName)
            elif conftype == "Filter":
                bydatasets[dataset][1].append((tree.internalName, val))

        if self.SHOW_FILTERS:
            pageconf = self.filterPagesTabWidget.currentWidget().widget()
            query = pageconf.query()

            for conftype, tree, val in query:
                dataset = self.dataset

                if conftype == "Attribute":
                    bydatasets[dataset][0].append(tree.internalName)
                elif conftype == "Filter":
                    bydatasets[dataset][1].append((tree.internalName, val))

        query = self.registry.query(
            format="TSV" if "tsv" in format.lower() else format.upper(),
            uniqueRows=self.uniqueRows,
            virtualSchema=dataset.virtualSchema,
            serverVirtualSchema=dataset.serverVirtualSchema
        )

        for dataset, (attributes, filters) in bydatasets.items():
            query.set_dataset(dataset if dataset else self.dataset)
            for attr in attributes:
                query.add_attribute(attr)
            for filter, value in filters:
                query.add_filter(filter, value)

        self.error(0)
        self.setEnabled(False)
        self._task = task = concurrent.Task(function=query.get_table)
        task.resultReady.connect(self.dataReady)
        task.exceptionReady.connect(self._handleException)
        self._executor.submit(task)

    def dataReady(self, data):
        self.setEnabled(True)
        self.send("Data", data)

    def pushAction(self, action):
        ref = action.ref
        ref_widget = self.findChild(QWidget, ref)
        if hasattr(ref_widget, "setOptions"):
            ref_widget.setOptions(action.subelements_top("Option"))

    def registerDelayedCall(self, call):
        self._afterInitQueue.append(call)

    def afterInit(self):
        while self._afterInitQueue:
            call = self._afterInitQueue.pop(0)
            call()

    def clearCache(self):
        self.registry.connection.clear_cache()


def test_main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv))

    ow = OWBioMart()
    ow.show()
    ow.raise_()

    r = app.exec_()
    ow.saveSettings()
    return r

if __name__ == "__main__":
    sys.exit(test_main())
