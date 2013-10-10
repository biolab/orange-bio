"""<name>BioMart</name>
<description>Query BioMart service</description>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si</contact>
<priority>2010</priority>
<icon>icons/BioMart.svg</icon>
"""

from __future__ import absolute_import

import sys, os
import traceback
import warnings
import socket
from collections import defaultdict
import itertools

import Orange

from Orange.OrangeWidgets import OWConcurrent
from Orange.OrangeWidgets.OWWidget import *

from .. import obiBioMart
from ..obiBioMart import *

NAME = "BioMart"
DESCRIPTION = "Query BioMart service"
ICON = "icons/BioMart.svg"
PRIORITY = 2010

INPUTS = [("Input ids", Orange.data.Table, "setData")]
OUTPUTS = [("Example Table", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWBioMart.OWBioMart"]


socket.setdefaulttimeout(60)


def is_hidden(tree):
    return getattr(tree, "hidden", "false") != "false" or getattr(tree, "hideDisplay", "false") != "false"

class Control(object):
    """ Base mixin class for query GUI widgets 
    """
    def __init__(self, tree, dataset, master):
        """ 
            :param tree: configuration element
            :type tree: obiBioMart.ConfigurationNode
            :param dataset: dataset
            :type dataset: obiBioMart.BioMartDataset
            :param master: main widget
            :type master: OWBioMart
            
        """
        self.tree = tree
        self.dataset = dataset
        self.master = master
        self.subControls = []
        
        if isinstance(self, QObject):
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
            controls = dict([(tree.internalName, control) for tree, control in self.subControls])
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
        self.textWidget = QPlainTextEdit() #TODO: Subclass to recieve drop events from item model views
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
                
        self.connect(self, SIGNAL("currentIndexChanged(int)"), self.onIndexChange)
        
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
                
    
#class GraphDropDownFilter(QComboBox, Control):
#    """ Tree menu filter
#    """
#    def __init__(self, tree, dataset, master, parent=None):
#        QComboBox.__init__(self, parent)
#        Control.__init__(self, tree, dataset, master)
#        
#        self.options = []
#        self.selectedIndex = 0
#
#        self.installEventFilter(self)
#        
#        self.setOptions(tree.subelements_top("Option"))
#                
##        self.connect(self, SIGNAL("currentIndexChanged(int)"), self.onIndexChange)
#        
#    def setOptions(self, options):
#        self.options = list(options)
#        self.clear()
#        self.addItem(self.options[0].displayName)
#        self.value = self.options[0].value
#        
#    def eventFilter(self, obj, event):
#        if event.type() == QEvent.MouseButtonPress:
#            self.showMenu()
#            return True
#        elif event.type() == QEvent.MouseButtonRelease:
#            return True
#        return False
#    
#    def showPopup(self):
#        return
#    
#    def closePopup(self):
#        return
#    
#    def buildMenu(self):
#        menu = QMenu(self)
#        
#        def addOption(menu, option):
#            suboptions = list(option.subelements_top("Option"))
#            if suboptions:
#                submenu = menu.addMenu(option.displayName)
#                self.connect(submenu, SIGNAL("hovered(QAction)"), self.updateLastHovered)
#                for op in suboptions:
#                    addOption(submenu, op)
#            else:
#                action = menu.addAction(option.displayName)
#                action._value = option.value
#                
#        for option in self.options:
#            addOption(menu, option)
#        
#        self.connect(menu, SIGNAL("hovered(QAction)"), self.updateLastHovered)
#        self.connect(menu, SIGNAL("triggered(QAction)"), lambda action: menu.close())
#        return menu 
#    
#    def updateLastHovered(self, action):
#        self.lastHovered = str(action.text())
#        print self.lastHovered
#    def showMenu(self):
#        menu = self.buildMenu()
#        action = menu.exec_(self.mapToGlobal(QPoint(0, 0)))
#        name = str(action.text())
#        self._value = value = action._value
#        self.setItemText(0, name)
#                    
#    def query(self):
#        return [("Filter", self.tree, self._value)]
    
    
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
        self.model().setStringList([option.displayName for option in self.options])
                
    def value(self):
        value = []
        for index in rows(self.selectedIndexes()):
            value.append(self.options[index].value)
        return ",".join(value)
    
    def query(self):
        return [("Filter", self.tree, self.value())]
    
    def setControlValue(self, name, value):
        if type(value) == str:
            value = value.split(",")
        selection = self.selectionModel()
        for i, option in enumerate(self.options):
            if option.value in values:
                selection.select(self.model().index(i), QItemSelectionModel.Select)
            
    
class DropDownRadioBooleanFilter(QWidget, Control):
    """ Container for multiple boolean filters
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
        filter = FilterDescription(self.tree.registry, "FilterDescription", filter.attributes, filter.children)
        return [("Filter", filter, self.value())] #"only" if self.radioButtons[0].isChecked() else "excluded")]
    
    def setControlValue(self, name, value):
        for i, option in enumerate(self.options):
            if option.internalName == name:
                self.cb.setCurrentIndex(i)
                if value == "Only":
                    self.radioButtons[0].setChecked(True)


class DropDownIdListFilter(QWidget, Control):
    """ Container for multiple id list filters
    """
    def __init__(self, tree, dataset, master, parent=None):
        QWidget.__init__(self, parent)
        Control.__init__(self, tree, dataset, master)
        
        self.setLayout(QVBoxLayout())
        self.setContentsMargins(0, 0, 0, 0)
        self.cb = QComboBox()
        self.idsEdit = QPlainTextEdit()
#        self.browseButton = QPushButton("Browse ...")
        
        self.layout().addWidget(self.cb)
        self.layout().addWidget(self.idsEdit)
#        self.layout().addWidget(self.browseButton)
#        self.layout().setAlignment(self.browseButton, Qt.AlignRight)
        
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
        filter = FilterDescription(self.tree.registry, "FilterDescription", filter.attributes, filter.children)
        return [("Filter", filter, self.value())]
    
    def setControlValue(self, name, value):
        if type(value) == list:
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
#        self.setLayout(QFormLayout())
#        self.enableCB = QCheckBox(getattr(tree, "displayName", ""), self)
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
    return all(["pointerFilter" in desc.attributes  for desc in tree.elements("AttributeDescription")])


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
            
        attributes = [attr for attr in collection.elements("AttributeDescription") if not is_hidden(attr)]
        
        numAttrs = max(len(attributes) / 2 + 1, 1)
        i = 0
        for attr in attributes:
            if attr.is_pointer():
                try:
                    dataset, newattr = attr.get_pointed()
                except ValueError, ex:
                    newattr = None
                    
                if not newattr:
                    continue
                attr = newattr
                
            attr_widget = AttributeWidget(attr, self.dataset, self.master, self)
            self.layout().addWidget(attr_widget, i % numAttrs, i / numAttrs)
            self.addSubControl(attr, attr_widget)
            i += 1
            
        
def filterType(filter):
    return getattr(filter, "displayType", ""), getattr(filter, "style", ""), getattr(filter, "multipleValues", "")

 
class FilterCollectionWidget(CollectionWidget):
    def __init__(self, tree, dataset, master, parent=None):
        CollectionWidget.__init__(self, tree, dataset, master, parent)
        self.setLayout(QGridLayout())
#        self.setLayout(QFormLayout())
#        self.layout().setLabelAlignment(Qt.AlignLeft | Qt.AlignVCenter)
#        self.layout().setFormAlignment(Qt.AlignHCenter | Qt.AlignTop)
        
        self.layout().setContentsMargins(5, 5, 5, 5)
        self.layout().setColumnStretch(0, 10)
        self.layout().setColumnStretch(1, 10)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self.setCollection(tree, dataset)
        
    def setCollection(self, collection, dataset, flags=0):
        CollectionWidget.setCollection(self, collection, dataset, flags)
        filters = [f for f in collection.elements("FilterDescription") if not is_hidden(f)] + \
                  [f for f in collection.elements("AttributeDescription") if not is_hidden(f)] # in case of pointers to filters (upstream/downstream flank)
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
#            self.layout().addRow(self.enableCB, QWidget(self))
            
        for filter in filters:
            fType = getattr(filter, "type", None)
            if filter.is_pointer():
                try:
                    dataset, newfilter  = filter.get_pointed()
                except ValueError, ex:
                    newfilter = None
                    
                if not newfilter:
                    continue
                filter = newfilter
                fType = getattr(filter, "type", None) if fType is None else fType
            label, filter_widget = self.buildFilter(filter, flags, fTypeHint=fType)
            if isinstance(label, basestring):
                label = QLabel(label) 
            
            self.layout().addWidget(label,i, 0)
            self.layout().addWidget(filter_widget, i, 1)
            i += 1
#            self.layout().addRow(label, filter_widget)
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
            
        elif fType == ("text", "", "1"):
            filter_widget = IdListFilter(filter, self.dataset, self.master, self)
            
        elif fType == ("list", "radio", ""):
            filter_widget = RadioBooleanFilter(filter, self.dataset, self.master, self)
            
        elif fType == ("list", "menu", ""):
#            if getattr(filter, "graph", "0") == "1":
#                filter_widget = GraphDropDownFilter(filter, self.dataset, self.master, self)
#            else:
#                filter_widget = DropDownFilter(filter, self.dataset, self.master, self)
            filter_widget = DropDownFilter(filter, self.dataset, self.master, self)
            
        elif fType == ("list", "menu", "1"):
            filter_widget = MultiSelectListFilter(filter, self.dataset, self.master, self) 
            
        elif fType == ("container", "", ""):
            fType = set(map(filterType, filter.elements_top("Option")))
            if len(fType) != 1:
                warnings.warn("Multiple filter types in a single container!" + str(fType))
            fType = fType.pop()
#            if fType == ("text", "", ""):
#                filter_widget = DropDownTextFilter(filter, self.dataset, self.master, self)
            if fType[0] == "text": #("text", "", "1"):
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
        

    
class OWBioMart(OWWidget):
    SHOW_FILTERS = True
    
    settingsList = ["selectedDataset"]
    def __init__(self, parent=None, signalManager=None, title="BioMart"):
        OWWidget.__init__(self, parent, signalManager, title)
        
        self.inputs = [("Input ids", ExampleTable, self.setData)]
        self.outputs = [("Example Table", ExampleTable)]
        
        self.selectedDatabase = None
        self.selectedDataset = None
        self.lastSelectedDatasets = (None, None)
        
        self.idAttr = 0
        self.useAttrNames = False
        self.uniqueRows = True
        
        self.loadSettings()
        
        OWGUI.button(OWGUI.widgetBox(self.controlArea, "Cache", addSpace=True),
                     self, "Clear cache",
                     tooltip="Clear saved query results",
                     callback = self.clearCache)
        
        
        self.martsCombo = OWGUI.comboBox(self.controlArea, self, "selectedDatabase", "Database",
                                         callback=self.setSelectedMart,
                                         addSpace=True)
        self.martsCombo.setMaximumWidth(250)
        
        self.datasetsCombo = OWGUI.comboBox(self.controlArea, self, "selectedDataset", "Dataset",
                                            callback=self.setSelectedDataset,
                                            addSpace=True)
        self.datasetsCombo.setMaximumWidth(250)
        
        box = OWGUI.widgetBox(self.controlArea, "Input Ids")
        
        cb = OWGUI.checkBox(box, self, "useAttrNames", "Use attribute names", 
                            tooltip="Use attribute names for ids",
                            callback=self.updateInputIds)
        
        self.idAttrCB = OWGUI.comboBox(box, self, "idAttr", 
                                       tooltip="Use attribute values from ...",
                                       callback=self.updateInputIds)
        
        cb.disables = [(-1, self.idAttrCB)]
        cb.makeConsistent()
        
#        self.idListView = QListView()
#        self.idListView.setSelectionMode(QListView.ExtendedSelection)
#        self.idListView.setModel(QStringListModel(self))
#        box.layout().addWidget(self.idListView)

        self.idText = QPlainTextEdit()
        self.idText.setReadOnly(True)
#        self.idText.setTextInteractionFlags(Qt.TextSelectableByKeyboard)
        box.layout().addWidget(self.idText)
        
        OWGUI.rubber(self.controlArea)
        
        box = OWGUI.widgetBox(self.controlArea, "Results")
        OWGUI.checkBox(box, self, "uniqueRows", "Unique results only",
                       tooltip="Return unique results only.",
                       )
        self.commitButton = OWGUI.button(box, self, "Get Results",
                                         callback=self.commit,
                                         tooltip="Query the BioMart server for results and output the results",
                                         autoDefault=True)
        self.commitButton.setEnabled(False)
        
        self.mainWidget = OWGUI.widgetBox(self.mainArea, orientation=QStackedLayout())
        
        self.mainTab = QTabWidget()
        
        self.mainWidget.layout().addWidget(self.mainTab)
        
        self.attributesConfigurationBox = OWGUI.createTabPage(self.mainTab, "Attributes")#, canScroll=True)
        if self.SHOW_FILTERS:
            self.filtersConfigurationBox = OWGUI.createTabPage(self.mainTab, "Filters")#, canScroll=True)

        self.myThread = OWConcurrent.WorkerThread()
        self.myThread.start()
        
        self.__thread = self.thread() # For assert(QThread.currentThread() is self.thread()) to work
        
        self.connect(self.myThread, SIGNAL("started()"), lambda :sys.stderr.write("Thread started\n"))
        self.connect(self.myThread, SIGNAL("finished()"), lambda :sys.stderr.write("Thread finished\n"))
        
        self.error(0)
        self.setEnabled(False)
        self.get_registry_async = OWConcurrent.createTask(self._get_registry,
                                      onResult=self.setBioMartRegistry,
                                      onFinished=self.onFinished,
                                      thread=self.myThread)
        
        
        self.resize(800, 600)
        
        self.candidateIdAttrs = []
        self._afterInitQueue = []
        self.data = None
        
        try:
            from Bio import SeqIO
            self.hasBiopython = True
        except Exception, ex:
            self.warning(100, "Biopython package not found.\nTo retrieve FASTA sequence data from BioMart install Biopython.")
            self.hasBiopython = False
        
        
    @staticmethod
    def _get_registry(url=None, precache=True):
        con = obiBioMart.BioMartConnection(url)
        reg = obiBioMart.BioMartRegistry(con)
        if precache:
            marts = reg.marts()
        return reg
        
    @pyqtSignature("onFinished(QString)")
    def onFinished(self, status):
        assert(QThread.currentThread() is self.thread())
        if str(status).lower() != "ok":
            print >> sys.stderr, "AsyncCall failed with message:", status
            self.error(0, str(status))
        self.setEnabled(True)
            
    @pyqtSignature("setBioMartRegistry(PyQt_PyObject)")
    def setBioMartRegistry(self, registry):
        assert(QThread.currentThread() is self.thread())
        self.registry = registry
        self.marts = [mart for mart in self.registry.marts() if getattr(mart, "visible", "0") != "0"]
        last, _ = self.lastSelectedDatasets 
        for mart in self.marts:
            self.martsCombo.addItem(mart.displayName)
            
    def setSelectedMart(self):
        self.mart = self.marts[self.selectedDatabase]
        self.error(0)
        self.setEnabled(False)
        self.get_datasets_async = OWConcurrent.createTask(self.mart.datasets,
                                    onResult=self.setBioMartDatasets,
                                    onFinished=self.onFinished,
                                    thread=self.myThread)
        

    @pyqtSignature("setBioMartDatasets(PyQt_PyObject)")
    def setBioMartDatasets(self, datasets):
        assert(QThread.currentThread() is self.thread())
        self.datasets = [data for data in datasets if getattr(data,"visible", "0") != "0"]
        self.datasetsCombo.clear()
        self.datasetsCombo.addItems([data.displayName for data in self.datasets])
            
    def setSelectedDataset(self):
        self.dataset = self.datasets[self.selectedDataset]
        self.error(0)
        
        def get_configuration(dataset):
            """ Only get the response, do not yet parse
            """
            connection = dataset.connection
            stream = connection.configuration(dataset=dataset.internalName, virtualSchema=dataset.virtualSchema)
            response = stream.read()
            return response
        
        self.setEnabled(False)
        
        self.get_configuration_async = OWConcurrent.createTask(get_configuration, (self.dataset,),
                                            onResult=self.setBioMartConfiguration,
                                            onFinished=self.onFinished,
                                            thread=self.myThread)
        
        
    @pyqtSignature("setBioMartConfiguration(PyQt_PyObject)")
    def setBioMartConfiguration(self, configuration):
        assert(QThread.currentThread() is self.thread())
        
        doc = obiBioMart.parseXML(configuration)
        config = list(doc.elements("DatasetConfig"))[0]
        configuration = obiBioMart.DatasetConfig(self.registry, config.tag, config.attributes, config.children)
        
        self.clearConfiguration()
        
        self.configuration = configuration
        
        def hidden(tree):
            return getattr(tree, "hidden", "false") != "false" or getattr(tree, "hideDisplay", "false") != "false"
        
        self.attributePagesTabWidget = tabs =  OWGUI.tabWidget(self.attributesConfigurationBox)
        
        for page in configuration.elements("AttributePage"):
            if not hidden(page):
                page_widget = PageWidget(page, self.dataset, self)
                OWGUI.createTabPage(tabs, getattr(page, "displayName", ""), widgetToAdd=page_widget, canScroll=True)
        
        if self.SHOW_FILTERS:
            self.filterPagesTabWidget = tabs = OWGUI.tabWidget(self.filtersConfigurationBox)
            for page in configuration.elements("FilterPage"):
                if not hidden(page):
                    page_widget = PageWidget(page, self.dataset, self)
                    OWGUI.createTabPage(tabs, getattr(page, "displayName", ""), widgetToAdd=page_widget, canScroll=True)
                    
        self.afterInit()
        
        self.commitButton.setEnabled(True)
    
    def clearConfiguration_(self):
        while True:
            w = self.attributesConfigurationBox.layout().takeAt(0)
            if isinstance(w, QWidgetItem):
                w = w.widget()
            if not w:
                break
            if isinstance(w, QWidget):
                w.setParent(None)
        
        while True:
            w = self.filtersConfigurationBox.layout().takeAt(0)
            if isinstance(w, QWidgetItem):
                w = w.widget()
            if not w:
                break
            if isinstance(w, QWidget):
                w.setParent(None)
            pass
        
    def clearConfiguration(self):
        self.mainTab = QTabWidget()
        self.mainWidget.layout().addWidget(self.mainTab)
        self.mainWidget.layout().setCurrentWidget(self.mainTab)
        
        self.attributesConfigurationBox = OWGUI.createTabPage(self.mainTab, "Attributes")#, canScroll=True)
        if self.SHOW_FILTERS:
            self.filtersConfigurationBox = OWGUI.createTabPage(self.mainTab, "Filters")#, canScroll=True)
            
            
    def count(self):
        filters = []
        count = self.dataset.get_count(filters = filters)
        
        
    def commit(self):
        pageconf = self.attributePagesTabWidget.currentWidget().widget()
        format = pageconf.outFormats
        
        self.error(100)
        if not self.hasBiopython and format.lower() == "fasta":
            self.error(100, "Cannot parse FASTA format")
            return
        
        query = pageconf.query()
        bydatasets = defaultdict(lambda : ([], []))

        version = getattr(self.configuration, "softwareVersion", "0.4")
        for conftype, tree, val in query:
            dataset = self.dataset
#            if version > "0.4":
#                dataset = self.dataset
#            else:
#                widget = self.findChild(QWidget, tree.internalName)
#                if widget:
#                    dataset = widget.dataset
                
            if conftype == "Attribute":
                bydatasets[dataset][0].append(tree.internalName)
            elif conftype == "Filter":
                bydatasets[dataset][1].append((tree.internalName, val))
                
        if self.SHOW_FILTERS:
            pageconf = self.filterPagesTabWidget.currentWidget().widget()
#            configuration =  pageconf.get_configuration()
            query = pageconf.query()
    
            for conftype, tree, val in query:
                dataset = self.dataset
#                if version > "0.4":
#                    dataset = self.dataset
#                else:
#                    widget = self.findChild(QWidget, tree.internalName)
#                    if widget:
#                        dataset = widget.dataset
                    
                if conftype == "Attribute":
                    bydatasets[dataset][0].append(tree.internalName)
                elif conftype == "Filter":
                    bydatasets[dataset][1].append((tree.internalName, val))
                
        query = self.registry.query(format="TSV" if "tsv" in format.lower() else format.upper(), uniqueRows=self.uniqueRows)
        
        for dataset, (attributes, filters) in bydatasets.items():
            query.set_dataset(dataset if dataset else self.dataset)
            for attr in attributes:
                query.add_attribute(attr)
            for filter, value in filters:
                query.add_filter(filter, value)
        
        self.error(0)
        self.setEnabled(False)
        self.run_query_async = OWConcurrent.createTask(query.get_example_table,
                                            onResult=self.dataReady,
                                            onFinished=self.onFinished,
                                            thread=self.myThread)

    
    def dataReady(self, data):
        self.setEnabled(True)
        self.send("Example Table", data)
        
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
            
    def setData(self, data =None):
        self.data = data
        self.idAttrCB.clear()
        attrs = data.domain.variables + data.domain.getmetas().values()
        attrs = [attr for attr in attrs if isinstance(attr, orange.StringVariable)]
        self.candidateIdAttrs = attrs
        self.idAttrCB.addItems([attr.name for attr in attrs])
        
        self.idAttr = min(self.candidateIdAttrs, len(self.candidateIdAttrs) - 1)
        
        self.updateInputIds()
    
    def updateInputIds(self):
        if self.data:
            if self.useAttrNames:
                names = [attr.name for attr in self.data.domain.attributes]
            elif self.candidateIdAttrs:
                attr = self.candidateIdAttrs[self.idAttr]
                names = [str(ex[attr]) for ex in self.data if not ex[attr].isSpecial()]
            else:
                names = []
        else:
            names = []
            
#        self.idListView.model().setStringList(names)
        self.idText.setPlainText("\n".join(names))
        
        
    def clearCache(self):
        self.registry.connection.clearCache()
    
if __name__ == "__main__":
    app = QApplication(sys.argv)
    ow = OWBioMart()
    ow.show()
    
    data = orange.ExampleTable("../../../doc/datasets/brown-selected")
    ow.setData(data)
    app.exec_()
    ow.saveSettings()
