"""<name>OWBioMart</name>
"""

from OWWidget import *
from PyQt4.QtNetwork import *
import obiBioMart

from obiBioMart import *

import sys, os
import traceback

from collections import defaultdict

#class QBioMartConnection(obiBioMart.BioMartConnection):
#    def __init__(self, url=None, manager=None):
#        obiBioMart.BioMartConnection.__init__(self, url)
#        self.manager = manager if manager else QNetworkAccessManager()
#        QObject.connect(self.manager, SIGNAL("finished(QNetworkReply)"), self.finished)
#        reply = self.manager.get(QNetworkRequest(QUrl("http://www.trolltech.com")))
#        self._f = False
#        
#    def finished(self, reply):
#        self._f = True
#        
#    def request(self, **kwargs):
#        print self.request_url(**kwargs)
#        reply = self.manager.get(QNetworkRequest(QUrl(self.request_url(**kwargs))))
#        QObject.connect(reply, SIGNAL("error(QNetworkReply::NetworkError)"), lambda error: sys.stdout.write(error))
#        QObject.connect(reply, SIGNAL("downloadProgress(qint64, qint64)"), lambda r, l: sys.stdout.write(str(r)))
#        data = ""
#        print reply.atEnd(), reply.error()
#        while not self._f:
#            qApp.processEvents()
#            reply.waitForReadyRead(100)
#            data += str(reply.readAll())
#        print data
#        import StringIO
#        return StringIO.StringIO(data)

def cached_function_instance(func):
    from functools import wraps
    cache = {}
    @wraps(func)
    def f(*args, **kwargs):
        sig = args + tuple(sorted(kwargs.items()))
        if sig not in cache:
            cache[sig] = func(*args, **kwargs)
        return cache[sig]
    return f

class BackgroundTask(QThread):
    def __init__(self, parent=None):
        QThread.__init__(self, parent)
        self.connect(self, SIGNAL("runTask()"), self.executeTask, Qt.QueuedConnection)
        self.start()
        
    def run(self):
        self.exec_()
        
    def runTask(self, callable):
        self.emit(SIGNAL("runTask()"), callable)
    
    @pyqtSignature("executeTask(PyQt_PyObject)")
    def executeTask(self, callable):
        self.usleep(10)
        try:
            results = callable()
        except Exception, ex:
            print >> sys.stderr, "Exeption occured in BackgroundTask.executeTask, running", callable
            traceback.print_exc(file=sys.stderr)
            self.emit(SIGNAL("executionFailed(PyQt_PyObject)"), ex)
            return
        self.emit(SIGNAL("executionFinished(PyQt_PyObject)"), results)
        
class RunnableWrapper(QRunnable):
    def __init__(self, callable):
        QRunnable.__init__(self)
        self.setAutoDelete(False)
        self.callable = callable
        
    def run(self):
        try:
            self.results = self.callable()
        except Exception, ex:
            print >> sys.stderr, ex
            
class BackgroundTask(QObject):
    def __init__(self, callable, parent=None):
        QObject.__init__(self, parent)
        self.callable = callable
        self.connect(self, SIGNAL("_ow_start()"), self.execute, Qt.QueuedConnection)
        
    @pyqtSignature("execute()")
    def execute(self):
        try:
            self.results = self.callable()
        except Exception, ex:
            import traceback
            traceback.print_exc()
            print >> sys.stderr, "Exception occured in BackgroundTask.execute, running", self.callable
            self.emit(SIGNAL("finished(QString)"), str(ex))
            return
        
        self.emit(SIGNAL("finished(QString)"), self.results)
        
    def __call__(self):
        print >> sys.stderr, QThread.currentThread()
        self.emit(SIGNAL("_ow_start()"))
        
    def run(self):
        self._runnable = RunnableWrapper(self)
        QThreadPool.globalInstance().reserveThread()
        QThreadPool.globalInstance().start(self._runnable)
        print QThreadPool.globalInstance().activeThreadCount()
        
class RunnableTask(QRunnable):
    def __init__(self, call):
        QRunnable.__init__(self)
        self.setAutoDelete(False)
        self._call = call
        
    def run(self):
        self._return = self._call()
        
class MartserviceTask(QObject):
    def __init__(self, callable, parent=None):
        QObject.__init__(self, parent)
        self.callable = callable
        
    def __call__(self):
        try:
            print >> sys.stderr, "Calling:", self.callable
            self.result = self.callable()
        except Exception, ex:
            self.emit(SIGNAL("finished(QString)"), QString(str(ex)))
            del self.callable
            print >> sys.stderr, ex
            return
        self.emit(SIGNAL("finished(QString)"), QString("Ok"))
        del self.callable
        
        print >> sys.stderr, self.result
        return self.result
    
    def run(self):
        QThreadPool.globalInstance().setExpiryTimeout(-1) #
        self._runnable = RunnableWrapper(self)
        QThreadPool.globalInstance().reserveThread()
        print >> sys.stderr, "Scheduling", self.callable
        QThreadPool.globalInstance().start(self._runnable)
        print "Active pool size:", QThreadPool.globalInstance().activeThreadCount()
                
def is_hidden(tree):
    return getattr(tree, "hidden", "false") != "false" or getattr(tree, "hideDisplay", "false") != "false"

class TextInputBox(QWidget):
    def __init__(self, parent=None, flags=0, **kwargs):
        QWidget.__init__(self, parent)
        layout = QVBoxLayout()
        self.textEdit = QTextEdit(self)
        layout.addWidget(self.textEdit)
        self.setLayout(layout)
        self.button = QPushButton("Browse", self)
        layout.addWidget(self.button)
        self.setLayout(layout)
        
        self.connect(self.button, SIGNAL("clicked()"), self._openFile)
        
    def _openFile(self):
        file = QFileDialog.getOpenFileName(self, "Open File")
        data = open(file, "rb").read()
        self.textEdit.setText(data)
        
    def text(self):
        return str(self.textEdit.toPlainText())
    
class IdListInputBox(TextInputBox):
    def __init__(self, *args):
        TextInputBox.__init__(self, *args)
        self.validator = QRegExpValidator(QRegExp(""), self)
         
    def list(self):
        return self.text().split()
#        return [id for id in self.text().split() if self.validator.validate(QString(id)) == QValidator.Acceptable]
    
class DropDownMenu(QComboBox):
    def __init__(self, tree, *args, **kwargs):
        QComboBox.__init__(self, *args, **kwargs)
        self.tree = tree
        self.model = model = [] #ListModel()
        self.setOptions(tree.options())
#        print "Menu options", tree.options(), tree._tree
#        for options in tree.options():
#            if not is_hidden(option):
#                model.append((option.displayName, option.field, option.value, option.push_actions()))
        self.connect(self, SIGNAL("activated(int)"), self.setSelected)
        
    def setOptions(self, options):
        self.clear()
        self.model = model = []
        for option in options:
            if not is_hidden(option):
                model.append((option.displayName, self.tree.field, option.value, option.push_actions()))
                self.addItem(option.displayName)
        # TODO: schedule call setSelected  
                
    def setSelected(self, index):
        self.index = index
        actions = self.model[self.index][-1]
#        print actions
        for action in actions:
            ref = action.ref
            control = self.configutation.get_by_name(ref)
            control.setModel(action.options())
        
    def get_configuration(self):
        return [("filter", "", self.tree, self.model[self.index][2])]
    
class BooleanFilter(QGroupBox):
    def __init__(self, tree, *args, **kwargs):
        self.tree = tree
        self.bg = bg = QButtonGroup(self)
        self.setOptions(tree.options())
#        self.model = model = []
#        for option in tree.options():
#            model.append((option.displayName, option.field, option.value, option.push_actions()))
#            rb = QRadioButton(option.displayName, self)
#            self.layout().addWidget(rb)
#            self.bg.addButton(rb)
        
        self.connect(self.bg, SIGNAL("buttonClicked(int)"), self.setSelected)
        
    def setOptions(self, options):
        
        self.model = model = []
        for option in options:
            model.append((option.displayName, option.field, option.value, option.push_actions()))
            rb = QRadioButton(option.displayName, self)
            self.layout().addWidget(rb)
            self.bg.addButton(rb)
            
    def setSelected(self, id):
        self.index = id
        
    def get_configuration(self):
        return [("filter", "", self.tree, self.model[self.index][2])]
            
class BooleanListFilter(QWidget):
    def __init__(self, tree, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.tree = tree
        
        layout = QHBoxLayout()
        self.setLayout(layout)
        combo = QComboBox(self)
        box = QVGroupBox(self)
        box.setFlat(True)
        bg = QButtonGroup(self)
        layout.addWidget(combo)
        layout.addWidget(box)
        self.setLayout(layout)
        self.model = model = []
        self.bool_widgets = []
        for option in tree.options():
            model.append((option.displayName, option.field, option.value, option.push_actions()))
            combo.addItem(option.displayName)
            self.bool_widgets.append(BooleanFilter(option))
        
        self.connect(combo, SIGNAL("activated(int)"), self.setSelected)
        
    def setOptions(self, options):
        self.model = model = []
        for option in options:
            pass
    def setSelected(self, id):
        self.index = id
        
    def get_configuration(self):
        w = self.bool_widgets[self.index]
        _, _, value, _ = w.get_configuration()[0]
        return [("filter", "", self.tree, value)]
        
class ConfigurationControl(object):
    def __init__(self, tree, registry):
        self.registry = registry
        self.__i = 0
        self.__tree = {}
        self.__sub_controls = []
        
    def new_name(self, tree):
        name = "_%i" % self.__i
        self.__i += 1
        self.__tree[name] = tree
        setattr(self, name, True if getattr(tree, "default", "false") == "true" else False)
        return name
    
    def new_sub_control(self, tree, control):
#        name = self.new_name(tree)
#        setattr(self, name, control)
        self.__sub_controls.append(control)
#        return getattr(self, name)
    
    def is_hidden(self, tree):
        return getattr(tree, "hidden", "false") != "false" or getattr(tree, "hideDisplay", "false") != "false"
    
    def is_pointer(self, tree):
        return hasattr(tree, "pointerDataset")
    
    def display_name(self, tree):
        return getattr(tree, "displayName", "")
    
    def get_configuration(self):
        conf = []
        for key, tree in sorted(self.__tree.items()):
            conftype = "attribute" if isinstance(tree, BioMartAttribute) else "filter"
            dataset = getattr(self, "dataset", getattr(tree, "dataset", ""))
            name = tree
            val = getattr(self, key)
#            print conftype, dataset, tree ,val
            if val:
                conf.append((conftype, dataset, tree, val))
            
#        print conf
            
        return reduce(list.__add__, [sub.get_configuration() for sub in self.__sub_controls], []) + conf
        
    def get_pointed(self, tree):
        if hasattr(tree, "pointerFilter"):
            dataset, name, getter = tree.pointerDataset, tree.pointerFilter, "filter"
        elif hasattr(tree, "pointerAttribute"):
            dataset, name, getter = tree.pointerDataset, tree.pointerAttribute, "attribute"
            
        conf = self.registry.connection.configuration(dataset=dataset)
        return dataset, getattr(conf, getter)(name)

class OWBioMartConfigurationControl(OWBaseWidget, ConfigurationControl):
    def __init__(self, configurationTree, registry, parent):
        OWBaseWidget.__init__(self, parent)
        ConfigurationControl.__init__(self, configurationTree, registry)
        self.setWindowFlags(Qt.Widget)
        
class OWBioMartFilterContainerText(OWBioMartConfigurationControl):
    def __init__(self, filter, register, parent):
        OWBioMartConfigurationControl.__init__(self, filter, register, parent)
        self.tree = filter
        layout = QVBoxLayout()
        self.optionIndex = 0
        self.cb = OWGUI.comboBox(self, self, "optionIndex",  orientation=layout, items=[], sendSelectedValue=False)
        self.text_filter = text_filter = IdListInputBox(self)
        text_filter.layout().setContentsMargins(0,0,0,0)
        layout.addWidget(text_filter)
        self.setLayout(layout)
        self.options = []
#        print filter, filter.__dict__
        for option in filter.options():
#            print option
            dname = getattr(option, "displayName", "")
            regexp = getattr(option, "regexp", "")
            field = getattr(option, "field", "")
            self.options.append((dname, regexp, field, option))
            
        self.cb.addItems([opt[0] for opt in self.options])
#        self.connect(self.cb, SIGNAL("selectionChanged()"), lambda :setattr(text_filter, "validator", None))

    def get_configuration(self):
        option = self.options[self.optionIndex]
        return [("filter", getattr(self, "dataset", getattr(option, "dataset", "")), option[2], self.text_filter.list())]
    
class OWBioMartFilterContainerBooleanList(OWBioMartConfigurationControl):
    def __init__(self, filter, register, parent):
        OWBioMartConfigurationControl.__init__(self, filter, register, parent)
        layout = QGridLayout()
        self.optionIndex = 0
        self.cb = OWGUI.comboBox(self, self, "optionIndex", addToLayout=False)
        layout.addWidget(self.cb, 0, 0, 2, 1)
        w = OWGUI.widgetBox(self, "", addToLayout=False)
        layout.addWidget(w, 0, 1, 2, 1)
#        print filter
        self.options = []
        for option in getattr(filter._tree, "Option", []):
            name = option.attributes.get("displayName","")
            field = option.attributes.get("field","")
            values = []
            for op in getattr(option, "Option", []):
                n = op.attributes.get("displayName","")
                v = op.attributes.get("value","")
                values.append((n, v))
            self.options.append((name, field, values))
        self.valueIndex = 0
        self.bg = OWGUI.radioButtonsInBox(w, self, "valueIndex", [], box=w)
        self.cb.addItems([opt[0] for opt in self.options])
        self.setLayout(layout) 
        self.setOption()
        
    def setOption(self):
        name, field, values = self.options[self.optionIndex]
        for n, value in values:
            OWGUI.appendRadioButton(self.bg, self, "valueIndex", n)# tooltip, insertInto, callback, addToLayout)
        
class OWBioMartFilterControl(OWBioMartConfigurationControl):
    def __init__(self, filter, register, parent):
        OWBioMartConfigurationControl.__init__(self, filter, register, parent)
        self.setLayout(QVBoxLayout())
        filter_type = self.filter_type(filter)
        
        if filter_type == ("text", "", ""):
            w = OWGUI.lineEdit(self, self, self.new_name(filter), "")
                 
        elif filter_type == ("text", "", "1"):
            w = IdListInputBox(self)
            self.layout().addWidget(w)
            
        elif filter_type == ("list", "radio", ""):
            options = [getattr(opt, "displayName", "") for opt in self.filter_options(filter)]
            w = OWGUI.radioButtonsInBox(self, self, self.new_name(filter), options)
            

        elif filter_type == ("list", "menu", ""):
            w = DropDownMenu(filter, self)
            self.layout().addWidget(w)
#            options = [getattr(opt, "displayName", "") for opt in self.filter_options(filter)]
#            w = OWGUI.comboBox(self, self, self.new_name(filter), items=options, sendSelectedValue=True, valueType=str)
            
        elif filter_type == ("list", "menu", "1"):
            name = self.new_name(filter)
            self.options = [getattr(opt, "displayName", "") for opt in self.filter_options(filter)]
            setattr(self, name, [])
            w = OWGUI.listBox(self, self, name, "options", selectionMode=QListWidget.MultiSelection)
        elif filter_type == ("container", "", ""):
#            print getattr(filter, "displayName", ""), filter._tree
#            print filter, filter._tree
            options = self.filter_options(filter)
            filter_type = self.filter_type(options[0]) if options else None
            type = getattr(filter, "type", "")
            
            if filter_type == ("text", "", "1") and type == "id_list":
                w = OWBioMartFilterContainerText(filter, register, self)
                self.new_sub_control(filter, w)
                self.layout().addWidget(w)
            elif filter_type == ("text", "", "1") and type == "list":
                options = [getattr(opt, "displayName", "") for opt in self.filter_options(filter)]
                w = OWGUI.comboBox(self, self, self.new_name(filter), items=options, sendSelectedValue=True, valueType=str)
            elif filter_type == ("text", "", ""):
                w = OWBioMartFilterContainerText(filter, register, self)
                self.new_sub_control(filter, w)
                self.layout().addWidget(w)
            elif filter_type == ("list", "radio", ""):
                w = OWBioMartFilterContainerBooleanList(filter, register, self)
                self.layout().addWidget(w)
                self.new_sub_control(filter, w)
            else:
                print "Unknown filter type", filter, getattr(filter, "internalName", "")
        else:
            print "Unknown filter type", filter, getattr(filter, "internalName", "")
            
    def filter_type(self, tree):
        return getattr(tree, "displayType", "text"), getattr(tree, "style", ""), getattr(tree, "multipleValues", "")
    
    def filter_options(self, tree):
        return [option for option in tree.options()]
    
    def new_name(self, filter, value=""):
        name = ConfigurationControl.new_name(self, filter)
        
        setattr(self, name, getattr(filter, "defaultValue", value)) #getattr(self, name)))
        return name
    
    def get_configuration(self):
        conf = OWBioMartConfigurationControl.get_configuration(self)
        return conf

class OWBioMartFilterTextField(OWBioMartFilterControl):
    def __init__(self, filter, registry, parent):
        OWBioMartFilterControl.__init__(self, filter, registry, parent)
        
class OWBioMartFilterCollection(OWBioMartConfigurationControl):
    def __init__(self, collection, registry, parent):
        OWBioMartConfigurationControl.__init__(self, collection, registry, parent)
        self.checked = False
        self.buildCollection(collection)
        
    def buildCollection(self, collection):
        collection_widget = self #OWGUI.widgetBox(widget, getattr(collection, "displayName", ""), flat=True, orientation=QGridLayout())
        
        self.checked = True if getattr(collection, "default", "false") == "true" else False
        cb = OWGUI.checkBox(collection_widget, self, "checked", self.display_name(collection))
        
        if isinstance(collection, BioMartFilterCollection):
            descriptions = [filter for filter in collection.filters() if not self.is_hidden(filter)]
        else:
            descriptions = [collection]
            
        grid = QGridLayout()
        if len(descriptions) > 1:
            self.setLayout(QVBoxLayout())
            self.layout().addWidget(cb)
            
            option = QStyleOptionButton()
            option.initFrom(cb)
        
            indent = qApp.style().subElementRect(QStyle.SE_CheckBoxIndicator, option, cb).width() + 3
            
            OWGUI.indentedBox(self, indent, grid, addSpace=False)
            for filter in descriptions:
                if self.is_pointer(filter):
                    (dataset, filter), pointer = self.get_pointed(filter), filter
                    if not filter:
                        continue
                w = OWBioMartFilterControl(filter, self.registry, self)
                self.layout().addWidget(w)
                self.new_sub_control(filter, w)
            
        else:
            self.setLayout(grid)
            grid.addWidget(cb, 0, 0)
            filter = descriptions[0]
            if self.is_pointer(filter):
                (dataset, filter), pointer = self.get_pointed(filter), filter
                if not filter:
                    return
            w = OWBioMartFilterControl(filter, self.registry, self)
            grid.addWidget(w, 0, 1)
            self.new_sub_control(filter, w)
            
        
    def get_configuration(self):
#        print "collection", self.checked
        if self.checked:
            conf = OWBioMartConfigurationControl.get_configuration(self)
#            if getattr(self, "dataset", None):
            return [(conftype, dataset if dataset else getattr(self, "dataset", dataset), tree, val) \
                    for  conftype, dataset, tree, val in conf]
        return []
    
class OWBioMartConfigurationPage(OWBioMartConfigurationControl):
    def __init__(self, pageTree, registry, parent=None,):
        OWBioMartConfigurationControl.__init__(self, pageTree, registry, parent)
        
        self.setLayout(QVBoxLayout(self))
        self.out_format = getattr(pageTree, "outFormats", "tsv")
         
        self.buildPage(pageTree)
    
    def buildPage(self, page):
            
        def buildAttributeControl(widget, attr):
            cb = OWGUI.checkBox(widget, self, self.new_name(attr), getattr(attr, "displayName", ""), addToLayout=False)
            return cb
        
        def buildFilterControl(widget, filter):
            w = OWBioMartFilterCollection(filter, self.registry, widget)
            self.new_sub_control(filter, w)
            return w
        
        def buildFilterCollection(widget, collection):
            w = OWBioMartFilterCollection(collection, self.registry, widget)
            self.new_sub_control(collection, w)
            widget.layout().addWidget(w, 0, 0)

        def buildAttributeCollection(widget, collection):
            descriptions = [attr for attr in collection.attributes() if not self.is_hidden(attr)]
            for i, desc in enumerate(descriptions):
                dataset, pointer = "", None
                if self.is_pointer(desc):
                    (dataset, desc), pointer = self.get_pointed(desc), desc
                    if not desc:
#                        print pointer.__dict__
                        continue
                    desc.dataset = dataset
                if isinstance(desc, BioMartAttribute):
                    control = buildAttributeControl(widget, desc)
                    widget.layout().addWidget(control, i % (max(len(descriptions) / 2 + 1, 1)), i / (len(descriptions) / 2 + 1))
                else:
                    control = buildFilterControl(widget, desc)
                    widget.layout().addWidget(control, i, 0, 1, 2)
                    control.dataset = dataset
                
                        
        def buildCollection(widget, collection):
            collection_widget = OWGUI.widgetBox(widget, getattr(collection, "displayName", ""), flat=True, orientation=QGridLayout())
            if isinstance(collection, BioMartAttributeCollection):
                buildAttributeCollection(collection_widget, collection)
            else:
                buildFilterCollection(collection_widget, collection)
                
        def buildGroup(widget, group):
            group_widget = OWGUI.collapsableWidgetBox(widget, getattr(group, "displayName", ""))
            for coll in group.attribute_collections() if isinstance(group, BioMartAttributeGroup) else \
                            group.filter_collections():
                if not self.is_hidden(coll):
                    buildCollection(group_widget, coll)
            
        for group in page.attribute_groups() if isinstance(page, BioMartAttributePage) else page.filter_groups():
            if not self.is_hidden(group):
                buildGroup(self, group)
                
        OWGUI.rubber(self) 
    
#    def get_attribute_configuration(self):
#        return [(getattr(self, name), desc) for name, desc in self.__tree.items() if isinstance(desc, BioMartAttribute)]
#    
#    def get_filter_configuration(self):
#        return [(getattr(self, name), desc) for name, desc in self.__tree.items() if isinstance(desc, BioMartFilter)]

       
class OWBioMart(OWWidget):
    def __init__(self, parent=None, signalManager=None, title="Bio Mart"):
        OWWidget.__init__(self, parent, signalManager, title)
        
        self.outputs = [("Example Table", ExampleTable)]
        
        self.selectedDatabase = None
        self.selectedDataset = None
        
        self.martsCombo = OWGUI.comboBox(self.controlArea, self, "selectedDatabase", "Database", callback=self.setSelectedMart)
        self.martsCombo.setMaximumWidth(250)
        self.datasetsCombo = OWGUI.comboBox(self.controlArea, self, "selectedDataset", "Dataset", callback=self.setSelectedDataset)
        self.datasetsCombo.setMaximumWidth(250)
        OWGUI.rubber(self.controlArea)
        
        OWGUI.button(self.controlArea, self, "Commit", callback = self.commit, tooltip="Commit the selected dataset to output")
        
        self.mainTab = OWGUI.tabWidget(self.mainArea)
        
        self.attributesConfigurationBox = OWGUI.createTabPage(self.mainTab, "Attributes", canScroll=True)
        #OWGUI.widgetBox(self.mainArea, "Attributes")
        self.filtersConfigurationBox = OWGUI.createTabPage(self.mainTab, "Filters", canScroll=True)
        #OWGUI.widgetBox(self.mainArea, "Filters")
#        OWGUI.rubber(self.mainArea)
        
        self.initTask = MartserviceTask(self._get_registry)
#        self.initTask.start()
        
#        self.connect(self.initTask, SIGNAL("executionFailed(PyQt_PyObject)"), self._error, Qt.QueuedConnection)
#        self.connect(self.initTask, SIGNAL("executionFinished(PyQt_PyObject)"), self.setBioMartRegistry, Qt.QueuedConnection)
        self.connect(self.initTask, SIGNAL("finished(QString)"), self.setBioMartRegistry, Qt.QueuedConnection)
        
#        self.datasetListTask = BackgroundTask()
#        self.datasetListTask.start()
        
#        self.connect(self.datasetListTask, SIGNAL("executionFailed"), self._error, Qt.QueuedConnection)
#        self.connect(self.datasetListTask, SIGNAL("executionFinished"), self.setBioMartDatasets, Qt.QueuedConnection)
        
#        self.configurationTask = BackgroundTask()
#        self.configurationTask.start()
        
#        self.connect(self.configurationTask, SIGNAL("executionFailed"), self._error, Qt.QueuedConnection)
#        self.connect(self.configurationTask, SIGNAL("executionFinished"), self.setBioMartConfiguration, Qt.QueuedConnection)
        
#        self.initTask.runTask(self._get_registry)
        self.initTask.run()
        
        self.resize(600, 400)
        
    @staticmethod
    def _get_registry(url=None, precache=True):
        con = obiBioMart.BioMartConnection(url)
        reg = con.registry()
        if precache:
            marts = reg.marts()
        return reg
    
    def _error(self, ex):
        import traceback
        print ex
    
    def setBioMartRegistry(self, state):
        if state == "Ok":
            self.registry = self.initTask.result
            self.marts = [mart for mart in self.registry.marts() if getattr(mart, "visible", "0") != "0"] 
            for mart in self.marts:
                self.martsCombo.addItem(mart.displayName)
        else:
            print state
            
    def setSelectedMart(self):
        self.mart = self.marts[self.selectedDatabase]
        
        self.datasetListTask = MartserviceTask(self.mart.datasets)
#        self.connect(self.datasetListTask, SIGNAL("executionFailed(PyQt_PyObject)"), self._error, Qt.QueuedConnection)
#        self.connect(self.datasetListTask, SIGNAL("executionFinished(PyQt_PyObject)"), self.setBioMartDatasets, Qt.QueuedConnection)
        self.connect(self.datasetListTask, SIGNAL("finished(QString)"), self.setBioMartDatasets, Qt.QueuedConnection)
        
        self.datasetListTask.run()
#        qApp.processEvents()
        
    def setBioMartDatasets(self, state):
        if state == "Ok":
            datasets = self.datasetListTask.result
            self.datasets = [data for data in datasets if getattr(data,"visible", "0") != "0"]
            self.datasetsCombo.clear()
            self.datasetsCombo.addItems([data.displayName for data in self.datasets])
            
    def setSelectedDataset(self):
        self.dataset = self.datasets[self.selectedDataset]
        
        self.configurationTask = MartserviceTask(self.dataset.configuration)
#        self.connect(self.configurationTask, SIGNAL("executionFailed(PyQt_PyObject)"), self._error, Qt.QueuedConnection)
#        self.connect(self.configurationTask, SIGNAL("executionFinished(PyQt_PyObject)"), self.setBioMartConfiguration, Qt.QueuedConnection)
        self.connect(self.configurationTask, SIGNAL("finished(QString)"), self.setBioMartConfiguration, Qt.QueuedConnection)
        
        self.configurationTask.run()
#        self.configurationTask.runTask(self.dataset.configuration)
        
    def setBioMartConfiguration(self, state):
        print QThread.currentThread()
        if state == "Ok":
            self.clearConfiguration()
        
            self.configuration = configuration = self.configurationTask.result
            
            def hidden(tree):
                return getattr(tree, "hidden", "false") != "false" or getattr(tree, "hideDisplay", "false") != "false"
            
            self.attributePagesTabWidget = tabs =  OWGUI.tabWidget(self.attributesConfigurationBox)
            self.registry.connection.configuration = cached_function_instance(self.registry.connection.configuration)
            for page in configuration.attribute_pages():
                if not hidden(page):
                    page_widget = OWBioMartConfigurationPage(page, self.registry, self)
                    OWGUI.createTabPage(tabs, getattr(page, "displayName", ""), widgetToAdd=page_widget )
                    
            self.filterPagesTabWidget = tabs = OWGUI.tabWidget(self.filtersConfigurationBox)
            for page in configuration.filter_pages():
                if not hidden(page):
                    page_widget = OWBioMartConfigurationPage(page, self.registry, self)
                    OWGUI.createTabPage(tabs, getattr(page, "displayName", ""), widgetToAdd=page_widget )
    
    def clearConfiguration(self):
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
        
    def count(self):
        filters = []
        count = self.dataset.get_count(filters = filters)
        
    def commit(self):
        pageconf = self.attributePagesTabWidget.currentWidget()
        format = pageconf.out_format
        configuration = pageconf.get_configuration()
        bydatasets = defaultdict(lambda : ([], []))
#        print configuration
        version = getattr(self.configuration, "softwareVersion", "0.4")
        for conftype, dataset, name, val in configuration:
            if version > "0.4":
                dataset = self.dataset
            if conftype == "attribute":
                bydatasets[dataset][0].append(name)
            elif conftype == "filter":
                bydatasets[dataset][1].append((name, val))
                
        pageconf = self.filterPagesTabWidget.currentWidget()
        configuration =  pageconf.get_configuration()
        print configuration
#        bydatasets = defaultdict(lambda : ([], []))
        for conftype, dataset, tree, val in configuration:
            if version > "0.4":
                dataset = self.dataset
            if conftype == "attribute":
                bydatasets[dataset][0].append(tree)
            elif conftype == "filter":
                bydatasets[dataset][1].append((tree, val))
            
        query = self.registry.query(format="TSV" if "tsv" in format.lower() else format.upper())
        
        for dataset, (attributes, filters) in bydatasets.items():
            query.set_dataset(dataset if dataset else self.dataset)
            for attr in attributes:
                query.add_attribute(attr)
            for filter, value in filters:
                query.add_filter(filter, value)
        
        print query.xml_query()
#        data = query.run()
        self.dataTask = MartserviceTask(query.get_example_table)
        self.connect(self.dataTask, SIGNAL("finished(QString)"), self.dataReady, Qt.QueuedConnection)
        self.dataTask.run()
        
#        data = query.get_example_table()
#        print data
        
#        self.send("Example Table", data)

    def dataReady(self, state):
        if state == "Ok":
            data = self.dataTask.result
            print "Data", data.domain
            self.send("Example Table", data)
        else:
            print state
        
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    ow = OWBioMart()
    ow.show()
    app.exec_()
    import gc
    print gc.garbage

        
        
        
        