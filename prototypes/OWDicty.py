"""
<name>Dicty database</name>
<description>Widget interface to dicty experaments</description>
"""
from __future__ import with_statement

from OWWidget import *
import obiDicty
import OWGUI

from functools import partial
from collections import defaultdict
from threading import Lock, Thread
from traceback import print_exception
import sys

class ShadedItemDelegate(QItemDelegate):
    def paint(self, painter, option, index):
        self.drawBackground(painter, option, index)
        value, ok = index.data(Qt.DisplayRole).toString()

class MyQThread(QThread):
    def __init__(self, parent, func=None):
        QThread.__init__(self, parent)
        self.func = func
        
    def run(self):
        try:
            self.returnValue = self.func()
        except Exception, ex:
            pass
##            print_exception(*sys.exc_info())
        
class OWDicty(OWWidget):
    settingsList = ["serverToken", "platform", "tables", "selectedTables"]
    def __init__(self, parent=None, signalManager=None, name="Dicty database"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.outputs = [("Example tables", ExampleTable, Multiple)]
        self.serverToken = ""
##        self.type = None
##        self.typesList = ["norms", "chips"] #[n for n, long in obiDicty.DatabaseConnection.obidPairs]

        self.platform = None
        self.platformList = []
        
        self.sampleList = []
        self.sampleSelected = []

        self.treatmentList = []
        self.treatmentSelected = []

        self.growthCondList = []
        self.growthCondSelected = []
        
        self.joinList = [n for n, long in obiDicty.DatabaseConnection.aoidPairs]
        self.joinSelected = []
        
        self.separateList = [n for n, long in obiDicty.DatabaseConnection.aoidPairs]
        self.separateSelected = []
        
        self.loadSettings()
        
        OWGUI.lineEdit(self.controlArea, self, "serverToken", box="Server Token", callback=self.Connect)
##        self.typeCombo = OWGUI.comboBox(self.controlArea, self, "type", box="Type", items=self.typesList, sendSelectedValue=True, callback=self.UpdateControls)
        self.platformCombo = OWGUI.comboBox(self.controlArea, self, "platform", box="Platform", items=self.platformList, sendSelectedValue=True, callback=partial(self.UpdateControls, ["sample", "treatment", "growthCond"]))
        self.sampleListBox = OWGUI.listBox(self.controlArea, self, "sampleSelected", "sampleList", box="Samples", selectionMode=QListWidget.ExtendedSelection, callback=partial(self.UpdateControls, ["treatment", "growthCond"]))
        self.treatmentListBox = OWGUI.listBox(self.controlArea, self, "treatmentSelected", "treatmentList", box="Treatment", selectionMode=QListWidget.ExtendedSelection, callback=partial(self.UpdateControls, ["growthCond"]))
        self.growthCondListBox = OWGUI.listBox(self.controlArea, self, "growthCondSelected", "growthCondList", box="Growth Condition", selectionMode=QListWidget.ExtendedSelection, callback=partial(self.UpdateControls, []))
        self.joinListBox = OWGUI.listBox(self.mainArea, self, "joinSelected", "joinList", box="Join", selectionMode=QListWidget.ExtendedSelection)
        self.separateListBox = OWGUI.listBox(self.mainArea, self, "separateSelected", "separateList", box="Separate", selectionMode=QListWidget.ExtendedSelection)
        
##        OWGUI.button(self.controlArea, self, "&Preview", callback=self.ShowPreview)
        OWGUI.button(self.controlArea, self, "&Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)

        self.dbc = None
        self.locks = defaultdict(Lock)
        self.InitControls()

    def __updateSelectionList(self, oldList, oldSelection, newList):
        oldList = [oldList[i] for i in oldSelection]
        return [ i for i, new in enumerate(newList) if new in oldList]

    def GetOptions(self):
        print self.platform, [self.sampleList[i] for i in self.sampleSelected], [self.joinList[i] for i in self.joinSelected], [self.separateList[i] for i in self.separateSelected]
        return self.platform, [self.sampleList[i] for i in self.sampleSelected], [self.joinList[i] for i in self.joinSelected], [self.separateList[i] for i in self.separateSelected]

    def GetQuery(self):
        query = {}
        if self.platform:
            query["platform"] = self.platform
        if self.sampleSelected:
            query["sample"] = [self.sampleList[i] for i in self.sampleSelected]
        if self.treatmentSelected:
            query["treatment"] = [self.treatmentList[i] for i in self.treatmentSelected]
        if self.growthCondSelected:
            query["growthCond"] = [self.growthCondList[i] for i in self.growthCondSelected]
        return {} #query
    
    def Connect(self):
        address = "http://www.ailab.si/dictyexpress/api/index.php?"
        if self.serverToken:
            address += "token="+self.serverToken+"&"
        try:
            self.dbc = obiDicty.DatabaseConnection(address)
        except Exception, ex:
            from traceback import print_exception
            print_exception(*sys.exc_info())
            self.error(0, "Error connecting to server" + str(ex))
            return
        self.error(0)

    def InitControls(self):
        if not self.dbc:
            self.Connect()
        self.platformThread = MyQThread(self, partial(self.dbc.annotationOptions, self.dbc.aoidt("platform")))
        self.platformCombo.setDisabled(True)
        def _set():
            self.platformCombo.clear()
            self.platformCombo.addItems(self.platformThread.returnValue["platform"])
            self.platformCombo.setDisabled(False)
        self.connect(self.platformThread, SIGNAL("finished()"), partial(_set))
        self.platformThread.start()
            
    def UpdateControls(self, controls=()):
        if not self.dbc:
            self.Connect()
        query = self.GetQuery()
        for control in controls:
            thread = MyQThread(self, partial(self.dbc.annotationOptions, self.dbc.aoidt(control), **query))
            setattr(self, control+"Thread", thread)
            listBox = getattr(self, control+"ListBox")
            listBox.setDisabled(True)
            def _set(control):
                setattr(self, control+"List", getattr(self, control+"Thread").returnValue[control])
                getattr(self, control+"ListBox").setDisabled(False)
            self.connect(thread, SIGNAL("finished()"), partial(_set, control))
            thread.start()

    def ShowPreview(self):
        pass

    def Commit(self):
        platform, sample, join, separate = self.GetOptions()
        tables = self.dbc.getData(platform=platform, sample=sample, join=join, separate=separate)
        self.send("Example tables", None)
        for i, table in enumerate(tables):
            self.send("Example tables", table, i)

if __name__ == "__main__":
    app  = QApplication(sys.argv)
    w = OWDicty()
    w.show()
    app.exec_()
    w.saveSettings()
            
        