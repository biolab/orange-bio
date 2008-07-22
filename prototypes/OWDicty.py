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
    settingsList = ["serverToken", "platform", "platformList", "experiments", "selectedExperiments"]
    def __init__(self, parent=None, signalManager=None, name="Dicty database"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.outputs = [("Example tables", ExampleTable, Multiple)]
        self.serverToken = ""

        self.platform = None
        self.platformList = []
        
        self.experiments = []
        self.selectedExperiments = []
        
        self.loadSettings()
        
        OWGUI.lineEdit(self.controlArea, self, "serverToken", box="Server Token", callback=self.Connect)
##        OWGUI.lineEdit(self.mainArea, self, "search", 
        self.experimentsWidget = QTreeWidget()
        self.experimentsWidget.setHeaderLabels(["Strain", "Treatment", "Growth condition"])
        self.experimentsWidget.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.experimentsWidget.setRootIsDecorated(False)
##        self.experimentsWidget.setAlternatingRowColors(True)
        self.mainArea.layout().addWidget(self.experimentsWidget)
##        OWGUI.button(self.controlArea, self, "&Preview", callback=self.ShowPreview)
        OWGUI.button(self.controlArea, self, "&Update", callback=self.UpdateExperiments)
        OWGUI.button(self.controlArea, self, "&Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)

        self.dbc = None        

        self.FillExperimentsWidget()

        self.resize(600, 400)

    def __updateSelectionList(self, oldList, oldSelection, newList):
        oldList = [oldList[i] for i in oldSelection]
        return [ i for i, new in enumerate(newList) if new in oldList]
    
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

    def UpdateExperiments(self):
        if not self.dbc:
            self.Connect()
        self.experiments = []
        self.experimentsWidget.clear()
        self.progressBarInit()
        strains = self.dbc.annotationOptions(self.dbc.aoidt("sample"))["sample"]
        for i, strain in enumerate(strains):
            opt = self.dbc.annotationOptions(sample=strain)
            treatments = opt["treatment"]
            growthConds = opt["growthCond"]
            for treatment in treatments:
                for cond in growthConds:
                    self.experiments.append([strain, treatment, cond])
                    QTreeWidgetItem(self.experimentsWidget, self.experiments[-1])
            self.progressBarSet((100.0 * i) / len(strains))
        self.progressBarFinished()

    def FillExperimentsWidget(self):
        self.experimentsWidget.clear()
        for strings in self.experiments:
            QTreeWidgetItem(self.experimentsWidget, strings)

    def ShowPreview(self):
        pass

    def Commit(self):
        if not self.dbc:
            self.Connect()
        allTables = []
##        print self.experimentsWidget.selectedIndexes()
        for item in self.experimentsWidget.selectedItems():
##            item = self.experimentsWidget.itemFromIndex(index)
            print str(item.text(0)), str(item.text(1)), str(item.text(2))
            tables = self.dbc.getData(sample=str(item.text(0)), treatment=str(item.text(1)), growthCond=str(item.text(2)))
##            print len(tables)
            allTables.extend(tables)
        self.send("Example tables", None)
        for i, table in enumerate(allTables):
            self.send("Example tables", table, i)

if __name__ == "__main__":
    app  = QApplication(sys.argv)
    w = OWDicty()
    w.show()
    app.exec_()
    w.saveSettings()
            
        