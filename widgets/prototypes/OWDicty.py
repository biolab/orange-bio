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

        self.searchString = ""
        
        self.loadSettings()
        
        OWGUI.lineEdit(self.controlArea, self, "serverToken", box="Server Token", callback=self.Connect)
        OWGUI.lineEdit(self.mainArea, self, "searchString", "Search", callbackOnType=True, callback=self.SearchUpdate)
        self.experimentsWidget = QTreeWidget()
        self.experimentsWidget.setHeaderLabels(["Strain", "Treatment", "Growth condition"])
        self.experimentsWidget.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.experimentsWidget.setRootIsDecorated(False)
        self.experimentsWidget.setSortingEnabled(True)
##        self.experimentsWidget.setAlternatingRowColors(True)
        self.mainArea.layout().addWidget(self.experimentsWidget)
##        OWGUI.button(self.controlArea, self, "&Preview", callback=self.ShowPreview)
        OWGUI.button(self.controlArea, self, "&Update list", callback=self.UpdateExperiments)
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
        self.items = []
        self.progressBarInit()
        strains = self.dbc.annotationOptions(self.dbc.aoidt("sample"))["sample"]
        for i, strain in enumerate(strains):
            opt = self.dbc.annotationOptions(sample=strain)
            treatments = opt["treatment"]
            growthConds = opt["growthCond"]
            for treatment in treatments:
                for cond in growthConds:
                    self.experiments.append([strain, treatment, cond])
                    self.items.append(QTreeWidgetItem(self.experimentsWidget, self.experiments[-1]))
            self.progressBarSet((100.0 * i) / len(strains))
        self.progressBarFinished()

    def FillExperimentsWidget(self):
        if not self.experiments:
            self.UpdateExperiments()
            return
        self.experimentsWidget.clear()
        self.items = []
        for strings in self.experiments:
            self.items.append(QTreeWidgetItem(self.experimentsWidget, strings))
        

    def ShowPreview(self):
        pass

    def SearchUpdate(self, string=""):
        print "s:",self.searchString
        for item in self.items:
            item.setHidden(not all(s in (item.text(0) + item.text(1) + item.text(2)) for s in self.searchString.split()))
            
        

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
            
        