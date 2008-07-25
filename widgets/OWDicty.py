"""
<name>Dicty database</name>
<description>Widget interface to dicty experaments</description>
"""

from OWWidget import *
import obiDicty
import OWGUI

import sys

class MyTreeWidgetItem(QTreeWidgetItem):
    def __contains__(self, text):
        return any(text.upper() in str(self.text(i)).upper() for i in range(self.columnCount()))    
    
class OWDicty(OWWidget):
    settingsList = ["serverToken", "platform", "platformList", "experiments", "selectedExperiments", "joinSelected", "separateSelected"]
    def __init__(self, parent=None, signalManager=None, name="Dicty database"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.outputs = [("Example tables", ExampleTable, Multiple)]
        self.serverToken = ""

        self.platform = None
        self.platformList = []

        self.joinList = [n for n, long in obiDicty.DatabaseConnection.aoidPairs]
        self.joinSelected = [0]
        
        self.separateList = [n for n, long in obiDicty.DatabaseConnection.aoidPairs]
        self.separateSelected = [4, 5]
        
        self.experiments = []
        self.selectedExperiments = []

        self.searchString = ""
        
        OWGUI.lineEdit(self.controlArea, self, "serverToken", box="Server Token", callback=self.Connect)
        
        self.joinListBox = OWGUI.listBox(self.controlArea, self, "joinSelected", "joinList", box="Join By", selectionMode=QListWidget.ExtendedSelection, callback=self.UpdateJoinSelected)
        self.separateListBox = OWGUI.listBox(self.controlArea, self, "separateSelected", "separateList", box="Separate By", selectionMode=QListWidget.ExtendedSelection, callback=self.UpdateJoinSelected)
        OWGUI.button(self.controlArea, self, "&Update list", callback=self.UpdateExperiments)
##        box = OWGUI.widgetBox(self.controlArea, "Data")
##        OWGUI.checkBox(box, self, "useCache", "Use cached data", 
        OWGUI.button(self.controlArea, self, "&Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)

        OWGUI.lineEdit(self.mainArea, self, "searchString", "Search", callbackOnType=True, callback=self.SearchUpdate)
        self.experimentsWidget = QTreeWidget()
        self.experimentsWidget.setHeaderLabels(["Strain", "Treatment", "Growth condition", "Platform"]) #, "Num. of replications", "Num. of tech replications"])
        self.experimentsWidget.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.experimentsWidget.setRootIsDecorated(False)
        self.experimentsWidget.setSortingEnabled(True)
##        self.experimentsWidget.setAlternatingRowColors(True)
        self.mainArea.layout().addWidget(self.experimentsWidget)
##        OWGUI.button(self.controlArea, self, "&Preview", callback=self.ShowPreview)

        self.loadSettings()
        self.dbc = None        

        self.UpdateJoinSelected()        
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

    def UpdateJoinSelected(self):
        for i, item in [(i, self.joinListBox.item(i)) for i in range(self.joinListBox.count())]:
            self.separateListBox.item(i).setHidden(item.isSelected())
        for i, item in [(i, self.separateListBox.item(i)) for i in range(self.separateListBox.count())]:
            self.joinListBox.item(i).setHidden(item.isSelected())

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
            platforms = opt["platform"]
            print opt
            for treatment in treatments:
                for cond in growthConds:
                    for platform in platforms:
                        self.experiments.append([strain, treatment, cond, platform]) #, str(len(opt["replicate"])), str(len(opt["techReplicate"]))])
                        self.items.append(MyTreeWidgetItem(self.experimentsWidget, self.experiments[-1]))
            self.progressBarSet((100.0 * i) / len(strains))
        self.progressBarFinished()

    def FillExperimentsWidget(self):
        if not self.experiments:
            self.UpdateExperiments()
            return
        self.experimentsWidget.clear()
        self.items = []
        for strings in self.experiments:
            self.items.append(MyTreeWidgetItem(self.experimentsWidget, strings))
        
    def SearchUpdate(self, string=""):
        for item in self.items:
            item.setHidden(not all(s in item for s in self.searchString.split()))
            
    
    def Commit(self):
        if not self.dbc:
            self.Connect()
        allTables = []
        join = [self.joinList[i] for i in self.joinSelected]
        separate = [self.separateList[i] for i in self.separateSelected]
        import time
        start = time.time()
        print "Start:", start
        for item in self.experimentsWidget.selectedItems():
            tables = self.dbc.getData(sample=str(item.text(0)), treatment=str(item.text(1)), growthCond=str(item.text(2)), join=join, separate=separate)
            for table in tables:
                table.name = ".".join([str(item.text(0)), str(item.text(1)), str(item.text(2)), str(item.text(3))])
            allTables.extend(tables)
        end = int(time.time()-start)
        print "End:","%ih:%im:%is" % (end/3600, (end/60)%60, end%60)
        
        self.send("Example tables", None)
        for i, table in enumerate(allTables):
            self.send("Example tables", table, i)

if __name__ == "__main__":
    app  = QApplication(sys.argv)
##    from pywin.debugger import set_trace
##    set_trace()
    w = OWDicty()
    w.show()
    app.exec_()
    w.saveSettings()
            
        