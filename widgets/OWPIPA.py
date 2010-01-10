"""
<name>PIPA database</name>
<description>Interface to PIPA database.</description>
<icon>icons/unknown.png</icon>
<priority>260</priority>
"""

from OWWidget import *
import obiDicty
import OWGUI
import orngEnviron
import sys
from collections import defaultdict

class MyTreeWidgetItem(QTreeWidgetItem):
    def __contains__(self, text):
        return any(text.upper() in str(self.text(i)).upper() for i in range(self.columnCount()))    

#set buffer file
bufferpath = os.path.join(orngEnviron.directoryNames["bufferDir"], "pipa")
try:
    os.makedirs(bufferpath)
except:
    pass
bufferfile = os.path.join(bufferpath, "database.sq3")

class OWPIPA(OWWidget):
    settingsList = [ "platform", "selectedExperiments", "server", "buffertime", "excludeconstant" ]
    def __init__(self, parent=None, signalManager=None, name="Dicty database"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.outputs = [("Example table", ExampleTable)]

        self.platform = None

        self.selectedExperiments = []
        self.buffer = obiDicty.BufferSQLite(bufferfile)

        self.searchString = ""
        self.excludeconstant = False
        
        box = OWGUI.widgetBox(self.controlArea, "Cache")
        OWGUI.button(box, self, "Clear cache", callback=self.clear_buffer)

        OWGUI.checkBox(self.controlArea, self, "excludeconstant", "Exclude labels with constant values" )

        OWGUI.button(self.controlArea, self, "&Commit", callback=self.Commit)
        box  = OWGUI.widgetBox(self.controlArea, "Authentification")
        OWGUI.rubber(self.controlArea)

        OWGUI.lineEdit(self.mainArea, self, "searchString", "Search", callbackOnType=True, callback=self.SearchUpdate)
        self.experimentsWidget = QTreeWidget()
        self.experimentsWidget.setHeaderLabels(["Species", "Strain", "Genotype", "Replicate", "Timepoint", "Treatment", "Growth", "ID"])
        self.experimentsWidget.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.experimentsWidget.setRootIsDecorated(False)
        self.experimentsWidget.setSortingEnabled(True)
##        self.experimentsWidget.setAlternatingRowColors(True)

        self.mainArea.layout().addWidget(self.experimentsWidget)

        self.loadSettings()
        self.dbc = None        

        QTimer.singleShot(0, self.UpdateExperiments)        

        self.resize(800, 600)

    def __updateSelectionList(self, oldList, oldSelection, newList):
        oldList = [oldList[i] for i in oldSelection]
        return [ i for i, new in enumerate(newList) if new in oldList]
    
    def Connect(self):

        address = self.server
        try:
            #obiDicty.verbose = 1
            self.dbc = obiDicty.PIPA(buffer=self.buffer)
        except Exception, ex:
            from traceback import print_exception
            print_exception(*sys.exc_info())
            self.error(0, "Error connecting to server" + str(ex))
            return
        self.error(0)

    def clear_buffer(self):
        self.buffer.clear()
        self.UpdateExperiments()

    def UpdateExperiments(self):
        self.chipsl = []
        self.experimentsWidget.clear()
        self.items = []

        self.progressBarInit()

        if not self.dbc:
            self.Connect()
 
        chips = self.dbc.list()
        annots = self.dbc.annotations(chips)
        elements = []

        pos = 0

        for chip,annot in zip(chips, annots):
            pos += 1
            d = defaultdict(lambda: "?", annot)
            print d.items()
            elements.append([d["species"], d["strain"], d["genotype"], d["replicate"], d["tp"], d["treatment"], d["growth"], chip])
            self.progressBarSet((100.0 * pos) / len(chips))
            self.items.append(MyTreeWidgetItem(self.experimentsWidget, elements[-1]))

        for i in range(7):
            self.experimentsWidget.resizeColumnToContents(i)


        self.progressBarFinished()

    def SearchUpdate(self, string=""):
        for item in self.items:
            item.setHidden(not all(s in item for s in self.searchString.split()))

    def Commit(self):
        if not self.dbc:
            self.Connect()
        allTables = []

        import time
        start = time.time()

        pb = OWGUI.ProgressBar(self, iterations=1000)

        table = None

        ids = []
        for item in self.experimentsWidget.selectedItems():
            ids += str(item.text(7)).split(",")

        table = self.dbc.get_data(ids=ids, callback=pb.advance, exclude_constant_labels=self.excludeconstant)

        end = int(time.time()-start)
        
        pb.finish()

        #self.send("Example table", None)
        self.send("Example table", table)

if __name__ == "__main__":
    app  = QApplication(sys.argv)
##    from pywin.debugger import set_trace
##    set_trace()
    w = OWDicty()
    w.show()
    app.exec_()
    w.saveSettings()
            
        
