"""
<name>PIPA Database</name>
<description>Access data from PIPA RNA-Seq database.</description>
<icon>icons/PIPA.png</icon>
<priority>30</priority>
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

CACHED_COLOR = Qt.darkGreen

class OWPIPA(OWWidget):
    settingsList = [ "platform", "selectedExperiments", "server", "buffertime", "excludeconstant", "username", "password","joinreplicates" ]
    def __init__(self, parent=None, signalManager=None, name="PIPA database"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.outputs = [("Example table", ExampleTable)]

        self.platform = None
        self.username = ""
        self.password = ""

        self.selectedExperiments = []
        self.buffer = obiDicty.BufferSQLite(bufferfile)

        self.searchString = ""
        self.excludeconstant = False
        self.joinreplicates = False

        self.chips = []
        self.annots = []
        

        self.controlArea.setMaximumWidth(250)
        self.controlArea.setMinimumWidth(250)

        OWGUI.button(self.controlArea, self, "Reload", callback=self.Reload)
        OWGUI.button(self.controlArea, self, "Clear cache", callback=self.clear_cache)

        OWGUI.rubber(self.controlArea)

        OWGUI.checkBox(self.controlArea, self, "excludeconstant", "Exclude labels with constant values" )
        OWGUI.checkBox(self.controlArea, self, "joinreplicates", "Average replicates (use median)" )

        OWGUI.button(self.controlArea, self, "&Commit", callback=self.Commit)

        OWGUI.rubber(self.controlArea)
        OWGUI.rubber(self.controlArea)
        OWGUI.rubber(self.controlArea)
        OWGUI.rubber(self.controlArea)

        box  = OWGUI.widgetBox(self.controlArea, "Authentication")


        OWGUI.lineEdit(box, self, "username", "Username:", labelWidth=100, orientation='horizontal', callback=self.ConnectAndUpdate)
        OWGUI.lineEdit(box, self, "password", "Password:", labelWidth=100, orientation='horizontal', callback=self.ConnectAndUpdate)



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

        QTimer.singleShot(100, self.UpdateExperiments)        

        self.resize(800, 600)

    def __updateSelectionList(self, oldList, oldSelection, newList):
        oldList = [oldList[i] for i in oldSelection]
        return [ i for i, new in enumerate(newList) if new in oldList]
    
    def ConnectAndUpdate(self):
        self.Connect()
        self.UpdateExperiments(reload=True)

    def Connect(self):
        try:
            #obiDicty.verbose = 1
            def en(x):
                return x if len(x) else None
            self.dbc = obiDicty.PIPA(buffer=self.buffer, username=en(self.username), password=self.password)
        except Exception, ex:
            from traceback import print_exception
            print_exception(*sys.exc_info())
            self.error(0, "Error connecting to server" + str(ex))
            return
        self.error(0)

    def Reload(self):
        #self.buffer.clear()
        self.UpdateExperiments(reload=True)

    def clear_cache(self):
        self.buffer.clear()
        self.Reload()

    def UpdateExperiments(self, reload=False):
        self.chipsl = []
        self.experimentsWidget.clear()
        self.items = []

        self.progressBarInit()

        if not self.dbc:
            self.Connect()
 
        #obiDicty.verbose = 1

        chips, annots = [], []
        
        sucind = False #success indicator for database index

        try:
            chips = self.dbc.list(reload=reload)
            annots = self.dbc.annotations(chips, reload=reload)
            sucind = True
        except Exception, ex:
            try:
                chips = self.dbc.list()
                annots = self.dbc.annotations(chips)
                self.warning(0, "Can not access database - using cached data.")
                sucind = True
            except Exception,ex:
                self.error(0, "Can not access database.")

        self.warning(0)
        self.error(0)

        elements = []
        pos = 0

        self.chips = list(chips)
        self.annots = list(annots)

        for chip,annot in zip(self.chips, self.annots):
            pos += 1
            d = defaultdict(lambda: "?", annot)
            elements.append([d["species"], d["strain"], d["genotype"], d["replicate"], d["tp"], d["treatment"], d["growth"], chip])
            self.progressBarSet((100.0 * pos) / len(chips))
            self.items.append(MyTreeWidgetItem(self.experimentsWidget, elements[-1]))

        for i in range(7):
            self.experimentsWidget.resizeColumnToContents(i)

        adic = dict(zip(self.chips, self.annots))
        #which is the ok buffer version
        self.wantbufver = lambda x,ad=adic: defaultdict(lambda: "?", ad[x])["map_stop1"]

        self.UpdateCached()

        self.progressBarFinished()

    def UpdateCached(self):
        if self.wantbufver:
            fn = self.dbc.chips_keynaming()
            for item in self.items:
                color = Qt.black
                c = str(item.text(7))
                if self.dbc.inBuffer(fn(c)) == self.wantbufver(c):
                    color = CACHED_COLOR
                brush = QBrush(color)
                for i in range(item.columnCount()):
                    item.setForeground(i, brush)

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
            ids += [ str(item.text(7)) ]

        table = self.dbc.get_data(ids=ids, callback=pb.advance, exclude_constant_labels=self.excludeconstant, bufver=self.wantbufver)

        if self.joinreplicates:
            table = obiDicty.join_replicates(table, ignorenames=["id", "replicate", "name", "map_stop1"], namefn=None, avg=obiDicty.median)

        end = int(time.time()-start)
        
        pb.finish()

        #self.send("Example table", None)
        self.send("Example table", table)

        self.UpdateCached()

if __name__ == "__main__":
    app  = QApplication(sys.argv)
##    from pywin.debugger import set_trace
##    set_trace()
    w = OWPIPA()
    w.show()
    app.exec_()
    w.saveSettings()
            
        
