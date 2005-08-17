"""
<name>Data Files Selector</name>
<description>Selects a subset of data files.</description>
<category>Genomics</category>
<icon>icons/ChipDataFiles.png</icon>
<priority>1060</priority>
"""

from OWWidget import *
import OWGUI
from OWDataFiles import DataFiles, ExampleSelection


class OWDataFilesSelector(OWWidget):
    settingsList  = ["applyOnChange"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Data Files Selector')

        self.callbackDeposit = []

        self.inputs = [("Structured Data", DataFiles, self.onDataInput, 1)]
        self.outputs = [("Structured Data", DataFiles), ("Examples", ExampleTable)]

        self.dataStructure = None;
        self.datasets = None
        self.lastSentIds = []

        # Settings
        self.applyOnChange = 0
        self.loadSettings()

        # GUI
        self.mainArea.setFixedWidth(0)
        ca=QFrame(self.controlArea)
        gl=QGridLayout(ca,3,1,5)

        # info
        box = QVGroupBox("Info", ca)
        gl.addWidget(box,0,0)
        self.infoa = QLabel('No data loaded.', box)
        self.infob = QLabel('', box)
        self.infoc = QLabel('', box)
            
        # LIST VIEW
        frmListView = QFrame(ca)
        gl.addWidget(frmListView,1,0)
        self.layout=QVBoxLayout(frmListView)
        self.splitter = QSplitter(QSplitter.Vertical, frmListView)
        self.layout.add(self.splitter)
        self.tree = QListView(self.splitter)
        self.tree.setAllColumnsShowFocus(1)
        self.tree.addColumn('Directory/Data File')
        self.tree.setColumnWidth(0, 379)
        self.tree.setColumnWidthMode(0, QListView.Manual)
        self.tree.setColumnAlignment(0, QListView.AlignLeft)

        # Output
        box = QVGroupBox("Output", ca)
        gl.addWidget(box,2,0)
        OWGUI.checkBox(box, self, 'applyOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.sendData, disabled=1)

        self.resize(425,425)


    def onDataInput(self, dataStructure):
        self.dataStructure = dataStructure
        self.datasets = []
        if dataStructure and len(dataStructure):
            for name, etList in dataStructure:
                for et in etList:
                    self.datasets.append(et)
            # enable commit, sumarize the data        
            self.commitBtn.setEnabled(len(dataStructure))
            numSets = len(self.dataStructure)
            numFiles = len(self.datasets)
            self.infoa.setText("%d set%s, total of %d data file%s." % (numSets, ["", "s"][numSets!=1], numFiles, ["","s"][numFiles!=1]))
            # construct lists that sumarize the data
            numExamplesList = []
            numAttrList = []
            hasClass = []
            attrNameList = []
            for et in self.datasets:
                numExamplesList.append(len(et))
                numAttrList.append(len(et.domain.attributes))
                hasClass.append(et.domain.classVar != None)
                for attr in et.domain.attributes:
                    if attr.name not in attrNameList:
                        attrNameList.append(attr.name)
            # report the number of attributes/class
            if len(numAttrList):
                minNumAttr = min(numAttrList)
                maxNumAttr = max(numAttrList)
                if minNumAttr != maxNumAttr:
                    infob = "From %d to %d attribute%s (%d in total)" % (minNumAttr, maxNumAttr, ["","s"][maxNumAttr!=1], len(attrNameList))
                else:
                    infob = "%d attribute%s" % (maxNumAttr, ["","s"][maxNumAttr!=1])
            else:
                infob = "No attributes"
            if sum(hasClass) == len(hasClass):
                infob += ", all files contain class variable."
            elif sum(hasClass) == 0:
                infob += ", no class variable."
            else:
                infob += ", some files contain class variable."
            self.infob.setText(infob)
            # report the number of examples
            if len(numExamplesList):
                infoc = "Files contain "
                minNumE = min(numExamplesList)
                maxNumE = max(numExamplesList)
                if minNumE == maxNumE:
                    infoc += "%d example%s, " % (maxNumE, ["","s"][maxNumE!=1])
                else:
                    infoc += "from %d to %d example%s, " % (minNumE, maxNumE, ["","s"][maxNumE!=1])
                infoc += "%d in total." % sum(numExamplesList)
            else:
                infoc = "Files contain no examples."
            self.infoc.setText(infoc)

            # read data
            self.setFileTree()
            self.okToCommit = 1
            if self.applyOnChange:
                self.sendData()
        else:
            self.infoa.setText('No data on input.')
            self.infob.setText('')
            self.infoc.setText('')


    def setFileTree(self):
        self.tree.clear()
        self.listitems = []
        for d in self.dataStructure:
            (dirname, files) = d
            diritem = myCheckListItem(self.tree, dirname, QCheckListItem.CheckBox)
            diritem.callback = self.selectionChanged
            self.listitems.append(diritem)
            diritem.setOpen(1)
            diritem.name = dirname
            for f in files:
                item = myCheckListItem(diritem, f.name, QCheckListItem.CheckBox)
                item.callback = self.selectionChanged
                self.listitems.append(item)
                item.data = f

    def selectionChanged(self):
        if self.applyOnChange and self.okToCommit:
            self.sendData()

    # checks which data has been selected, builds a chip data structure, and sends it out
    def sendData(self):
        data = []
        dir = self.tree.firstChild()
        # clear what has been previously sent
        for id in self.lastSentIds:
            self.send("Examples", None, id)
        self.lastSentIds = []
        # send new data
        id = 0
        while dir:
            if dir.isOn():
                files = []
                f = dir.firstChild()
                while f:
                    if f.isOn():
                        files.append(f.data)
                        self.send("Examples", f.data, id)
                        self.lastSentIds.append(id)
                        id += 1
                    f = f.nextSibling()
##                if len(files):
##                    data.append((dir.name, files))
                # it should also be possible to send out (dir.name, [])
                data.append((dir.name, files))
            dir = dir.nextSibling()
        self.send("Structured Data", data)


class myCheckListItem(QCheckListItem):
    def __init__(self, *args):
        self.callback = None
        QCheckListItem.__init__(self, *args)
        self.setOn(1)

    def stateChange(self, b):
        if self.callback:
            self.callback()


if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWDataFilesSelector()
    a.setMainWidget(ow)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
