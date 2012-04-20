"""
<name>Data Files Selector</name>
<description>Selects a subset of data files.</description>
<icon>icons/DataFilesSelector.png</icon>
<priority>1060</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

from __future__ import absolute_import

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .OWDataFiles import DataFiles, ExampleSelection

class OWDataFilesSelector(OWWidget):
    settingsList  = ["applyOnChange"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Data Files Selector', wantMainArea = 0, resizingEnabled = 1)

        self.callbackDeposit = []

        self.inputs = [("Structured Data", DataFiles, self.onDataInput)]
        self.outputs = [("Examples", ExampleTable), ("Structured Data", DataFiles)]

        self.dataStructure = None
        self.datasets = None
        self.lastSentIds = []

        # Settings
        self.applyOnChange = 0
        self.loadSettings()

        # GUI
        # info
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace = True)
        self.infoa = OWGUI.widgetLabel(box, 'No data loaded.')
        self.infob = OWGUI.widgetLabel(box, '')
        self.infoc = OWGUI.widgetLabel(box, '')
            
        # LIST VIEW
        frmListView = OWGUI.widgetBox(self.controlArea, None, addSpace = True)
        self.tree = QTreeWidget(frmListView)
        self.tree.setSelectionMode(QAbstractItemView.MultiSelection)
        self.tree.setHeaderLabel("Directory/Data File")
        frmListView.layout().addWidget(self.tree)
        self.connect(self.tree,SIGNAL('itemSelectionChanged()'),self.selectionChanged)

        # Output
        box = OWGUI.widgetBox(self.controlArea, "Output", addSpace = True)
        OWGUI.checkBox(box, self, 'applyOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.sendData, disabled=1)
        self.resize(300,600)

    def setFileTree(self):
        self.disconnect(self.tree,SIGNAL('itemSelectionChanged()'),self.selectionChanged)
        self.tree.clear()
        self.listitems = []
        if self.dataStructure:
            for d in self.dataStructure:
                (dirname, files) = d
                diritem = QTreeWidgetItem(self.tree, [dirname], QTreeWidgetItem.UserType)
                diritem.setSelected(1)
                self.listitems.append(diritem)
                diritem.setExpanded(1)
                diritem.name = dirname
                for f in files:
                    item = QTreeWidgetItem(diritem, [f.name], QTreeWidgetItem.UserType)
                    item.setSelected(1)
                    self.listitems.append(item)
                    item.data = f
        self.connect(self.tree,SIGNAL('itemSelectionChanged()'),self.selectionChanged)

    def selectionChanged(self):
        if self.applyOnChange and self.okToCommit:
            self.sendData()

    # checks which data has been selected, builds a data structure, and sends it out
    def sendData(self):
        data = []
        # clear what has been previously sent
        for id in self.lastSentIds:
            self.send("Examples", None, id)
        self.lastSentIds = []
        # send new data
        id = 0
        for tIdx in range(self.tree.topLevelItemCount()):
            dir = self.tree.topLevelItem(tIdx)
            if dir in self.tree.selectedItems():
                files = []
                for cIdx in range(dir.childCount()):
                    f = dir.child(cIdx)
                    if f in self.tree.selectedItems():
                        files.append(f.data)
                        self.send("Examples", f.data, id)
                        self.lastSentIds.append(id)
                        id += 1
                data.append((dir.name, files))
        self.send("Structured Data", data)

    def onDataInput(self, dataStructure):
        self.dataStructure = dataStructure
        self.datasets = []
        if dataStructure and len(dataStructure):
            for name, etList in dataStructure:
                for et in etList:
                    self.datasets.append(et)
            # sumarize the data
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
            self.okToCommit = 1
            self.commitBtn.setEnabled(1)
        else:
            self.infoa.setText('No data on input.')
            self.infob.setText('')
            self.infoc.setText('')
            self.okToCommit = 0
            self.commitBtn.setEnabled(0)
        # enable/disable commit
        # read data
        self.setFileTree()
        if self.applyOnChange:
            self.sendData()


if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWDataFilesSelector()
    ow.show()
    a.exec_()
    ow.saveSettings()
