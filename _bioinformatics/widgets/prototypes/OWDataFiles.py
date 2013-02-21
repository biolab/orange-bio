"""
<name>Data Files</name>
<description>Reads data from designated directory.</description>
<icon>icons/ChipDataFiles.png</icon>
<priority>1050</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

import os, os.path

import orange
from Orange.OrangeWidgets.OWWidget import *
from Orange.OrangeWidgets import OWGUI

import warnings
warnings.filterwarnings("ignore", "'strain'", orange.AttributeWarning)

class DataFiles(orange.OrangeBase):
    """Structure for communicating multiple ExampleTables:
    [(name1,[exampleTable1a,exampleTable1b,...]), (name2,[exampleTable2a,...]), ...]
    """
    pass


class ExampleSelection(orange.OrangeBase):
    """Structure for selection of examples: (selectorName, [0,1,...])
    """
    pass


class OWDataFiles(OWWidget):
    settingsList  = ["recentDirs", "selectedDirName", "applyOnChange"]

    def __init__(self, parent=None, signalManager = None, loaddata=1):
        OWWidget.__init__(self, parent, signalManager, 'Data Files', wantMainArea = 0, resizingEnabled = 1)

        self.callbackDeposit = []

        self.inputs = []
        self.outputs = [("Examples", ExampleTable), ("Structured Data", DataFiles)]

        self.dataStructure = []
        self.datasets = None
        self.lastSentIds = []

        # Settings
        self.recentDirs=[]  
        self.selectedDirName = "None"
        self.applyOnChange = 0
        self.loadSettings()

        # CONTROLS
        box = OWGUI.widgetBox(self.controlArea, "Directory", addSpace = True, orientation=0)
        self.dircombo=QComboBox(box)
        box.layout().addWidget(self.dircombo)
        button = OWGUI.button(box, self, '...', callback = self.browseDirectory, disabled=0)
        button.setMaximumWidth(25)
        # connecting GUI to code
        self.connect(self.dircombo,SIGNAL('activated(int)'),self.selectDir)

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

        # initial settings            
        self.recentDirs=filter(os.path.exists,self.recentDirs)
        self.setDirlist()
        self.dircombo.setCurrentIndex(0)
        if self.recentDirs!=[] and loaddata:
            self.loadData(self.recentDirs[0])
            
    def setFileTree(self):
        self.disconnect(self.tree,SIGNAL('itemSelectionChanged()'),self.selectionChanged)
        self.tree.clear()
        self.listitems = []
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
        if self.applyOnChange:
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
        #print "DATAS", data
        self.send("Structured Data", data)

    # Loads the data from a root directory, sends the data to the output channels
    def loadData(self, root):
        dataStructure = [] # structured [(dirname0, [d00, d01, ...]), ...]
        datasets = []  # flat list containing all the data sets
        dirs = os.listdir(root)
        lenDirs = len(dirs)
        if lenDirs:
            self.progressBarInit()
            pbStep = 100./lenDirs
        for d in dirs:
            dirname = os.path.join(root, d)
            if os.path.isdir(dirname):
                dirdata = []   
                files  = os.listdir(dirname)
                for f in files:
                    name, ext = os.path.splitext(f)
                    if ext in ['.tab', '.txt', '.data']:
                        try:
                            data = None
                            data = orange.ExampleTable(os.path.join(dirname,f))
                            data.name = name
                            data.strain = os.path.basename(dirname)
                            dirdata.append(data)
                        except orange.KernelException:
                            print 'Warning: file %s\\%s not in appropriate format' %(dirname, f)
                if len(dirdata):
                    dataStructure.append((os.path.split(dirname)[1], dirdata))
                    datasets = datasets + dirdata
            self.progressBarAdvance(pbStep)
        if lenDirs:
            self.progressBarFinished()

        # enable commit, set file tree
        self.commitBtn.setEnabled(len(dataStructure))
        self.dataStructure = dataStructure
        self.datasets = datasets
        self.setFileTree()

        # set infos (sumarize the data)
        if len(dataStructure):
            numSets = len(self.dataStructure)
            numFiles = len(self.datasets)
            self.infoa.setText("%d set%s, total of %d data file%s." % (numSets, ["", "s"][numSets!=1], numFiles, ["","s"][numFiles!=1]))
            # construct lists that sumarize the data
            numExamplesList = []
            numAttrList = []
            hasClass = []
            attrNameList = []
            for et in datasets:
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

            # Add a directory to the start of the file list. 
            # If it exists, move it to the start of the list
            if root in self.recentDirs:
                self.recentDirs.remove(root)
            self.recentDirs.insert(0, root)
            self.setDirlist()
            self.selectedDirName = root
                
        else:
            self.infoa.setText('No data on input.')
            self.infob.setText('')
            self.infoc.setText('')

        # read data
        if self.applyOnChange:
            self.sendData()

    # displays a file dialog and selects a directory
    def browseDirectory(self):
        if len(self.recentDirs):
            startdir=self.recentDirs[0]
        else:
            startdir ="."
        dirname=str(QFileDialog.getExistingDirectory(self, 'Data Directory', startdir, QFileDialog.ShowDirsOnly))
        if len(dirname):
            self.loadData(str(dirname))

    def setDirlist(self):
        self.dircombo.clear()
        if len(self.recentDirs):
            for dir in self.recentDirs:
                (upperdir,dirname)=os.path.split(dir)
                # leave out the path
                self.dircombo.insertItem(self.dircombo.count(), dirname)
        else:
            self.dircombo.insertItem(0,"(none)")

    # called when user makes a selection from the drop-down menu
    def selectDir(self, n):
        if self.recentDirs:
            self.loadData(self.recentDirs[n])
        else:
            self.loadData("(none)")



if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWDataFiles()
    ow.show()
    ow.loadData("/home/marko/anova/smallchipdata")
    a.exec_()
    ow.saveSettings()
