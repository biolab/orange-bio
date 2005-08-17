"""
<name>Data Files</name>
<description>Reads data from designated directory.</description>
<category>Genomics</category>
<icon>icons/ChipDataFiles.png</icon>
<priority>1050</priority>
"""

from OWWidget import *
import OWGUI
import os, os.path, orange


class DataFiles(orange.Orange):
    """Structure for communicating multiple ExampleTables:
    [(name1,[exampleTable1a,exampleTable1b,...]), (name2,[exampleTable2a,...]), ...]
    """
    pass


class ExampleSelection(orange.Orange):
    """Structure for selection of examples: (selectorName, [0,1,...])
    """
    pass


class OWDataFiles(OWWidget):
    settingsList  = ["recentDirs", "selectedDirName", "applyOnChange"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Data Files')

        self.callbackDeposit = []

        self.inputs = []
        self.outputs = [("Structured Data", DataFiles), ("Examples", ExampleTable)]

        self.dataStructure = None;
        self.datasets = None
        self.lastSentIds = []

        # Settings
        self.recentDirs=[]
        self.selectedDirName = "None"
        self.applyOnChange = 0
        self.loadSettings()

        # GUI
        self.mainArea.setFixedWidth(0)
        ca=QFrame(self.controlArea)
        gl=QGridLayout(ca,4,1,5)
        # CONTROLS
        box = QHGroupBox("Directory", ca)
        gl.addWidget(box,0,0)
        self.dircombo=QComboBox(box)
        button = OWGUI.button(box, self, '...', callback = self.browseDirectory, disabled=0)
        button.setMaximumWidth(25)
        # connecting GUI to code
        self.connect(self.dircombo,SIGNAL('activated(int)'),self.selectDir)

        # info
        box = QVGroupBox("Info", ca)
        gl.addWidget(box,1,0)
        self.infoa = QLabel('No data loaded.', box)
        self.infob = QLabel('', box)
        self.infoc = QLabel('', box)
            
        # LIST VIEW
        frmListView = QFrame(ca)
        gl.addWidget(frmListView,2,0)
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
        gl.addWidget(box,3,0)
        OWGUI.checkBox(box, self, 'applyOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.sendData, disabled=1)

        self.resize(425,425)
            
        # initial settings            
        self.recentDirs=filter(os.path.exists,self.recentDirs)
        self.setDirlist()
        self.dircombo.setCurrentItem(0)

        if self.recentDirs!=[]:
            self.loadData(self.recentDirs[0])

##        self.commitBtn = OWGUI.button(box, self, "Test", callback=self.test)
            
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

    # Loads the chip data from a root directory, sends the data to the output channels
    def loadData(self, root):
        self.okToCommit = 0
        if root == "(none)":
            self.send("Structured Data", None)
            self.send("Examples", None)
        dataStructure = [] # structured [(dirname0, [d00, d01, ...]), ...]
        datasets = []  # flat list containing all the data sets
        dirs = os.listdir(root)
        lenDirs = len(dirs)
        if lenDirs:
            self.progressBarInit()
            pbStep = 100./lenDirs
        for d in dirs:
            dirname = root+'\\'+d
            if os.path.isdir(dirname):
                dirdata = []   
                files  = os.listdir(dirname)
                for f in files:
                    name, ext = os.path.splitext(f)
                    if ext in ['.tab', '.txt', '.data']:
                        try:
                            data = None
                            data = orange.ExampleTable(dirname+'\\'+f)
                            data.name = name
                            dirdata.append(data)
                        except orange.KernelException:
                            print 'Warning: file %s\\%s not in appropriate format' %(dirname, f)
                if len(dirdata):
                    dataStructure.append((os.path.split(dirname)[1], dirdata))
                    datasets = datasets + dirdata
            self.progressBarAdvance(pbStep)
        if lenDirs:
            self.progressBarFinished()
        # enable commit, sumarize the data        
        self.commitBtn.setEnabled(len(dataStructure))
        if len(dataStructure):
            self.dataStructure = dataStructure
            self.datasets = datasets
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

            # read data
            self.setFileTree()
            self.okToCommit = 1
            if self.applyOnChange:
                self.sendData()
        else:
            self.infoa.setText('No data on input.')
            self.infob.setText('')
            self.infoc.setText('')


    # displays a file dialog and selects a directory
    def browseDirectory(self):
        if len(self.recentDirs):
            startdir=os.path.split(self.recentDirs[0][:-1])[0]
        else:
            startdir ="."
        dirname=str(QFileDialog.getExistingDirectory(startdir, None, '', 'Microarray Data Directory', 1))
        if len(dirname):
            self.loadData(str(dirname))
            self.addDirToList(dirname) # XXX do this only if loadData successfull

    def setDirlist(self):
        self.dircombo.clear()
        if len(self.recentDirs):
            for dir in self.recentDirs:
                (upperdir,dirname)=os.path.split(dir[:-1]) #:-1 removes the trailing '\'
                #leave out the path
                self.dircombo.insertItem(dirname)
        else:
            self.dircombo.insertItem("(none)")
        self.dircombo.adjustSize() #doesn't work properly :(

    def addDirToList(self, dir):
        # Add a directory to the start of the file list. 
        # If it exists, move it to the start of the list
        if dir in self.recentDirs:
            self.recentDirs.remove(dir)
        self.recentDirs.insert(0, str(dir))
        self.setDirlist()
        self.selectedDirName = dir

    # called when user makes a selection from the drop-down menu
    def selectDir(self, n):
        if self.recentDirs:
            self.loadData(self.recentDirs[n])
            self.addDirToList(self.recentDirs[n])
        else:
            self.loadData("(none)")


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
    ow=OWDataFiles()
    a.setMainWidget(ow)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
