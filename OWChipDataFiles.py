"""
<name>Chip Data Files</name>
<description>Loads microarray data sets from designated directory.</description>
<category>Genomics</category>
<icon>icons/ChipDataFiles.png</icon>
<priority>1100</priority>
"""

from OWWidget import *
import OWGUI
import os, os.path, orange

class ChipData(orange.Orange):
    pass

class Mico(orange.Orange):
    pass

class myCheckListItem(QCheckListItem):
    def __init__(self, *args):
        self.callback = None
        QCheckListItem.__init__(self, *args)
        self.setOn(1)

    def stateChange(self, b):
        if self.callback:
            self.callback()

class OWChipDataFiles(OWWidget):
    settingsList  = ["recentDirs","selectedDirName", "applyOnChange"]

    def __init__(self, parent=None):
        OWWidget.__init__(self, parent, 'Chip Data Files')

        self.callbackDeposit = []

        self.inputs = []
        self.outputs = [("Structured Chip Data", ChipData)]

        self.chipdata = None; self.datasets = None
        # Settings
        self.recentDirs=[]
        self.selectedDirName = "None"
        self.applyOnChange = 0
        self.loadSettings()

        # GUI
        # CONTROLS
        box = QHGroupBox("Directory", self.controlArea)
        self.dircombo=QComboBox(box)
        self.dircombo.setMinimumWidth(250)
        button = OWGUI.button(box, self, '...', callback = self.browseDirectory, disabled=0)
        button.setMaximumWidth(25)

        # info
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel('No data loaded.', box)
        self.infob = QLabel('', box)
            
        # Output
        out = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(out, self, 'applyOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(out, self, "Commit", callback=self.sendData, disabled=1)

        # connecting GUI to code
        self.connect(self.dircombo,SIGNAL('activated(int)'),self.selectDir)
        self.resize(500,200)
            
        # LIST VIEW
        self.layout=QVBoxLayout(self.mainArea)
        self.splitter = QSplitter(QSplitter.Vertical, self.mainArea)
        self.layout.add(self.splitter)

        self.tree = QListView(self.splitter)
        self.tree.setAllColumnsShowFocus(1)
        self.tree.addColumn('Directory/Data File')
        self.tree.setColumnWidth(0, 300)
        self.tree.setColumnWidthMode(0, QListView.Manual)
        self.tree.setColumnAlignment(0, QListView.AlignLeft)
            
        self.recentDirs=filter(os.path.exists,self.recentDirs)
        self.setDirlist()
        self.dircombo.setCurrentItem(0)

        if self.recentDirs!=[]:
            self.loadData(self.recentDirs[0])

##        self.commitBtn = OWGUI.button(box, self, "Test", callback=self.test)
            
    def setFileTree(self):
        self.tree.clear()
        self.listitems = []
        for d in self.chipdata:
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
        while dir:
            if dir.isOn():
                files = []
                f = dir.firstChild()
                while f:
                    if f.isOn():
                        files.append(f.data)
                    f = f.nextSibling()
                if len(files):
                    data.append((dir.name, files))
            dir = dir.nextSibling()
        self.send("Structured Chip Data", data)

    # Loads the chip data from a root directory, sends the data to the output channels
    def loadData(self, root):
        self.okToCommit = 0
        if root == "(none)":
            self.send("Structured Chip Data", None)
        chipdata = [] # structured [(dirname0, [d00, d01, ...]), ...]
        datasets = []  # flat list containing all the data sets
        dirs = os.listdir(root)
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
                            print 'eee Exception, file not in appropriate format'
                if len(dirdata):
                    chipdata.append((os.path.split(dirname)[1], dirdata))
                    datasets = datasets + dirdata

        self.commitBtn.setEnabled(len(chipdata))
        if len(chipdata):
            self.chipdata = chipdata
            self.datasets = datasets

##            self.addDirToList(root)
##            self.selectedDirName = root

            self.infoa.setText("Microarray data, %d strains, total of %d data files" % (len(self.chipdata), len(self.datasets)))
            d = self.datasets[0]
            self.infob.setText("Each data file contains %d measurements of %d genes" % (len(d.domain.attributes), len(d)))
            
            self.setFileTree()
            self.okToCommit = 1
            if self.applyOnChange:
                self.sendData()

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

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWChipDataFiles()
    a.setMainWidget(ow)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
