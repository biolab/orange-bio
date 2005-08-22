"""
<name>Save Data Files</name>
<description>Saves data to selected directory.</description>
<category>Genomics</category>
<icon>icons/DataFilesSave.png</icon>
<priority>1055</priority>
"""

from OWWidget import *
import OWGUI
import os, orange
from OWDataFiles import DataFiles


class OWDataFilesSave(OWWidget):
    settingsList  = ["recentDirs", "selectedDirName"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "Save Data Files")

        self.inputs = [("Structured Data", DataFiles, self.structuredData)]
        self.outputs = []

        self.dataStructure = None;

        # Settings
        self.recentDirs=[]
        self.selectedDirName = "(none)"
        self.loadSettings()

        # GUI
        vb = OWGUI.widgetBox(self.space, orientation="horizontal")
        
        rfbox = OWGUI.widgetBox(vb, "Directory", orientation="horizontal")
        self.dircombo = QComboBox(rfbox)
        browse = OWGUI.button(rfbox, self, '&Browse...', callback = self.browseDirectory, disabled=0)

        fbox = OWGUI.widgetBox(vb, "Directory")
        self.save = OWGUI.button(fbox, self, '&Save', callback = self.saveData, disabled=1)
        self.adjustSize()
           
        # initial settings            
        self.recentDirs=filter(os.path.exists, self.recentDirs)
        self.setDirlist()
        self.dircombo.setCurrentItem(0)

    def structuredData(self, data):
        self.dataStructure = data
        self.save.setDisabled(self.dataStructure == None)
            
    # Saves data into root directory
    def saveData(self):
        rootDir = self.recentDirs[self.dircombo.currentItem()]
        if rootDir == "(none)" or self.dataStructure == None:
            return

        # count number of files to save
        n = sum([sum([1 for d in ds]) for (sdn, ds) in self.dataStructure])	
        if n == 0:
            return
        pbStep = 100./n
        self.progressBarInit()

        for (subDirName, datasets) in self.dataStructure:
            targetDir = rootDir + subDirName
            if not os.path.exists(targetDir):
                try:
                    os.mkdir(targetDir)
                except:
                    self.error("Could not create target directory: " + targetDir)

            for data in datasets:
                fname = targetDir + '/' + data.name + '.tab'
                orange.saveTabDelimited(fname, data)
                self.progressBarAdvance(pbStep)
        self.progressBarFinished()

    # displays a file dialog and selects a directory
    def browseDirectory(self):
        if len(self.recentDirs):
            startdir=os.path.split(self.recentDirs[0][:-1])[0]
        else:
            startdir ="."
        dirname=str(QFileDialog.getExistingDirectory(startdir, None, '', 'Select Directory to Save into', 1))
        if len(dirname):
            self.addDirToList(dirname)
            self.saveData()

    def setDirlist(self):
        self.dircombo.clear()
        if len(self.recentDirs):
            for dir in self.recentDirs:
                (upperdir,dirname)=os.path.split(dir[:-1]) #:-1 removes the trailing '\'
                #leave out the path
                self.dircombo.insertItem(dirname)
        else:
            self.dircombo.insertItem("(none)")

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
    ow=OWDataFilesSave()
    a.setMainWidget(ow)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
