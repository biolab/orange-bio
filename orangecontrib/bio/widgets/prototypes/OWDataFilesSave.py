"""
<name>Save Data Files</name>
<description>Saves data to selected directory.</description>
<icon>icons/DataFilesSave.png</icon>
<priority>1055</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

from __future__ import absolute_import

import os, orange

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .OWDataFiles import DataFiles

class OWDataFilesSave(OWWidget):
    settingsList  = ["recentDirs", "selectedDirName"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "Save Data Files", wantMainArea = 0, resizingEnabled = 1)

        self.inputs = [("Structured Data", DataFiles, self.structuredData, Default)]
        self.outputs = []
        self.dataStructure = None

        # Settings
        self.recentDirs=[]
        self.selectedDirName = "(none)"
        self.loadSettings()

        # GUI
        rfbox = OWGUI.widgetBox(self.controlArea, "Directory", orientation="horizontal", addSpace = True)
        self.dircombo = QComboBox(rfbox)
        rfbox.layout().addWidget(self.dircombo)
        browse = OWGUI.button(rfbox, self, '&...', callback = self.browseDirectory, disabled=0)
        browse.setMaximumWidth(25)
        # info
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace = True)
        self.infoa = OWGUI.widgetLabel(box, 'No data on input.')
        # Output
        box = OWGUI.widgetBox(self.controlArea, "Output", addSpace = True)
        self.save = OWGUI.button(box, self, '&Save', callback = self.saveData, disabled=1)
        self.adjustSize()
           
        # initial settings            
        self.recentDirs=filter(os.path.exists, self.recentDirs)
        self.setDirlist()
        self.dircombo.setCurrentIndex(0)
        self.resize(300,self.height())

    def structuredData(self, data):
        self.dataStructure = data
        self.save.setDisabled(self.dataStructure == None)
        if self.dataStructure == None:
            self.infoa.setText("No data on input.")
        else:
            self.infoa.setText("New data on input.")
            
    def saveData(self):
        # Saves data into root directory
        rootDir = self.recentDirs[self.dircombo.currentIndex()]
        if rootDir == "(none)" or self.dataStructure == None:
            self.infoa.setText("Select a directory first.")
            return
        # count number of files to save
        n = sum([sum([1 for d in ds]) for (sdn, ds) in self.dataStructure])	
        if n == 0:
            self.infoa.setText("No files to save.")
            return
        pbStep = 100./n
        self.progressBarInit()

        for (subDirName, datasets) in self.dataStructure:
            targetDir = os.path.join(rootDir, subDirName)
            if not os.path.exists(targetDir):
                try:
                    os.mkdir(targetDir)
                except:
                    self.infoa.setText("Could not create target directory: " + targetDir)
                    self.error("Could not create target directory: " + targetDir)

            for data in datasets:
                fname = os.path.join(targetDir, data.name + '.tab')
                orange.saveTabDelimited(fname, data)
                self.progressBarAdvance(pbStep)
        self.infoa.setText("Data saved to %s" % rootDir)
        self.progressBarFinished()

    def browseDirectory(self):
        # displays a file dialog and selects a directory
        if len(self.recentDirs):
            startdir=self.recentDirs[0]
        else:
            startdir ="."
        dirname=str(QFileDialog.getExistingDirectory(self, 'Select Directory to Save into', startdir, QFileDialog.ShowDirsOnly))
        if len(dirname):
            self.addDirToList(dirname)

    def setDirlist(self):
        self.dircombo.clear()
        if len(self.recentDirs):
            for dir in self.recentDirs:
                ### leave out the path
                ##(upperdir,dirname)=os.path.split(dir)
                self.dircombo.insertItem(self.dircombo.count(), dir)
        else:
            self.dircombo.insertItem(0, "(none)")

    def addDirToList(self, dir):
        # Add a directory to the start of the file list. 
        # If it exists, move it to the start of the list
        dir = os.path.normpath(dir)
        if dir in self.recentDirs:
            self.recentDirs.remove(dir)
        self.recentDirs.insert(0, str(dir))
        self.setDirlist()
        self.selectedDirName = dir

    def selectDir(self, n):
        # called when user makes a selection from the drop-down menu
        if self.recentDirs:
            self.loadData(self.recentDirs[n])
            self.addDirToList(self.recentDirs[n])
        else:
            self.loadData("(none)")


if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWDataFilesSave()
    ow.show()
    a.exec_()
    ow.saveSettings()
