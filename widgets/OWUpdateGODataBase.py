"""
<name>Update GO Database</name>
<description>Update GO Database</description>
<contact>Ales Erjavec</contact>
<priority>101</priority>
"""

import go
import sys
import OWGUI
from OWWidget import OWWidget

from qt import *

class OWUpdateGODataBase(OWWidget):
    settingsList = ["dataDir", "aspect", "organism"] 
    def __init__(self, parent=None, signalManager=None, name="Update GO data base", **kwds):
        OWWidget.__init__(self, parent, signalManager, name, **kwds)
        self.dataDir = go.getDataDir()
        self.aspect = 0
        self.organism = 0
        self.space = box = OWGUI.widgetBox(self.controlArea)
        #self.layout = QVBoxLayout(self, 4)
        #self.layout.addWidget(self.space)
        OWGUI.lineEdit(box, self, "dataDir", "Custom Data Base Location", valueType = str)
        #self.aspectnSpin = OWGUI.spinBox(self, self, "aspect", "Aspect")
        self.organsmNames = go.listOrganisms()
        self.organismCombo = OWGUI.comboBox(box, self, "organism", "Organism", items = self.organsmNames)
        OWGUI.button(box, self, "Update Annotation", callback=self.UpdateAnnotation)
        OWGUI.button(box, self, "Updata Ontology", callback=self.UpdateOntology)
        self.loadSettings()
        go.setDataDir(self.dataDir)
        #self.UpdateAspectCombo()
        self.UpdateOrganismCombo()
        self.resize(100, 100)
        #self.updateGUI()

    def UpdateAspectCombo(self):
        self.aspectnCombo.clear()
        self.aspectnCombo.insertStrList(go.listDownloadedAspects())

    def UpdateOrganismCombo(self):
        self.organsmNames = go.listOrganisms()
        self.organismCombo.clear()
        self.organismCombo.insertStrList(self.organsmNames)

    def UpdateAnnotation(self):
        go.setDataDir(self.dataDir)
        org =  self.organsmNames[self.organism]
        self.space.setDisabled(True)
        self.progressBarInit()
        print org
        go.downloadAnnotation(org, self.progressBarSet)
        self.progressBarFinished()
        self.space.setDisabled(False)

    def UpdateOntology(self):
        go.setDataDir(self.dataDir)
        self.space.setDisabled(True)
        self.progressBarInit()
        go.downloadGO(self.progressBarSet)
        self.progressBarFinished()
        self.space.setDisabled(False)

    def close(self, destroy):
        self.emit(PYSIGNAL("closed()"), ())
        self.saveSettings()
        return OWWidget.close(self, destroy)
        

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWUpdateGODataBase()
    app.setMainWidget(w)
    w.show()
    app.exec_loop()
    w.saveSettings()