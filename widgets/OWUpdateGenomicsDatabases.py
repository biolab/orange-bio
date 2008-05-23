"""
<name>Update Genomics Databases</name>
<description>Central widget for updating genomics databases</description>
<contact>Ales Erjavec</contact>
<priority>101</priority>
<icon>icons/UpdateDatabases.png</icon>
"""

import go
import obiKEGG
import obiData
import sys
import OWGUI
from OWWidget import OWWidget

from qt import *

class OWUpdateGenomicsDatabases(OWWidget):
    settingsList = ["goDataDir", "goAspect", "goOrganism"] 
    def __init__(self, parent=None, signalManager=None, name="Update Genomics Databases", **kwds):
        OWWidget.__init__(self, parent, signalManager, name, **kwds)
        self.goDataDir = go.getDataDir()
        self.goAspect = 0
        self.goOrganism = 0
        self.keggDataDir = obiKEGG.default_database_path
        self.keggOrganism = 0
        self.tabWidget = QTabWidget(self.controlArea)
        ##GO
        ###############
        self.goTab = box = OWGUI.widgetBox(self.controlArea)
        self.tabWidget.addTab(self.goTab, "GO")
        OWGUI.lineEdit(box, self, "goDataDir", "Custom Data Base Location", valueType = str)
        self.goOrgansmNames = go.listOrganisms()
        self.goOrganismCombo = OWGUI.comboBox(box, self, "goOrganism", "Organism", items = self.goOrgansmNames)
        OWGUI.button(box, self, "Update Annotation", callback=self.UpdateAnnotation)
        OWGUI.button(box, self, "Updata Ontology", callback=self.UpdateOntology)
        OWGUI.rubber(box)
        ##KEGG
        ##############
        self.keggTab = box = OWGUI.widgetBox(self.controlArea)
        self.tabWidget.addTab(self.keggTab, "KEGG")
        OWGUI.lineEdit(box, self, "keggDataDir", "Custom Data Base Location", valueType = str)
        self.keggOrganisms = obiKEGG.KEGGInterfaceLocal().list_organisms().items()
        self.keggOrganisms.sort()
        cb = OWGUI.comboBox(box, self, "keggOrganism", "Organisms", items=["%s: %s" % t for t in self.keggOrganisms])
        cb.setMaximumWidth(250)
        OWGUI.button(box, self, "Update organism data", callback=self.UpdateKEGGOrganism)
        OWGUI.rubber(box)
        #OWGUI.button(box, self, "Update enzyme and compound data", callback=self.UpdateKEGGEnzymeAndCompounds)
        
        self.loadSettings()
        go.setDataDir(self.goDataDir)
        #self.UpdateAspectCombo()
        self.UpdateOrganismCombo()
        self.resize(100, 200)
        #self.updateGUI()

    def UpdateOrganismCombo(self):
        self.goOrgansmNames = go.listOrganisms()
        self.goOrganismCombo.clear()
        self.goOrganismCombo.insertStrList(self.goOrgansmNames)

    def UpdateAnnotation(self):
        go.setDataDir(self.goDataDir)
        org =  self.goOrgansmNames[self.goOrganism]
        self.tabWidget.setDisabled(True)
        self.progressBarInit()
        go.downloadAnnotation(org, self.progressBarSet)
        self.progressBarFinished()
        self.tabWidget.setDisabled(False)

    def UpdateOntology(self):
        go.setDataDir(self.goDataDir)
        self.tabWidget.setDisabled(True)
        self.progressBarInit()
        go.downloadGO(self.progressBarSet)
        self.progressBarFinished()
        self.tabWidget.setDisabled(False)

    def UpdateKEGGOrganism(self):
        org = self.keggOrganisms[self.keggOrganism][0]
        self.tabWidget.setDisabled(True)
        self.progressBarInit()
        obiKEGG.KEGGInterfaceLocal(True, self.keggDataDir, threads=10, download_progress_callback=self.progressBarSet).download_organism_data(org)
        self.progressBarFinished()
        self.tabWidget.setDisabled(False)

    def UpdateKEGGEnzymeAndCompounds(self):
        self.progressBarInit()
        f = obiData.FtpDownloader("ftp.genome.jp", self.keggDataDir, "/pub/kegg/")
        self.progressBarInit()
        f.retrieve("brite/ko/ko00001.kegg", True, progressCallback=self.progressBarSet)
        f.retrieve("ligand/compound/compound",True, progressCallback=self.progressBarSet)
        f.retrieve("ligand/enzyme/enzyme", True, progressCallback=self.progressBarSet)
        self.progressBarFinished()

    def close(self, destroy):
        self.emit(PYSIGNAL("closed()"), ())
        self.saveSettings()
        return OWWidget.close(self, destroy)
        

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWUpdateGenomicsDatabases()
    app.setMainWidget(w)
    w.show()
    app.exec_loop()
    w.saveSettings()