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
from OWWidget import *

from functools import partial
        
class MultiUpdateItem(QFrame):
    def __init__(self, items=[], parent=None):
        QFrame.__init__(self, parent)
        self.setItems(items)

    def setItems(self, items=[]):
        layout = QGridLayout()
        self.items = items
        for i, (label, time, buttonText, callback) in enumerate(items):
            l1 = QLabel(label)
            l2 = QLabel(time)
            button = QPushButton(buttonText)
            self.connect(button, SIGNAL("clicked()"), callback)
            layout.addWidget(l1, i, 0)
            layout.addWidget(l2, i, 1)
            layout.addWidget(button, i, 2)
        self.setLayout(layout)
        
class KEGGUpdateWidget(QGroupBox):
    def __init__(self, master, parent):
        QGroupBox.__init__(self, "KEGG Update", parent)
        self.master = master
        self._id = 1
        self.setLayout(QVBoxLayout())
        self.update = obiKEGG.Update.getinstance(obiKEGG.default_database_path, self.master.progressBarSet)
##        updatable = dict([(func, (desc, args)) for func, desc, args in self.update.GetUpdatable()])
##        downloadable = dict([(func, (desc, args)) for func, desc, args in self.update.GetDownloadable()])
##        organisms = obiKEGG.KEGGInterfaceLocal(update=False).list_organisms()
        
        ###Organism data update
##        desc, args = updatable.get(obiKEGG.Update.UpdateOrganism, ("", []))
##        args.sort()
##        items = [(arg+" : "+organisms.get(arg, "")+" last updated on: "+self.update.GetLastUpdateTime(obiKEGG.Update.UpdateOrganism, (arg,)).strftime("%Y %B %d %H:%M:%S"),
##                  "Update", partial(self.master._update_helper, self.update.UpdateOrganism, arg)) for arg in args]
        self.updateBox = OWGUI.widgetBox(self, "Update Organisms")
        self.updateArea = QScrollArea()
        self.updateBox.layout().addWidget(self.updateArea)
        self.updateStack = QStackedWidget()
        self.updateArea.setWidget(self.updateStack)
##        self.updateStack.setSizePolicy(QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum))
##        self.updateView.addChild(self.updateStack)
##        self.orgUpdate = MultiUpdateItem(items, self.updateView)
##        self.updateView.addChild(self.orgUpdate)
##        scrollView = QScrollView(box)
##        self.orgUpdate = MultiUpdateItem(items, scrollView)
##        scrollView.addChild(self.orgUpdate)

        ###Organism data download
##        desc, args = downloadable.get(obiKEGG.Update.UpdateOrganism, ("", []))
##        args.sort()
##        items = [(arg+" : "+organisms.get(arg, ""), "Download", partial(self.master._update_helper, self.update.UpdateOrganism, arg)) for arg in args]
        self.downloadBox = OWGUI.widgetBox(self, "Download Organisms")
        self.downloadArea = QScrollArea()
        self.downloadBox.layout().addWidget(self.downloadArea)
        self.downloadStack = QStackedWidget()
        self.downloadArea.setWidget(self.downloadStack)
##        self.orgDownload = MultiUpdateItem(items, self.downloadView)
##        self.downloadView.addChild(self.orgDownload)
##        scrollView = QScrollView(box)
##        self.orgDownload = MultiUpdateItem(items, scrollView)
##        scrollView.addChild(self.orgDownload)

        ###Reference data update
        box = OWGUI.widgetBox(self, "Reference Pathways")
        text = "" #obiKEGG.Update.UpdateReference in updatable and "Update" or "Download"
        self.label = QLabel("")
        self.button = QPushButton(text)
        box.layout().addWidget(self.label)
        box.layout().addWidget(self.button)
        self.connect(self.button, SIGNAL("clicked()"), partial(self._update_helper, self.update.UpdateReference))
        self.SetItems()

    def SetItems(self):
        self._id+=1
        updatable = dict([(func, (desc, args)) for func, desc, args in self.update.GetUpdatable()])
        downloadable = dict([(func, (desc, args)) for func, desc, args in self.update.GetDownloadable()])
        organisms = obiKEGG.KEGGInterfaceLocal(update=False).list_organisms()
        print downloadable
        ###Organism data update
        desc, args = updatable.get(obiKEGG.Update.UpdateOrganism, ("", []))
        args.sort()
        items = [(arg+" : "+organisms.get(arg, ""), "last updated on: "+self.update.GetLastUpdateTime(obiKEGG.Update.UpdateOrganism, (arg,)).strftime("%Y %B %d %H:%M:%S"),
                  "Update", partial(self._update_helper, self.update.UpdateOrganism, arg)) for arg in args]
##        if hasattr(self, "orgUpdate"):
##            self.updateView.removeChild(self.orgUpdate)
        self.orgUpdate = MultiUpdateItem(items, self.updateStack)
##        self.updateView.addChild(self.orgUpdate)
        self.updateStack.addWidget(self.orgUpdate)
        self.updateStack.setCurrentWidget(self.orgUpdate)
        self.updateStack.resize(self.orgUpdate.sizeHint())
##        self.orgUpdate.hide()
##        self.orgUpdate.show()

        ###Organism data download
        desc, args = downloadable.get(obiKEGG.Update.UpdateOrganism, ("", []))
        args.sort()
        items = [(arg+" : "+organisms.get(arg, ""), " ", "Download", partial(self._update_helper, self.update.UpdateOrganism, arg)) for arg in args]
##        if hasattr(self, "orgDownload"):
##            self.downloadView.removeChild(self.orgDownload)
        self.orgDownload = MultiUpdateItem(items, self.downloadStack)
##        self.downloadView.addChild(self.orgDownload)
        self.downloadStack.addWidget(self.orgDownload)
        self.downloadStack.setCurrentWidget(self.orgDownload)
##        self.downloadStack.setChildGeometries()
        self.downloadStack.resize(self.orgDownload.sizeHint())
##        self.orgDownload.hide()
##        self.orgDownload.show()

        ###Reference data update
        text = obiKEGG.Update.UpdateReference in updatable and "Update" or "Download"
        if text=="Update":
            self.label.setText("last updated on: "+self.update.GetLastUpdateTime(obiKEGG.Update.UpdateReference, ()).strftime("%Y %B %d %H:%M:%S"))
        self.button.setText(text)

    def _update_helper(self, func, *args):
        self.master.progressBarInit()
        func(*args)
        self.master.progressBarFinished()
        self.SetItems()

from obiGeneMatch import GeneMatchMk2
import obiGO

class GOUpdateWidget(QGroupBox):
    def __init__(self, master, parent):
        QGroupBox.__init__(self, "GO Update", parent)
        self.master = master
        self.update = obiGO.Update.getinstance(obiGO.getDataDir(), self.master.progressBarSet)
        self.setLayout(QVBoxLayout())
        self._id = 1
        box = OWGUI.widgetBox(self, "Update Annotation")
        area = QScrollArea()
        box.layout().addWidget(area)
        self.updateStack = QStackedWidget()
        area.setWidget(self.updateStack)

        box = OWGUI.widgetBox(self, "Download Annotation")
        area = QScrollArea()
        box.layout().addWidget(area)
        self.downloadStack = QStackedWidget()
        area.setWidget(self.downloadStack)

        box = OWGUI.widgetBox(self, "Ontology")
        self.label = QLabel("")
        self.button = QPushButton("")
        box.layout().addWidget(self.label)
        box.layout().addWidget(self.button)
        self.connect(self.button, SIGNAL("clicked()"), partial(self.update.UpdateOntology))
        self.SetItems()
        
    def SetItems(self):
        self._id+=1
        updatable = dict([(func, (desc, args)) for func, desc, args in self.update.GetUpdatable()])
        downloadable = dict([(func, (desc, args)) for func, desc, args in self.update.GetDownloadable()])
        organisms = obiKEGG.KEGGInterfaceLocal(update=False).list_organisms()

        desc, args = updatable.get(obiGO.Update.UpdateAnnotation, ("", []))
        items = [(arg+" : "+organisms.get(GeneMatchMk2.dbOrgMap.get(arg, ""), ""), "last updated on: "+self.update.GetLastUpdateTime(obiGO.Update.UpdateAnnotation, (arg,)).strftime("%Y %B %d %H:%M:%S"),
                 "Update", partial(self._update_helper, self.update.UpdateAnnotation, arg)) for arg in args]
        self.annoUpdate = MultiUpdateItem(items, self.updateStack)
        self.updateStack.addWidget(self.annoUpdate)
        self.updateStack.setCurrentWidget(self.annoUpdate)
##        self.updateStack.setChildGeometries()
        self.updateStack.resize(self.annoUpdate.sizeHint())

        desc, args = downloadable.get(obiGO.Update.UpdateAnnotation, ("", []))
        items = [(arg+" : "+organisms.get(GeneMatchMk2.dbOrgMap.get(arg, ""), ""), "",
                 "Download", partial(self._update_helper, self.update.UpdateAnnotation, arg)) for arg in args]
        self.annoDownload = MultiUpdateItem(items, self.downloadStack)
        self.downloadStack.addWidget(self.annoDownload)
        self.downloadStack.setCurrentWidget(self.annoDownload)
##        self.downloadStack.setChildGeometries()
        self.downloadStack.resize(self.annoDownload.sizeHint())

        text = obiGO.Update.UpdateOntology in updatable and "Update" or "Download"
        if text=="Update":
            self.label.setText("last updated on: "+self.update.GetLastUpdateTime(obiGO.Update.UpdateOntology, ()).strftime("%Y %B %d %H:%M:%S"))
        self.button.setText(text)

    def _update_helper(self, func, *args):
        self.master.progressBarInit()
        func(*args)
        self.master.progressBarFinished()
        self.SetItems()
        
class OWUpdateGenomicsDatabases(OWWidget):
    def __init__(self, parent=None, signalManager=None, name="Update Genomics Databases", **kwds):
        OWWidget.__init__(self, parent, signalManager, name, **kwds)
##        self.layout = QVBoxLayout(self.mainArea)
##        self.layout.setAutoAdd(True)
        self.tabWidget = OWGUI.tabWidget(self.mainArea)
        self.keggUpdate = KEGGUpdateWidget(self, self.tabWidget)
        OWGUI.createTabPage(self.tabWidget, "KEGG", self.keggUpdate)
##        self.tabWidget.addTab(self.keggUpdate, "KEGG")
        self.goUpdate = GOUpdateWidget(self, self.tabWidget)
        OWGUI.createTabPage(self.tabWidget, "GO", self.goUpdate)
##        self.tabWidget.addTab(self.goUpdate, "GO")

    def saveSettings(self, *args, **kw):
        OWWidget.saveSettings(self, *args, **kw)
        self.keggUpdate.update.shelve.sync()
        self.goUpdate.update.shelve.sync()
        
if __name__ == "__main__":
##    from pywin.debugger import set_trace
    app = QApplication(sys.argv)
##    set_trace()
    w = OWUpdateGenomicsDatabases()
    w.show()
    app.exec_()
    w.saveSettings()