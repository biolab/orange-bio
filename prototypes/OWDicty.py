"""
<name>Dicty database</name>
<description>Widget interface to dicty experaments</description>
"""

from OWWidget import *
import obiDicty
import OWGUI
import sys

class ShadedItemDelegate(QItemDelegate):
    def paint(self, painter, option, index):
        self.drawBackground(painter, option, index)
        value, ok = index.data(Qt.DisplayRole).toString()
class OWDicty(OWWidget):
    settingsList = ["typeIndex", "joinSelected", "separateSelected"]
    def __init__(self, parent=None, signalManager=None, name="Dicty database"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.serverToken = ""
        self.typeIndex = 0
        self.typesList = [n for n, long in obiDicty.DatabaseConnection.obidPairs]

        self.platformIndex = 0
        self.platformList = []
        
        self.sampleList = []
        self.sampleSelected = []
        
        self.joinList = [n for n, long in obiDicty.DatabaseConnection.aoidPairs]
        self.joinSelected = []
        
        self.separateList = [n for n, long in obiDicty.DatabaseConnection.aoidPairs]
        self.separateSelected = []
        
        self.loadSettings()        
        
        OWGUI.lineEdit(self.controlArea, self, "serverToken", box="Server Token", callback=self.UpdateControls)
        self.typeCombo = OWGUI.comboBox(self.controlArea, self, "typeIndex", box="Type", items=self.typesList, sendSelectedValue=True, callback=self.UpdateControls)
        self.platformCombo = OWGUI.comboBox(self.controlArea, self, "platformIndex", box="Platform", items=self.platformList, sendSelectedValue=True, callback=self.UpdateControls)
        self.sampleListBox = OWGUI.listBox(self.controlArea, self, "sampleSelected", "sampleList", box="Samples", selectionMode=QListWidget.ExtendedSelection, callback=self.UpdateControls)
        self.joinListBox = OWGUI.listBox(self.mainArea, self, "joinSelected", "joinList", box="Join", selectionMode=QListWidget.ExtendedSelection, callback=self.UpdateControls)
        self.separateListBox = OWGUI.listBox(self.mainArea, self, "separateSelected", "separateList", box="Separate", selectionMode=QListWidget.ExtendedSelection, callback=self.UpdateControls)
        
        OWGUI.button(self.controlArea, self, "&Preview", callback=self.ShowPreview)
        OWGUI.button(self.controlArea, self, "&Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)

        self.dbc = None
##        self.UpdateControls()

    def __updateSelectionList(self, oldList, oldSelection, newList):
        oldList = [oldList[i] for i in oldSelection]
        return [ i for i, new in enumerate(newList) if new in oldList]

    def GetOptions(self):
        print self.typeIndex, self.platformIndex, self.sampleSelected, self.joinSelected, self.separateSelected
        return self.typeIndex, self.platformIndex, self.sampleSelected, self.joinSelected, self.separateSelected

    def UpdateControls(self):
        if not self.dbc:
            try:
                self.dbc = obiDicty.DatabaseConnection("http://purple.bioch.bcm.tmc.edu/~anup/index.php?token="+self.serverToken+"&")
            except Exception, ex:
                from traceback import print_exception
                print_exception(*sys.exc_info())
                self.error(0, "Error connecting to server")
                return
        else:
            self.error(0)
        type, platform, sample, join, separate = self.GetOptions()
        try:
            opt = self.dbc.annotationOptions(type=type, platform=platform) #, sample=list(sample), join=list(join), separate=list(separate))
        except Exception, ex:
            from traceback import print_exception
            print_exception(*sys.exc_info())
            self.error(1, "Error retrieving data from server")
            return
        else:
            self.error(1)
        self.platformCombo.clear()
        self.platformCombo.addItems(opt.get("platform", []))
        self.sampleList = opt.get("sample", [])


    def ShowPreview(self):
        pass

    def Commit(self):
        pass

if __name__ == "__main__":
    app  = QApplication(sys.argv)
    w = OWDicty()
    w.show()
    app.exec_()
    w.saveSettings()
            
        