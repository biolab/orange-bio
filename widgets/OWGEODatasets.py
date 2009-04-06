"""<name>GEO DataSets</name>
<description>Access to Gene Expression Omnibus data sets.</description>
<priority>251</priority>
<contact>Ales Erjavec (ales.erjavec(@at@)fri.uni-lj.si)</contact>
<icon>icons/GEODataSets.png</icon>
"""

import sys, os, glob
from OWWidget import *
import OWGUI
import obiGEO
import orngServerFiles

class SortableItem(QTreeWidgetItem):
    def __lt__(self ,other):
        widget = self.treeWidget()
        column = widget.sortColumn()
        if column in [2, 3, 4, 5]:
            return int(self.text(column)) < int(other.text(column))
        return QTreeWidgetItem.__lt__(self, other)
    
class LinkItem(QWidget):
    def __init__(self, pubmed_id, parent=None):
        QWidget.__init__(self, parent)
        layout = QHBoxLayout()
        if pubmed_id:
            self.link = QLabel('<a href="http://www.ncbi.nlm.nih.gov/pubmed/%s">%s</a>' % (pubmed_id, pubmed_id), self)
        else:
            self.link = QLabel()
        self.link.setOpenExternalLinks(True)
        layout.addWidget(self.link)
        self.setLayout(layout)
        
class OWGEODatasets(OWWidget):
    settingsList = ["outputRows", "minSamples", "includeIf", "mergeSpots"]

    def __init__(self, parent=None ,signalManager=None, name=" GEO Data sets"):
        OWWidget.__init__(self, parent ,signalManager, name)

        self.outputs = [("Example Table", ExampleTable)]

        ## Settings
        self.selectedSubsets = []
        self.sampleSubsets = []
        self.includeIf = False
        self.minSamples = 3
        self.autoCommit = False
        self.outputRows = 0
        self.mergeSpots = True
        self.filterString = ""

        self.loadSettings()

        ## GUI
        self.infoBox = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Info"), "\n")
        box = OWGUI.widgetBox(self.controlArea, "Sample Subset")
        OWGUI.listBox(box, self, "selectedSubsets", "sampleSubsets", selectionMode=QListWidget.ExtendedSelection)
##        OWGUI.button(box, self, "Clear selection", callback=self.clearSubsetSelection)
##        c = OWGUI.checkBox(box, self, "includeIf", "Include if at least", callback=self.commitIf)
##        OWGUI.spin(OWGUI.indentedBox(box), self, "minSamples", 2, 100, posttext="samples", callback=self.commitIf)

        box = OWGUI.widgetBox(self.controlArea, "Output")
        OWGUI.radioButtonsInBox(box, self, "outputRows", ["Genes or spots", "Samples"], "Rows") ##, callback=self.commitIf)
        OWGUI.checkBox(box, self, "mergeSpots", "Merge spots of same gene") ##, callback=self.commitIf)

        box = OWGUI.widgetBox(self.controlArea, "Output")
        OWGUI.button(box, self, "Commit", callback=self.commit)
##        OWGUI.checkBox(box, self, "autoCommit", "Commit automatically")
        OWGUI.rubber(self.controlArea)

        OWGUI.lineEdit(self.mainArea, self, "filterString", "Filter", callbackOnType=True, callback=self.filter)
        self.treeWidget = QTreeWidget(self.mainArea)
        self.treeWidget.setHeaderLabels(["ID", "Organism", "Samples", "Features", "Genes", "Subsets", "PubMedID"])
        self.treeWidget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSortingEnabled(True)
        self.mainArea.layout().addWidget(self.treeWidget)
        self.connect(self.treeWidget, SIGNAL("itemSelectionChanged ()"), self.updateSelection)
##        self.connect(self.treeWidget, SIGNAL("currentItemChanged(QTreeWidgetItem*, QTreeWidgetItem*))"), self.updateSelection)
        self.infoGDS = OWGUI.widgetLabel(OWGUI.widgetBox(self.mainArea, "Description"), "")
        self.infoGDS.setWordWrap(True)

        QTimer.singleShot(50, self.updateTable)
        self.resize(700, 500)

    def updateInfo(self):
        gds_info = obiGEO.GDSInfo()
        self.infoBox.setText("%i datasets\n%i datasets cached" %(len(gds_info), len(glob.glob(orngServerFiles.localpath("GEO") + "/GDS*"))))
        
    def updateTable(self):
        self.treeWidget.clear()
        self.treeItems = []
        info = obiGEO.GDSInfo()
        self.progressBarInit()
        milestones = set(range(0, len(info), max(len(info)/100, 1)))
        for i, (name, gds) in enumerate(info.items()):
            item = SortableItem(None, [gds["dataset_id"], gds["platform_organism"], str(len(gds["samples"])), str(gds["feature_count"]),
                                                     str(gds["gene_count"]), str(len(gds["subsets"])), ""])
            item.link = LinkItem(gds.get("pubmed_id"), self.treeWidget)
##            self.treeWidget.setItemWidget(item, 6, link)
            item.gdsName = name
            item.gds = gds
            self.treeItems.append(item)
            if i in milestones:
                self.progressBarSet(100.0*i/len(info)/2)

        self.treeWidget.addTopLevelItems(self.treeItems)
        for i, item in enumerate(self.treeItems):
            self.treeWidget.setItemWidget(item, 6, item.link)
            if i in milestones:
                self.progressBarSet(50.0 + 100.0*i/len(self.treeItems)/2)
        self.progressBarFinished()                

        self.updateInfo()

    def updateSelection(self):
        current = self.treeWidget.selectedItems()
        if current:
            self.currentItem = current[0]
            self.setSubsets(current[0].gds)
            self.infoGDS.setText(current[0].gds.get("description", ""))
        else:
            self.currentItem = None
        
    def setSubsets(self, gds):
        self.sampleSubsets = ["%s (%d)" % (s["description"], len(s["sample_id"])) for s in gds["subsets"]]

    def clearSubsetSelection(self):
        pass

    def filter(self):
        filterStrings = self.filterString.lower().split()
        searchKeys = ["dataset_id", "platform_organism", "description"]
        for item in self.treeItems:
            item.setHidden(not all([any([s in unicode(item.gds.get(key, "").lower(), errors="ignore") for key in searchKeys]) for s in filterStrings]))

    def commit(self):
        if self.currentItem:
            classes = [s["description"] for s in self.currentItem.gds["subsets"]]
            classes = [classes[i] for i in self.selectedSubsets] or None
            self.progressBarInit()
            self.progressBarSet(10)
            gds = obiGEO.GDS(self.currentItem.gdsName)
            data = gds.getdata(report_genes=self.mergeSpots, transpose=self.outputRows, classes=classes)
            self.progressBarFinished()
            self.send("Example Table", data)
        else:
            pass

    def commitIf(self):
        if self.commitIf:
            self.commit()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWGEODatasets()
    w.show()
    app.exec_()
    w.saveSettings()
