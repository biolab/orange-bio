"""<name>Gene Info</name>
"""

import obiGene, obiTaxonomy
import orange
import orngServerFiles

from OWWidget import *
import OWGUI

class LinkItem(QWidget):
    def __init__(self, gene_id, parent=None):
        QWidget.__init__(self, parent)
        layout = QHBoxLayout()
        self.link = QLabel('<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=%s">%s</a>' % (gene_id, gene_id), self)
        self.link.setOpenExternalLinks(True)
        layout.addWidget(self.link)
        self.setLayout(layout)
            
class LinkItemDelegate(QItemDelegate):
    def sizeHint(self, option, index):
        size = QItemDelegate.sizeHint(self, option, index)
        parent = self.parent()
        item = parent.itemFromIndex(index)
        widget = parent.itemWidget(item, 0)
        if widget:
            size = QSize(size.width(), widget.sizeHint().height())
        return size

class OWGeneInfo(OWWidget):
    settingsList = ["organismIndex", "geneAttr", "useAttr", "autoCommit"]
    def __init__(self, parent=None, signalManager=None, name="Gene info"):
        OWWidget.__init__(self, parent, signalManager, name)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Selected Examples", ExampleTable)]

        self.organismIndex = 0
        self.geneAttr = 0
        self.useAttr = False
        self.autoCommit = False
        self.searchString = ""
        self.loadSettings()
        
        self.infoLabel = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Info"), "")
        self.organisms = [name.split(".")[-2] for name in orngServerFiles.listfiles("NCBI_geneinfo")]
        self.orgaismsComboBox = OWGUI.comboBox(self.controlArea, self, "organismIndex", "Organism", items=[obiTaxonomy.name(id) for id in self.organisms], callback=self.setItems)
        box = OWGUI.widgetBox(self.controlArea, "Gene names")
        self.geneAttrComboBox = OWGUI.comboBox(box, self, "geneAttr", "Gene atttibute", callback=self.setItems)
        c = OWGUI.checkBox(box, self, "useAttr", "Use attribute names", callback=self.setItems, disables=[(-1, self.geneAttrComboBox)])
        self.geneAttrComboBox.setDisabled(bool(self.useAttr))

        box = OWGUI.widgetBox(self.controlArea, "Commit")
        OWGUI.button(box, self, "Commit", callback=self.commit)
        OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.rubber(self.controlArea)

        OWGUI.lineEdit(self.mainArea, self, "searchString", "Filter", callbackOnType=True, callback=self.searchUpdate)
        self.treeWidget = QTreeWidget(self.mainArea)
        self.treeWidget.setHeaderLabels(["NCBI ID", "Symbol", "Locus Tag", "Chromosome", "Description", "Synonyms", "Nomenclature"])
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)
##        self.treeWidget.setItemDelegate(LinkItemDelegate(self.treeWidget))
        self.connect(self.treeWidget, SIGNAL("itemSelectionChanged()"), self.commitIf)
        self.mainArea.layout().addWidget(self.treeWidget)
        
        box = OWGUI.widgetBox(self.mainArea, "", orientation="horizontal")
        OWGUI.button(box, self, "Select Filtered", callback=self.selectFiltered)
        OWGUI.button(box, self, "Clear Selection", callback=self.clearSelection)

        self.geneinfo = []
        
    def setData(self, data=None):
        self.data = data
        if data:
            self.geneAttrComboBox.clear()
            self.attributes = [attr for attr in self.data.domain.variables + self.data.domain.getmetas().values() if attr.varType in [orange.VarTypes.String, orange.VarTypes.Discrete]]
            self.geneAttrComboBox.addItems([attr.name for attr in self.attributes])
            self.geneAttr = min(self.geneAttr, len(self.attributes) - 1)
            self.setItems()
        else:
            self.clear()

    def setItems(self):
        if self.useAttr:
            genes = [attr.name for attr in self.data.domain.attributes]
        else:
            attr = self.attributes[self.geneAttr]
            genes = [str(ex[attr]) for ex in self.data if not ex[attr].isSpecial()]
        info = obiGene.NCBIGeneInfo(self.organisms[self.organismIndex])
        self.geneinfo = geneinfo = [(gene, info.get_info(gene, None)) for gene in genes]
        print genes[:10]
        self.treeWidget.clear()
        self.widgetItems = []
        self.progressBarInit()
##        milestones = set([i for i in range(0, len(geneinfo), max(len(geneinfo)/100, 1))])
        milestones = set([i for i in range(0, len(geneinfo), max(len(geneinfo)/100, 1))])
        for i, (gene, gi) in enumerate(geneinfo):
            if gi:
                item = QTreeWidgetItem(self.treeWidget,
                                       ["", gi.symbol, gi.locus_tag or "",
                                        gi.chromosome or "", gi.description or "", ", ".join(gi.synonyms),
                                        gi.symbol_from_nomenclature_authority or ""])
                item.info = gi
                link = LinkItem(gi.gene_id, self.treeWidget)
                self.treeWidget.setItemWidget(item, 0, link)
##                link.show()
                self.widgetItems.append(item)
            if i in milestones:
                self.progressBarSet(100.0*i/len(geneinfo))
##        self.treeWidget.addTopLevelItems(self.widgetItems)
        self.progressBarFinished()
##        self.widgetItems[-1].setText(0, "")
##        self.treeWidget.update(self.treeWidget.indexFromItem(self.widgetItems[-1]))
        self.treeWidget.viewport().update()
        self.infoLabel.setText("%i genes\n%i matched NCBI's IDs" % (len(genes), len(self.widgetItems)))

    def clear(self):
        self.infoLabel.setText("No data on input\n")
        self.treeWidget.clear()
        self.send("Selected Examples", None)

    def commitIf(self):
        if self.autoCommit:
            self.commit()

    def commit(self):
        selected = [self.treeWidget.itemFromIndex(index).info for index in self.treeWidget.selectedIndexes()]
##        selected = [item.info for item in self.widgetItems if not item.isHidden() and item.isSelected()]
        if self.useAttr:
            attrs = [attr for attr, (name, gi) in zip(self.data.domain.attributes, self.geneinfo) if gi in selected]
            domain = orange.Domain(attrs, self.data.domain.classVar)
            domain.addmetas(self.data.domain.getmetas())
            self.send("Selected Examples", orange.ExampleTable(domain, self.data))
        else:
            attr = self.attributes[self.geneAttr]
            geneinfo = dict(self.geneinfo)
            examples = [ex for ex in self.data if geneinfo.get(str(ex[attr])) in selected]
            self.send("Selected Examples", orange.ExampleTable(examples) if examples else None)
        
    def searchUpdate(self):
        searchStrings = self.searchString.lower().split()
        for item in self.widgetItems:
            item.setHidden(not all(any(string in str(item.text(i)).lower() for i in range(7)) for string in searchStrings))

    def selectFiltered(self):
        for item in self.widgetItems:
            item.setSelected(not item.isHidden())

    def clearSelection(self):
        for item in self.widgetItems:
            item.setSelected(False)
                

if __name__ == "__main__":
    app = QApplication(sys.argv)
    data = orange.ExampleTable("../../orange/doc/datasets/brown-selected.tab")
    w = OWGeneInfo()
    w.show()
    w.setData(data)
    app.exec_()
    w.saveSettings()
        
        
        
        