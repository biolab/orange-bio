"""
<name>Gene Info</name>
<description>Displays gene information from NCBI.</description>
<priority>2010</priority>
<contact>Ales Erjavec (ales.erjevec(@at@)fri.uni-lj.si)</contact>
<icon>icons/GeneInfo.png</icon>
"""
from __future__ import with_statement

import obiGene, obiTaxonomy
import orange
import orngServerFiles

from OWWidget import *
import OWGUI

from orngDataCaching import data_hints

from collections import defaultdict
from functools import partial

class TreeModel(QAbstractItemModel):
    def __init__(self, data, header, parent):
        QAbstractItemModel.__init__(self, parent)
        self._data = [[QVariant(s) for s in row] for row in data]
        self._dataDict = {}
        self._header = header
        self._roleData = {Qt.DisplayRole:self._data}
        self._roleData = partial(defaultdict, partial(defaultdict, partial(defaultdict, QVariant)))(self._roleData)
    
    def setColumnLinks(self, column, links):
        font =QFont()
        font.setUnderline(True)
        font = QVariant(font)
        for i, link in enumerate(links):
            self._roleData[LinkRole][i][column] = QVariant(link)
            self._roleData[Qt.FontRole][i][column] = font
            self._roleData[Qt.ForegroundRole][i][column] = QVariant(QColor(Qt.blue))
    
    def setRoleData(self, role, row, col, data):
        self._roleData[role][row][col] = data
        
    def data(self, index, role):
        row, col = index.row(), index.column()
        return self._roleData[role][row][col]
        
    def index(self, row, col, parent=QModelIndex()):
        return self.createIndex(row, col, 0)
    
    def parent(self, index):
        return QModelIndex()
    
    def rowCount(self, index):
        if index.isValid():
            return 0
        else:
            return len(self._data)
        
    def columnCount(self, index):
        return len(self._header)

    def headerData(self, section, orientation, role):
        if role==Qt.DisplayRole:
            return QVariant(self._header[section])
        return QVariant()

from OWGUI import LinkStyledItemDelegate, LinkRole
        
class OWGeneInfo(OWWidget):
    settingsList = ["organismIndex", "geneAttr", "useAttr", "autoCommit"]
    contextHandlers = {"":DomainContextHandler("", ["organismIndex", "geneAttr", "useAttr"])}
    def __init__(self, parent=None, signalManager=None, name="Gene Info"):
        OWWidget.__init__(self, parent, signalManager, name)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Selected Examples", ExampleTable)]

        self.organismIndex = 0
        self.geneAttr = 0
        self.useAttr = False
        self.autoCommit = False
        self.searchString = ""
        self.selectionChangedFlag = False
        self.loadSettings()
        
        self.infoLabel = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Info", addSpace=True), "No data on input\n")
        self.organisms = sorted(set([name.split(".")[-2] for name in orngServerFiles.listfiles("NCBI_geneinfo")] + obiGene.NCBIGeneInfo.essential_taxids()))
    
        self.orgaismsComboBox = OWGUI.comboBox(self.controlArea, self, "organismIndex", "Organism", items=[obiTaxonomy.name(id) for id in self.organisms], callback=self.setItems, debuggingEnabled=0)
        OWGUI.separator(self.controlArea)
        box = OWGUI.widgetBox(self.controlArea, "Gene names", addSpace=True)
        self.geneAttrComboBox = OWGUI.comboBox(box, self, "geneAttr", "Gene atttibute", callback=self.setItems)
        c = OWGUI.checkBox(box, self, "useAttr", "Use attribute names", callback=self.setItems, disables=[(-1, self.geneAttrComboBox)])
        self.geneAttrComboBox.setDisabled(bool(self.useAttr))

        box = OWGUI.widgetBox(self.controlArea, "Commit", addSpace=True)
        b = OWGUI.button(box, self, "Commit", callback=self.commit)
        c = OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.setStopper(self, b, c, "selectionChangedFlag", callback=self.commit)
        
        ##A label for dictyExpress link
        self.dictyExpressBox = OWGUI.widgetBox(self.controlArea, "Dicty Express")
        self.linkLabel = OWGUI.widgetLabel(self.dictyExpressBox, "")
        self.linkLabel.setOpenExternalLinks(True)
        self.dictyExpressBox.hide()
        
        OWGUI.rubber(self.controlArea)

        OWGUI.lineEdit(self.mainArea, self, "searchString", "Filter", callbackOnType=True, callback=self.searchUpdate)
#        self.treeWidget = QTreeWidget(self.mainArea)
        self.treeWidget = QTreeView(self.mainArea)
        #self.treeWidget.setHeaderLabels(["NCBI ID", "Symbol", "Locus Tag", "Chromosome", "Description", "Synonyms", "Nomenclature"])
        self.treeWidget.setRootIsDecorated(False)
        self.treeWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.treeWidget.setItemDelegate(LinkStyledItemDelegate(self.treeWidget))
        #self.connect(self.treeWidget, SIGNAL("itemSelectionChanged()"), self.commitIf)
        self.treeWidget.viewport().setMouseTracking(True)
        self.treeWidget.setSortingEnabled(True)
        self.mainArea.layout().addWidget(self.treeWidget)
        
        box = OWGUI.widgetBox(self.mainArea, "", orientation="horizontal")
        OWGUI.button(box, self, "Select Filtered", callback=self.selectFiltered)
        OWGUI.button(box, self, "Clear Selection", callback=self.treeWidget.clearSelection)
        
        self.resize(1000, 700)        

        self.geneinfo = []
        self.cells = []
        self.data = None
        self.currentLoaded = None, None
        self.selectionUpdateInProgress = False
        
    def setData(self, data=None):
        self.closeContext()
        self.data = data
        if data:
            self.geneAttrComboBox.clear()
            self.attributes = [attr for attr in self.data.domain.variables + self.data.domain.getmetas().values() if attr.varType in [orange.VarTypes.String, orange.VarTypes.Discrete]]
            self.geneAttrComboBox.addItems([attr.name for attr in self.attributes])
            self.openContext("", data)
            self.geneAttr = min(self.geneAttr, len(self.attributes) - 1)
            
            taxid = data_hints.get_hint(self.data, "taxid", "")
            if taxid in self.organisms:
                self.organismIndex = self.organisms.index(taxid)
                
            self.useAttr = data_hints.get_hint(self.data, "genesinrows",  self.useAttr)
            
            self.setItems()
        else:
            self.clear()

    def setItems(self):
        self.warning(0)
        if not self.data:
            return
        if self.useAttr:
            genes = [attr.name for attr in self.data.domain.attributes]
        elif self.attributes:
            attr = self.attributes[self.geneAttr]
            genes = [str(ex[attr]) for ex in self.data if not ex[attr].isSpecial()]
        else:
            genes = []
        if not genes:
            self.warning(0, "Could not extract genes from input dataset.")
        self.warning(1)
        org = self.organisms[min(self.organismIndex, len(self.organisms) - 1)]
        info , currorg = self.currentLoaded
        self.error(0)
        if currorg != org:
            self.progressBarInit()
            
            ## Load the gene info in a worker thread
            call = self.asyncCall(obiGene.NCBIGeneInfo, (org,), name="Loading NCBI Gene info", blocking=False)
            call.connect(call, SIGNAL("progressChanged(float)"), self.progressBarSet, Qt.QueuedConnection)
            with orngServerFiles.DownloadProgress.setredirect(call.emitProgressChanged):
                call.__call__()
                info = call.get_result()

#            with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
#                info = obiGene.NCBIGeneInfo(self.organisms[min(self.organismIndex, len(self.organisms) - 1)])
            self.progressBarFinished()
            self.currentLoaded = info, org
            
        self.updateDictyExpressLink(genes, show=org == "352472")
            
        self.geneinfo = geneinfo = [(gene, info.get_info(gene, None)) for gene in genes]

        self.progressBarInit()
        milestones = set([i for i in range(0, len(geneinfo), max(len(geneinfo)/100, 1))])
        self.cells = cells = []
        links = []
        for i, (gene, gi) in enumerate(geneinfo):
            if gi:
                cells.append([gi.gene_id, gi.symbol + " (%s)" % gene if gene != gi.symbol else gi.symbol,
                            gi.locus_tag or "", gi.chromosome or "", gi.description or "",
                            ", ".join(gi.synonyms), gi.symbol_from_nomenclature_authority or ""])
                links.append("http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=%s" % gi.gene_id)

            if i in milestones:
                self.progressBarSet(100.0*i/len(geneinfo))
        model = TreeModel(cells, ["NCBI ID", "Symbol", "Locus Tag", "Chromosome", "Description", "Synonyms", "Nomenclature"], self.treeWidget)
        model.setColumnLinks(0, links)
        proxyModel = QSortFilterProxyModel(self)
        proxyModel.setSourceModel(model)
        self.treeWidget.setModel(proxyModel)
        self.connect(self.treeWidget.selectionModel(), SIGNAL("selectionChanged(QItemSelection , QItemSelection )"), self.commitIf)
        for i in range(7):
            self.treeWidget.resizeColumnToContents(i)
            self.treeWidget.setColumnWidth(i, min(self.treeWidget.columnWidth(i), 200))
        self.treeWidget.update()
        self.progressBarFinished()

        self.infoLabel.setText("%i genes\n%i matched NCBI's IDs" % (len(genes), len(cells)))
        self.matchedInfo = len(genes), len(cells)

    def clear(self):
        self.infoLabel.setText("No data on input\n")
        self.treeWidget.setModel(TreeModel([], ["NCBI ID", "Symbol", "Locus Tag", "Chromosome", "Description", "Synonyms", "Nomenclature"], self.treeWidget))
        self.geneAttrComboBox.clear()
        self.send("Selected Examples", None)

    def commitIf(self, *args):
        if self.autoCommit and not self.selectionUpdateInProgress:
            self.commit()
        else:
            self.selectionChangedFlag = True

    def commit(self):
        if not self.data:
            return
        
        mapToSource = self.treeWidget.model().mapToSource
        selectedIds = [self.cells[mapToSource(index).row()][0] for index in self.treeWidget.selectedIndexes()]
        
        selected = [gi for gene, gi in self.geneinfo if gi and gi.gene_id in selectedIds]
        if self.useAttr:
            attrs = [attr for attr, (name, gi) in zip(self.data.domain.attributes, self.geneinfo) if gi in selected]
            domain = orange.Domain(attrs, self.data.domain.classVar)
            domain.addmetas(self.data.domain.getmetas())
            newdata = orange.ExampleTable(domain, self.data)
            self.send("Selected Examples", newdata)
        else:
            attr = self.attributes[self.geneAttr]
            geneinfo = dict(self.geneinfo)
            examples = [ex for ex in self.data if geneinfo.get(str(ex[attr])) in selected]
            if examples:
                newdata = orange.ExampleTable(examples)
            else:
                newdata = None
            self.send("Selected Examples", newdata)
            
    def rowFiltered(self, row):
        searchStrings = self.searchString.lower().split()
        row = unicode(" ".join(self.cells[row]).lower(), errors="ignore")
        return not all([s in row for s in searchStrings])
    
    def searchUpdate(self):
        if not self.data:
            return
        searchStrings = self.searchString.lower().split()
        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            row = unicode(" ".join(row).lower(), errors="ignore")
            self.treeWidget.setRowHidden(mapFromSource(index(i, 0)).row(), QModelIndex(), not all([s in row for s in searchStrings]))
        #self.treeWidget.model().setFilterRegExp(QRegExp(self.searchString, Qt.CaseInsensitive, QRegExp.FixedString))
            
    def selectFiltered(self):
        if not self.data:
            return
        itemSelection = QItemSelection()
        
        index = self.treeWidget.model().sourceModel().index
        mapFromSource = self.treeWidget.model().mapFromSource
        for i, row in enumerate(self.cells):
            if not self.rowFiltered(i):
                itemSelection.select(mapFromSource(index(i, 0)), mapFromSource(index(i, 0)))
        self.treeWidget.selectionModel().select(itemSelection, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        
    def sendReport(self):
        import OWReport
        genes, matched = self.matchedInfo
        info, org = self.currentLoaded
        self.reportRaw("<p>Input: %i genes of which %i (%.1f%%) matched NCBI synonyms<br>Organism: %s<br>Filter: %s</p>" % (genes, matched, 100.0 * matched / genes, obiTaxonomy.name(org), self.searchString))
        self.reportSubsection("Gene list")
        self.reportRaw(reportItemView(self.treeWidget))
        
    def updateDictyExpressLink(self, genes, show=False):
        def fix(ddb):
            if ddb.startswith("DDB"): 
                if not ddb.startswith("DDB_G"):
                    ddb = ddb.replace("DDB", "DDB_G")
                return ddb
            return None 
        if show:
            genes = [fix(gene) for gene in genes if fix(gene)]
            link1 = '<a href="http://www.ailab.si/dictyexpress/run/index.php?gene=%s">Microarray profile</a>' % (" ".join(genes))
            link2 = '<a href="http://www.ailab.si/dictyexpress/run/index.php?gene=%s&db=rnaseq">RNA-Seq profile</a>' % (" ".join(genes))
            self.linkLabel.setText(link1 + "<br/>" + link2)
            
            show = any(genes)
                
        if show:
            self.dictyExpressBox.show()
        else:
            self.dictyExpressBox.hide()
        
def reportItemView(view):
    model = view.model()
    return reportItemModel(view, model)
    
def reportItemModel(view, model, index=QModelIndex()):
    if not index.isValid() or model.hasChildren(index):
        columnCount, rowCount = model.columnCount(index), model.rowCount(index)
        if not index.isValid():
            text = '<table>\n<tr>' + ''.join('<th>%s</th>' % model.headerData(i, Qt.Horizontal, Qt.DisplayRole).toString() for i in range(columnCount)) +'</tr>\n'
        else:
#            variant = model.data(index, Qt.DisplayRole)
#            text = '<table' + (' caption="%s"' % variant.toString() if variant.isValid() else '') + '>\n'
            pass
        text += ''.join('<tr>' + ''.join('<td>' + reportItemModel(view, model, model.index(row, column, index)) + '</td>' for column in range(columnCount)) + '</tr>\n' for row in range(rowCount) if not view.isRowHidden(row, index))
        text += '</table>'
        return text
    else:
        variant = model.data(index, Qt.DisplayRole)
        return str(variant.toString()) if variant.isValid() else ""
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    data = orange.ExampleTable("../../orange/doc/datasets/brown-selected.tab")
    w = OWGeneInfo()
    w.show()
    w.setData(data)
    app.exec_()
    w.saveSettings()
        
        
        
        
