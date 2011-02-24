"""
<name>KEGG Pathway Browser</name>
<description>Browse KEGG pathways that include an input set of genes.</description>
<priority>2030</priority>
<icon>icons/KEGG.png</icon>
"""
from __future__ import with_statement 

import sys
import orange
import obiKEGG, orngServerFiles
import obiTaxonomy

import webbrowser

from OWWidget import *
import OWGUI
from collections import defaultdict

from orngDataCaching import data_hints

def split_and_strip(string, sep=None):
    return [s.strip() for s in string.split(sep)]

class EntryGraphicsItem(QGraphicsPathItem):
    def __init__(self, graphics, *args):
        QGraphicsPathItem.__init__(self, *args)
        path = QPainterPath()
        x, y, w, h = [int(graphics.get(c, 0)) for c in ["x", "y", "width", "height"]]
        type = graphics.get("type", "rectangle")
        if type == "rectangle":
            path.addRect(QRectF(x - w/2, y - h/2, w, h))
        elif type == "roundrectangle":
            path.addRoundedRect(QRectF(x - w/2, y - h/2, w, h), 10, 10)
        elif type == "circle":
            path.addEllipse(QRectF(x - w/2, y - h/2, w, h))
            
        self.setPath(path)
        self.setAcceptHoverEvents(True)
        self._actions = []
        self.link = None
        
    def hoverEnterEvent(self, event):
        self.setBrush(QBrush(QColor(0, 100, 0, 100)))
        
    def hoverLeaveEvent(self, event):
        self.setBrush(QBrush(Qt.NoBrush))
        
    def contextMenuEvent(self, event):
        if self._actions:
            self._menu = menu = QMenu()
            for action in self._actions:
                menu.addAction(action)
            menu.popup(event.screenPos())
            
    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemSelectedHasChanged:
            self.setPen(QPen(Qt.red if self.isSelected() else Qt.blue, 2))
            
        return QGraphicsPathItem.itemChange(self, change, value)
            
class GraphicsPathwayItem(QGraphicsPixmapItem):
    def __init__(self, pathway, objects, *args, **kwargs):
        QGraphicsPixmapItem.__init__(self, *args)
        self.setTransformationMode(Qt.SmoothTransformation)
        self.setPathway(pathway)
        self.setMarkedObjects(objects, name_mapper=kwargs.get("name_mapper", {}))
        
    def setPathway(self, pathway):
        self.pathway = pathway
        if pathway:
            image_filename = pathway.get_image()
            self._pixmap = QPixmap(image_filename)
        else:
            self._pixmap = QPixmap()
        self.setPixmap(self._pixmap)
        
    def setMarkedObjects(self, objects, name_mapper={}):
        for entry in self.pathway.entrys() if self.pathway else []:
            if entry.type == "group":
                continue
            graphics = entry.graphics
            contained_objects = [obj for obj in objects if obj in entry.name]
            item = EntryGraphicsItem(graphics, self, self.scene())
            item.setToolTip(self.tooltip(entry, contained_objects, name_mapper))
            item._actions = self.actions(entry, contained_objects)
            item.marked_objects = contained_objects
            if contained_objects:
                item.setPen(QPen(Qt.blue, 2))
                item.setFlag(QGraphicsItem.ItemIsSelectable, True)
    
    def actions(self, entry, marked_objects=[]):
        actions = []
        type = entry.type
        if marked_objects:
            action = QAction("View genes on kegg website", None)
            org = set([s.split(":")[0] for s in marked_objects]).pop()
            genes = [s.split(":")[-1] for s in marked_objects]
            address = "http://www.genome.jp/dbget-bin/www_bget?" + "+".join([org] + genes)
            QObject.connect(action, SIGNAL("triggered()"), lambda toggled=False, address=address: webbrowser.open(address))
            actions.append(action)
        elif hasattr(entry, "link"):
            action = QAction("View %s on KEGG website" % str(type), None)
            QObject.connect(action, SIGNAL("triggered()"), lambda toggled=False, address=entry.link: webbrowser.open(address))
            actions.append(action)
        return actions
    
    def tooltip(self, entry, objects, name_mapper={}):
        names = [obj for obj in objects if obj in entry.name]
        names = [name_mapper.get(name, name) for name in names]
        text = "<p>%s</p>" % (entry.name[:16] + " ..." if len(entry.name) > 20 else entry.name)
        if names:
            text += "<br>".join(names)
        return text
    
    def contextMenuEvent(self, event):
        self._menu = menu = QMenu()
        action = menu.addAction("View this pathway on KEGG website")
        address ="http://www.kegg.jp/kegg-bin/show_pathway?%s%s" % (self.pathway.org, self.pathway.number)
        QObject.connect(action, SIGNAL("triggered()"), lambda : webbrowser.open(address))
        menu.popup(event.screenPos())

class PathwayView(QGraphicsView):
    def __init__(self, master, *args):
        QGraphicsView.__init__(self, *args)
        self.master = master
        
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        
        self.setRenderHints(QPainter.Antialiasing)
        scene = QGraphicsScene(self)
        self.pixmapGraphicsItem = QGraphicsPixmapItem(None, scene) 
        self.setScene(scene)
        
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        
        self.setFocusPolicy(Qt.WheelFocus)
        
    def SetPathway(self, pathway=None, objects=[]):
        self.scene().clear()
        self.pathway = pathway
        self.objects = objects
        self.pathwayItem = GraphicsPathwayItem(pathway, objects, None, self.scene(), name_mapper=getattr(self.master, "uniqueGenesDict", {}))
        self.scene().setSceneRect(self.pathwayItem.boundingRect())
        self.updateTransform()
    
    def resizeEvent(self, event):
        self.updateTransform()
        return QGraphicsView.resizeEvent(self, event)
            
    def updateTransform(self):
        if self.master.autoResize:
            self.fitInView(self.scene().sceneRect().adjusted(-1, -1, 1, 1), Qt.KeepAspectRatio)
        else:
            self.setTransform(QTransform())
            
    def paintEvent(self, event):
        QGraphicsView.paintEvent(self, event)
        if getattr(self, "_userMessage", None):
            painter = QPainter(self.viewport())
            font = QFont(self.font())
            font.setPointSize(15)
            painter.setFont(font)
            painter.drawText(self.viewport().geometry(), Qt.AlignCenter, self._userMessage)
            painter.end()
        
    
class OWKEGGPathwayBrowser(OWWidget):
    settingsList = ["organismIndex", "geneAttrIndex", "autoCommit", "autoResize", "useReference", "useAttrNames", "caseSensitive", "showOrthology"]
    contextHandlers = {"":DomainContextHandler("",[ContextField("organismIndex", DomainContextHandler.Required + DomainContextHandler.IncludeMetaAttributes),
                                                   ContextField("geneAttrIndex", DomainContextHandler.Required + DomainContextHandler.IncludeMetaAttributes),
                                                   ContextField("useAttrNames", DomainContextHandler.Required + DomainContextHandler.IncludeMetaAttributes)])}
    def __init__(self, parent=None, signalManager=None, name="KEGG Pathway Browser"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Examples", ExampleTable, self.SetData), ("Reference", ExampleTable, self.SetRefData)]
        self.outputs = [("Selected Examples", ExampleTable), ("Unselected Examples", ExampleTable)]
        self.organismIndex = 0
        self.geneAttrIndex = 0
        self.autoCommit = False
        self.autoResize = True
        self.useReference = False
        self.useAttrNames = 0
        self.caseSensitive = True
        self.showOrthology = True
        self.autoFindBestOrg = False
        self.loadSettings()

        self.controlArea.setMaximumWidth(250)
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.infoLabel = OWGUI.widgetLabel(box, "No data on input\n")
        
        self.allOrganismCodes = {} 

        self.organismCodes = []

        self.organismComboBox = cb = OWGUI.comboBox(self.controlArea, self, "organismIndex", box="Organism", items=[], callback=self.Update, addSpace=True, debuggingEnabled=0)
        cb.setMaximumWidth(200)
        
        self.signalManager.freeze(self).push() #setFreeze(1)
        QTimer.singleShot(100, self.UpdateOrganismComboBox)
        
        box = OWGUI.widgetBox(self.controlArea, "Gene attribute")
        self.geneAttrCombo = OWGUI.comboBox(box, self, "geneAttrIndex", callback=self.Update)
        OWGUI.checkBox(box, self, "useAttrNames", "Use variable names", disables=[(-1, self.geneAttrCombo)], callback=self.UseAttrNamesCallback)
        self.geneAttrCombo.setDisabled(bool(self.useAttrNames))
        
#        OWGUI.checkBox(box, self, "caseSensitive", "Case sensitive gene matching", callback=self.Update)
        OWGUI.separator(self.controlArea)
        
        OWGUI.checkBox(self.controlArea, self, "useReference", "From signal", box="Reference", callback=self.Update)
        OWGUI.separator(self.controlArea)

        OWGUI.checkBox(self.controlArea, self, "showOrthology", "Show pathways in full orthology", box="Orthology", callback=self.UpdateListView)
        
        OWGUI.checkBox(self.controlArea, self, "autoResize", "Resize to fit", box="Image", callback=lambda :self.pathwayView.updateTransform())
        OWGUI.separator(self.controlArea)

        box = OWGUI.widgetBox(self.controlArea, "Selection")
        OWGUI.checkBox(box, self, "autoCommit", "Commit on update")
        OWGUI.button(box, self, "Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)
        
        spliter = QSplitter(Qt.Vertical, self.mainArea)
        self.pathwayView = PathwayView(self, spliter)
        self.mainArea.layout().addWidget(spliter)

        self.listView = QTreeWidget(spliter)
        spliter.addWidget(self.listView)
        
        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(4)
        self.listView.setHeaderLabels(["Pathway", "P value", "Genes", "Reference"])

        self.listView.setSelectionMode(QAbstractItemView.SingleSelection)
            
        self.listView.setSortingEnabled(True)
        #self.listView.setAllColumnsShowFocus(1)
        self.listView.setMaximumHeight(200)
        
        self.connect(self.listView, SIGNAL("itemSelectionChanged()"), self.UpdatePathwayView)
        
        self.ctrlPressed = False
        self.selectedObjects = defaultdict(list)
        self.data = None
        self.refData = None
        self.loadedOrganism = None
        
        self.resize(800, 600)
        
        self.connect(self, SIGNAL("widgetStateChanged(QString, int, QString)"), self.onStateChange)
        
    def UpdateOrganismComboBox(self):
        try:
            self.progressBarInit()
            with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
                genome = obiKEGG.KEGGGenome()
            self.progressBarFinished()
            
            self.allOrganismCodes = genome 
    
            essential = genome.essential_organisms()
            
            local = [name.split(".")[0].split("_")[-1] for name in orngServerFiles.listfiles("KEGG") if "kegg_genes" in name]
            self.organismCodes = [(code, organism.definition) for code, organism in self.allOrganismCodes.items() if code in local or code in essential]
            self.organismCodes.sort()
            items = [desc for code, desc in self.organismCodes]
            self.organismCodes = [code for code, desc in self.organismCodes]
            
            self.organismComboBox.addItems(items)
        finally:
            self.signalManager.freeze(self).pop() #setFreeze(0)

        
    def SetData(self, data=None):
        if not self.organismCodes: ## delay this call until we retrieve organism codes from the server files 
            QTimer.singleShot(200, lambda: self.SetData(data))
            return

        self.closeContext()
        self.data = data
        self.warning(0)
        if data:
            self.SetGeneAttrCombo()
            taxid = data_hints.get_hint(data, "taxid", None)
            if taxid:
                try:
                    code = obiKEGG.from_taxid(taxid)
                    self.organismIndex = self.organismCodes.index(code)
                except Exception, ex:
                    print ex, taxid
            
            self.useAttrNames = data_hints.get_hint(data, "genesinrows", self.useAttrNames)
            
            self.openContext("", data)
            self.Update()
        else:
            self.infoLabel.setText("No data on input\n")
            self.listView.clear()
            self.selectedObjects = defaultdict(list)
            self.pathwayView.SetPathway(None)
            self.send("Selected Examples", None)
            self.send("Unselected Examples", None)

    def SetRefData(self, data=None):
        self.refData = data
        if self.useReference and self.data and self.organismCodes:
            self.Update()

    def UseAttrNamesCallback(self):
##        self.geneAttrCombo.setDisabled(bool(self.useAttrNames))
        self.Update()

    def SetGeneAttrCombo(self):
        self.geneAttrCandidates = self.data.domain.attributes + self.data.domain.getmetas().values()
        self.geneAttrCandidates = filter(lambda v:v.varType in [orange.VarTypes.Discrete ,orange.VarTypes.String], self.geneAttrCandidates)
        self.geneAttrCombo.clear()
        #print 'geneAttrCandidates', self.geneAttrCandidates
        self.geneAttrCombo.addItems([var.name for var in self.geneAttrCandidates])
        return
                
    def PreDownload(self, org=None, pb=None):
        pb, finish = (OWGUI.ProgressBar(self, 0), True) if pb is None else (pb, False)
        files = ["kegg_brite.tar.gz", "kegg_pathways_map.tar.gz", "kegg_genome.tar.gz"]
        if org:
            files += ["kegg_genes_%s.tar.gz" % org, "kegg_pathways_%s.tar.gz" % org]
        files = [file for file in files if file not in orngServerFiles.listfiles("KEGG")]
        pb.iter += len(files) * 100
        for i, filename in enumerate(files):
#            print filename
            orngServerFiles.download("KEGG", filename, callback=pb.advance)
        if finish:
            pb.finish()
            
    def UpdateListView(self):
        self.bestPValueItem = None
        self.listView.clear()
        if not self.data:
            return
        allPathways = self.org.pathways()
        allRefPathways = obiKEGG.pathways("map")
        self.progressBarFinished()
        items = []
        if self.showOrthology:
            self.koOrthology = obiKEGG.KEGGBrite("ko00001")
            self.listView.setRootIsDecorated(True)
            path_ids = set([s[-5:] for s in self.pathways.keys()])
            def _walkCollect(koEntry):
                num = koEntry.title[:5] if koEntry.title else None
                if num  in path_ids:
                    return [koEntry] + reduce(lambda li,c:li+_walkCollect(c), [child for child in koEntry.entrys], [])
                else:
                    c = reduce(lambda li,c:li+_walkCollect(c), [child for child in koEntry.entrys], [])
                    return c + (c and [koEntry] or [])
            allClasses = reduce(lambda li1, li2: li1+li2, [_walkCollect(c) for c in self.koOrthology], [])
            def _walkCreate(koEntry, lvItem):
                item = QTreeWidgetItem(lvItem)
                id = "path:"+self.organismCodes[min(self.organismIndex, len(self.organismCodes)-1)] + koEntry.title[:5]
                if koEntry.title[:5] in path_ids:
                    genes, p_value, ref = self.pathways[id]
                    item.setText(0, obiKEGG.KEGGPathway(id).title)
#                    print id, obiKEGG.KEGGPathway(id).title
                    item.setText(1, "%.5f" % p_value)
                    item.setText(2, "%i of %i" %(len(genes), len(self.genes)))
                    item.setText(3, "%i of %i" %(ref, len(self.referenceGenes)))
                    item.pathway_id = id
                else:
                    item.setText(0, obiKEGG.KEGGPathway(id).title if id in allPathways else koEntry.title)
                    if id in allPathways:
                        item.pathway_id = id
                    elif "path:map" + koEntry.title[:5] in allRefPathways:
                        item.pathway_id = "path:map" + koEntry.title[:5]
                    else:
                        item.pathway_id = None
                
                for child in koEntry.entrys:
                    if child in allClasses:
                        _walkCreate(child, item)
            
            for koEntry in self.koOrthology:
                if koEntry in allClasses:
                    _walkCreate(koEntry, self.listView)
                    
            self.listView.update()
        else:
            self.listView.setRootIsDecorated(False)
            pathways = self.pathways.items()
            pathways.sort(lambda a,b:cmp(a[1][1], b[1][1]))
            for id, (genes, p_value, ref) in pathways:
                item = QTreeWidgetItem(self.listView)
                item.setText(0, obiKEGG.KEGGPathway(id).title)
                item.setText(1, "%.5f" % p_value)
                item.setText(2, "%i of %i" %(len(genes), len(self.genes)))
                item.setText(3, "%i of %i" %(ref, len(self.referenceGenes)))
                item.pathway_id = id
                items.append(item)
                
        self.bestPValueItem = items and items[0] or None
        self.listView.expandAll()
        for i in range(4):
            self.listView.resizeColumnToContents(i)
            
        if self.bestPValueItem:
            self.listView.selectionModel().select(self.listView.indexFromItem(self.bestPValueItem), QItemSelectionModel.ClearAndSelect)

    def UpdatePathwayView(self):
        items = self.listView.selectedItems()
        
        if len(items) > 0:
            item = items[0]
        else:
            item = None
            
        self.selectedObjects = defaultdict(list)
        self.Commit()
        item = item or self.bestPValueItem
        if not item or not item.pathway_id:
            self.pathwayView.SetPathway(None)
            return
        self.pathway = obiKEGG.KEGGPathway(item.pathway_id)
        self.pathwayView.SetPathway(self.pathway, self.pathways.get(item.pathway_id, [[]])[0])
            
            
    def Update(self):
        if not self.data:
            return
        self.error(0)
        self.information(0)
        pb = OWGUI.ProgressBar(self, 100)
        if self.useAttrNames:
            genes = [str(v.name).strip() for v in self.data.domain.attributes]
        elif self.geneAttrCandidates:
            geneAttr = self.geneAttrCandidates[min(self.geneAttrIndex, len(self.geneAttrCandidates)-1)]
            genes = [str(e[geneAttr]) for e in self.data if not e[geneAttr].isSpecial()]
            if any("," in gene for gene in genes):
                genes = reduce(list.__add__, (split_and_strip(gene, ",") for gene in genes), [])
                self.information(0, "Separators detected in input gene names. Assuming multiple genes per example.")
        else:
            self.error(0, "Cannot extact gene names from input")
            genes = []
        org_code = self.organismCodes[min(self.organismIndex, len(self.organismCodes)-1)]
        if self.loadedOrganism != org_code:
            self.PreDownload(org_code, pb=pb)
            self.org = obiKEGG.KEGGOrganism(org_code)
            self.loadedOrganism = org_code
        uniqueGenes, conflicting, unknown = self.org.get_unique_gene_ids(set(genes), self.caseSensitive)
        genesCount = len(set(genes))
        self.infoLabel.setText("%i unique gene names on input\n%i (%.1f%%) genes names matched" % (genesCount, len(uniqueGenes), 100.0*len(uniqueGenes)/genesCount if genes else 0.0))  
#        if conflicting:
#            print >> sys.stderr, "Conflicting genes:", conflicting
#        if unknown:
#            print >> sys.stderr, "Unknown genes:", unknown
        self.information(1)
        if self.useReference and self.refData:
            if self.useAttrNames:
                reference = [str(v.name).strip() for v in self.refData]
            else:
                geneAttr = self.geneAttrCandidates[min(self.geneAttrIndex, len(self.geneAttrCandidates)-1)]
                reference = [str(e[geneAttr]) for e in self.refData if not e[geneAttr].isSpecial()]
                if any("," in gene for gene in reference):
                    reference = reduce(list.__add__, (split_and_strip(gene, ",") for gene in reference), [])
                    self.information(1, "Separators detected in reference gene names. Assuming multiple genes per example.")
            uniqueRefGenes, conflicting, unknown = self.org.get_unique_gene_ids(set(reference), self.caseSensitive)
            self.referenceGenes = reference = uniqueRefGenes.keys()
        else:
            self.referenceGenes = reference = self.org.get_genes()
        self.uniqueGenesDict = uniqueGenes
        self.genes = uniqueGenes.keys()
        self.revUniqueGenesDict = dict([(val, key) for key, val in self.uniqueGenesDict.items()])
#        self.progressBarInit()
#        with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
        self.pathways = self.org.get_enriched_pathways(self.genes, reference, callback=lambda value: pb.advance()) #self.progressBarSet)
        if not self.pathways:
            self.warning(0, "No enriched pathways found.")
        else:
            self.warning(0)
            
#        self.progressBarFinished()
        self.UpdateListView()
        pb.finish()
##        print self.bestPValueItem
        #self.bestPValueItem.setSelected(True)
        #self.UpdatePathwayView()

    def SelectObjects(self, objs):
        if (not self.selectedObjects or self.ctrlPressed) and not objs:
            return
        if self.ctrlPressed:
            for id, graphics in objs:
                graphics = tuple(sorted(graphics.items()))
                if id in self.selectedObjects[graphics]:
                    self.selectedObjects[graphics].pop(self.selectedObjects[graphics].index(id))
                    if not self.selectedObjects[graphics]:
                        del self.selectedObjects[graphics]
                else:
                    self.selectedObjects[graphics].append(id)
        else:
            self.selectedObjects.clear()
            for id, graphics in objs:
                graphics = tuple(sorted(graphics.items()))
                self.selectedObjects[graphics].append(id)
        if self.autoCommit:
            self.Commit()
            

    def Commit(self):
        if self.data:
            selectedItems = self.pathwayView.scene().selectedItems()
            selectedGenes = reduce(set.union, [item.marked_objects for item in selectedItems], set())
            
            if self.useAttrNames:
#                selectedGenes = reduce(set.union, self.selectedObjects.values(), set())
                selectedVars = [self.data.domain[self.uniqueGenesDict[gene]] for gene in selectedGenes]
                newDomain = orange.Domain(selectedVars ,0)
                data = orange.ExampleTable(newDomain, self.data)
                self.send("Selected Examples", data)
            elif self.geneAttrCandidates:
                geneAttr = self.geneAttrCandidates[min(self.geneAttrIndex, len(self.geneAttrCandidates)-1)]
                selectedExamples = []
                otherExamples = []
#                selectedGenes = reduce(set.union, self.selectedObjects.values(), set())
                for ex in self.data:
                    names = [self.revUniqueGenesDict.get(name, None) for name in split_and_strip(str(ex[geneAttr]), ",")]
                    if any(name and name in selectedGenes for name in names):
                        selectedExamples.append(ex)
                    else:
                        otherExamples.append(ex)
                        
                if selectedExamples:
                    selectedExamples = orange.ExampleTable(selectedExamples)
                else:
                    selectedExamples = None
                    
                if otherExamples:
                    otherExamples = orange.ExampleTable(otherExamples)
                else:
                    otherExamples = None
                    
                self.send("Selected Examples", selectedExamples)
                self.send("Unselected Examples", otherExamples)
        else:
            self.send("Selected Examples", None)
            self.send("Unselected Examples", None)
        
    def keyPressEvent(self, key):
        if key.key()==Qt.Key_Control:
            self.ctrlPressed=True
        else:
            OWWidget.keyPressEvent(self, key)

    def keyReleaseEvent(self, key):
        if key.key()==Qt.Key_Control:
            self.ctrlPressed=False
        else:
            OWWidget.keyReleaseEvent(self, key)
            
    def onStateChange(self, stateType, id, text):
        if stateType == "Warning":
            self.pathwayView._userMessage = text
            self.pathwayView.viewport().update()
            

if __name__=="__main__":
    app = QApplication(sys.argv)
    data = orange.ExampleTable("../../../../orange/doc/datasets/brown-selected.tab")
    w = OWKEGGPathwayBrowser()
    w.UpdateOrganismComboBox()
##    app.setMainWidget(w)
    w.show()
    w.SetData(data)
    app.exec_()
    w.saveSettings()
