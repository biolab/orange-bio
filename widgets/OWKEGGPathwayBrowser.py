"""<name>KEGG Pathway browser</name>
"""

import sys
import orange
import orngKEGG

from OWWidget import *
from qt import *

import OWGUI

from collections import defaultdict

class PathwayToolTip(QToolTip):
    def __init__(self, parent):
        QToolTip.__init__(self, parent)
        self.parent = parent

    def maybeTip(self, p):
        objs = self.parent.GetObjects(p.x() ,p.y())
        if objs:
            genes = map(self.parent.master.uniqueGenesDict.get, dict(objs).keys())
            text = "<br>".join(genes)
            self.tip(QRect(p.x()-2, p.y()-2, 4, 4), text)

class PathwayView(QScrollView):
    def __init__(self, master, *args, **argkw):
        QScrollView.__init__(self, *args, **argkw)
        self.master = master
        self.toolTip = PathwayToolTip(self)
        self.setHScrollBarMode(QScrollView.Auto)
        self.setVScrollBarMode(QScrollView.Auto)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        self.bbDict = {}
        self.pixmap = None
        self.image = None
        self.popup = QPopupMenu()
        self.popup.insertItem("View genes on KEGG website", 0, 0)
        self.popup.insertItem("View pathway on KEGG website", 1, 1)
        self.connect(self.popup, SIGNAL("activated ( int ) "), self.PopupAction)
        
    def SetPathway(self, pathway=None, objects=[]):
        self.pathway = pathway
        self.objects = objects
        if pathway:
            self.image = image = self.pathway.get_image()
            self.bbDict = self.pathway.get_bounding_box_dict()
            self.ShowImage()
##            image.save(self.pathway.local_database_path+"TmpPathwayImage.png")
##            self.pixmap = QPixmap(orngKEGG.default_database_path+"TmpPathwayImage.png")
##            w, h = image.size
##            self.resizeContents(w, h)
##            self.updateContents(self.contentsX(), self.contentsY() ,self.viewport().width(), self.viewport().height())
        else:
            self.bbDict = {}
            self.pixmap = None
            self.resizeContents(0,0)

    def ShowImage(self):
        if self.master.autoResize:
            import Image
            w, h = self.image.size
            self.resizeFactor = factor = min(self.viewport().width()/float(w), self.viewport().height()/float(h))
            image = self.image.resize((int(w*factor), int(h*factor)), Image.ANTIALIAS)
        else:
            image = self.image
            self.resizeFactor = 1
        image.save(self.pathway.local_database_path+"TmpPathwayImage.png")
        self.pixmap = QPixmap(self.pathway.local_database_path+"TmpPathwayImage.png")
        w, h = image.size
        self.resizeContents(w, h)
        self.updateContents(self.contentsX(), self.contentsY() ,self.viewport().width(), self.viewport().height())

    def drawContents(self, painter, cx=0, cy=0, cw=-1, ch=-1):
        QScrollView.drawContents(self, painter, cx, cy, cw, ch)
        if self.pixmap:
            cw = cw!=-1 and cw or self.viewport().width()
            ch = ch!=-1 and ch or self.viewport().height()
            painter.drawPixmap(cx, cy, self.pixmap, cx, cy, min(cw, self.pixmap.width()-cx), min(ch, self.pixmap.height()-cy))
            painter.save()

            painter.setPen(QPen(Qt.blue, 2, Qt.SolidLine))
            painter.setBrush(QBrush(Qt.NoBrush))
            for rect in reduce(lambda a,b:a.union(b), [bbList for id, bbList in self.bbDict.items() if id in self.objects], set()):
                x1, y1, x2, y2 = map(lambda x:int(self.resizeFactor*x), rect)
                painter.drawRect(x1+1, y1+1, x2-x1, y2-y1)
                
            painter.setPen(QPen(Qt.red, 2, Qt.SolidLine))
            for rect in self.master.selectedObjects.keys():
                x1, y1, x2, y2 = map(lambda x:int(self.resizeFactor*x), rect)
                painter.drawRect(x1+1, y1+1, x2-x1, y2-y1)
            painter.restore()

    def GetObjects(self, x, y):
        def _in(x, y, bb):
##            if bb[0]=="rect":
            x1, y1, x2, y2 = map(lambda x:int(self.resizeFactor*x), bb)
            return x>=x1 and y>=y1 and x<x2 and y<y2
##            else:
##                x1, y1, r = map(lambda x:int(self.resizeFactor*x), bb[1:])
##                return abs(x1-x)<=r and abs(y1-y)<=r
        x, y = self.viewportToContents(x, y)
        objs = []
        for id, bbList in self.bbDict.items():
            if id in self.objects:
                for bb in bbList:
                    if _in(x, y, bb):
                        objs.append((id, bb))
        return objs
    
    def viewportMouseMoveEvent(self, event):
        x, y = event.x(), event.y()
        objs = self.GetObjects(x, y)

    def viewportMousePressEvent(self, event):
        x, y = event.x(), event.y()
        old = set(self.master.selectedObjects.keys())
        objs = self.GetObjects(x, y)
        if event.button()==Qt.LeftButton:
            self.master.SelectObjects(objs)
            for rect in set(self.master.selectedObjects.keys()).union(old):
                x1, y1, x2, y2 = map(lambda x:int(self.resizeFactor*x), rect)
                self.updateContents(x1-1, y1-1, x2-x1+2, y2-y1+2)
        elif event.button()==Qt.RightButton:
            self.popup.objs = objs
            self.popup.setItemEnabled(0, bool(objs))
            self.popup.popup(self.mapToGlobal(event.pos()))
        else:
            QScrollView.viewportMousePressEvent(self, event)

    def resizeEvent(self, event):
        QScrollView.resizeEvent(self, event)
        if self.master.autoResize and self.image:
            self.ShowImage()

    def PopupAction(self, id):
        import webbrowser
        if id==0:
            genes = [s.split(":")[-1].strip() for s, t in self.popup.objs]
            address = "http://www.genome.jp/dbget-bin/www_bget?"+self.pathway.org+"+"+"+".join(genes)
        elif id==1:
##            genes = [s for s, t in self.popup.objs]
##            s = reduce(lambda s,g:s.union(self.master.org.get_enzymes_by_gene(g)), genes, set())
##            address = "http://www.genome.jp/dbget-bin/www_bget?enzyme+"+"+".join([e.split(":")[-1] for e in s])
            genes = [s.split(":")[-1].strip() for s, t in self.popup.objs]
            address = "http://www.genome.jp/dbget-bin/show_pathway?"+self.pathway.pathway_id.split(":")[-1]+"+"+"+".join(genes)
        try:
            webbrowser.open(address)
        except:
            pass
    
class OWKEGGPathwayBrowser(OWWidget):
    settingsList = ["organismIndex", "geneAttrIndex", "autoCommit", "autoResize", "useReference"]
    contextHandlers = {"":DomainContextHandler("",[ContextField("organismIndex", DomainContextHandler.Required + DomainContextHandler.IncludeMetaAttributes),
                                                   ContextField("geneAttrIndex", DomainContextHandler.Required + DomainContextHandler.IncludeMetaAttributes),
                                                   ContextField("useAttrNames", DomainContextHandler.Required + DomainContextHandler.IncludeMetaAttributes)])}
    def __init__(self, parent=None, signalManager=None, name="KEGG Pathway browser"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Examples", ExampleTable, self.SetData), ("Reference", ExampleTable, self.SetRefData)]
        self.outputs = [("Selected Examples", ExampleTable), ("Unselected Examples", ExampleTable)]
        self.organismIndex = 0
        self.geneAttrIndex = 0
        self.autoCommit = False
        self.autoResize = True
        self.useReference = False
        self.useAttrNames = False
        self.loadSettings()

        self.controlArea.setMaximumWidth(250)
        self.organismCodes = orngKEGG.KEGGInterfaceLocal().list_organisms().items()
        self.organismCodes.sort()
        items = [code+": "+desc for code, desc in self.organismCodes]
        self.organismCodes = [code for code, desc in self.organismCodes]
        cb = OWGUI.comboBox(self.controlArea, self, "organismIndex", box="Organism", items=items, callback=self.Update, addSpace=True)
        cb.setMaximumWidth(200)
        box = OWGUI.widgetBox(self.controlArea, "Gene attribure")
        self.geneAttrCombo = OWGUI.comboBox(box, self, "geneAttrIndex", callback=self.Update, addSpace=True)
        OWGUI.checkBox(box, self, "useAttrNames", "Use variable names", callback=self.UseAttrNamesCallback)
        self.geneAttrCombo.setDisabled(bool(self.useAttrNames))
        OWGUI.checkBox(self.controlArea, self, "useReference", "From signal", box="Reference", callback=self.Update)

        self.listView = QListView(self.controlArea)
        for header in ["Pathway", "P value", "Genes", "Reference"]:
            self.listView.addColumn(header)
        self.listView.setSelectionMode(QListView.Single)
        self.connect(self.listView, SIGNAL("selectionChanged ( QListViewItem * )"), self.UpdatePathwayView)

        self.pathwayLayout = QVBoxLayout(self.mainArea, QVBoxLayout.TopToBottom)
        self.pathwayView = PathwayView(self, self.mainArea)
        self.pathwayLayout.addWidget(self.pathwayView)

        OWGUI.checkBox(self.controlArea, self, "autoResize", "Resize to fit", box="Image", callback=lambda :self.pathwayView.image and self.pathwayView.ShowImage())

        box = OWGUI.widgetBox(self.controlArea, "Selection")
        OWGUI.checkBox(box, self, "autoCommit", "Commit on update")
        OWGUI.button(box, self, "Commit", callback=self.Commit)

        self.ctrlPressed=False
        self.selectedObjects = defaultdict(list)
        self.refData = None
        
    def SetData(self, data=None):
        self.closeContext()
        self.data = data
        if data:
            self.SetBestGeneAttrAndOrganism()
            self.openContext("", data)
            self.Update()
        else:
            self.listView.clear()
            self.selectedObjects = defaultdict(list)
            self.pathwayView.SetPathway(None)
            self.send("Selected Examples", None)
            self.send("Unselected Examples", None)

    def SetRefData(self, data=None):
        self.refData = data
        if self.useReference and self.data:
            self.Update()

    def UseAttrNamesCallback(self):
        self.geneAttrCombo.setDisabled(bool(self.useAttrNames))
        self.Update()

    def SetBestGeneAttrAndOrganism(self):
        self.geneAttrCandidates = self.data.domain.attributes + self.data.domain.getmetas().values()
        self.geneAttrCandidates = filter(lambda v:v.varType in [orange.VarTypes.Discrete ,orange.VarTypes.String], self.geneAttrCandidates)
        self.geneAttrCombo.clear()
        self.geneAttrCombo.insertStrList([var.name for var in self.geneAttrCandidates])
        data = self.data
        if len(data)>20:
            data = data.select(orange.MakeRandomIndices2(data, 20))
        from cPickle import load
        score = {}
        self.progressBarInit()
        attrNames = [str(v.name).strip() for v in self.data.domain.attributes]
        for i, org in enumerate(self.organismCodes):
            try:
                geneNames = load(open(orngKEGG.default_database_path+org+"_genenames.pickle"))
            except:
                continue
            for attr in self.geneAttrCandidates:
                vals = [str(e[attr]).strip() for e in data if not e[attr].isSpecial()]
                match = filter(lambda v:v in geneNames, vals)
                score[(attr, org)] = len(match)
            match = [v for v in attrNames if v in geneNames]
            score[("_var_names_", org)] = len(match)
            self.progressBarSet(i*100.0/len(self.organismCodes))
        self.progressBarFinished()
        score = [(s, attr, org) for (attr, org), s in score.items()]
        score.sort()
        if score[-1][1]=="_var_names_":
            self.useAttrNames = True
            self.geneAttrIndex = 0 #self.geneAttrCandidates.index(score[-2][1])
            self.geneAttrCombo.setDisabled(bool(self.useAttrNames))
        else:
            self.geneAttrIndex = self.geneAttrCandidates.index(score[-1][1])
        self.organismIndex = self.organismCodes.index(score[-1][2])
                
    def UpdateListView(self):
        self.listView.clear()
        pathways = self.pathways.items()
        allPathways = self.org.list_pathways()
        pathways.sort(lambda a,b:-cmp(a[1][1], b[1][1]))
        for id, (genes, p_value, ref) in pathways:
            item = QListViewItem(self.listView)
            item.setText(0, allPathways.get(id, id))
            item.setText(1, "%.4f" % p_value)
            item.setText(2, "%i of %i" %(len(genes), len(self.genes)))
            item.setText(3, "%i of %i" %(ref, len(self.referenceGenes)))
            item.pathway_id = id

    def UpdatePathwayView(self, item=None):
        self.selectedObjects = defaultdict(list)
        self.Commit()
        item = item and self.listView.selectedItem()
        if not item:
            self.pathwayView.SetPathway(None)
            return
        self.pathway = orngKEGG.KEGGPathway(item.pathway_id)
        self.pathwayView.SetPathway(self.pathway, self.pathways[item.pathway_id][0])
        
    def Update(self):
        self.error(0)
        if self.useAttrNames:
            genes = [str(v.name).strip() for v in self.data.domain.attributes]
        elif self.geneAttrCandidates:
            geneAttr = self.geneAttrCandidates[min(self.geneAttrIndex, len(self.geneAttrCandidates)-1)]
            genes = [str(e[geneAttr]) for e in self.data if not e[geneAttr].isSpecial()]
        else:
            self.error(0, "Cannot extact gene names from input")
            genes = []
        self.org = orngKEGG.KEGGOrganism(self.organismCodes[self.organismIndex])
        uniqueGenes, conflicting, unknown = self.org.get_unique_gene_ids(set(genes))
        if conflicting:
            print "Conflicting genes:", conflicting
        if unknown:
            print "Unknown genes:", unknown
        if self.useReference and self.refData:
            if self.useAttrNames:
                reference = [str(v.name).strip() for v in self.refData]
            else:
                geneAttr = self.geneAttrCandidates[min(self.geneAttrIndex, len(self.geneAttrCandidates)-1)]
                reference = [str(e[geneAttr]) for e in self.refData if not e[geneAttr].isSpecial()]
            uniqueRefGenes, conflicting, unknown = self.org.get_unique_gene_ids(set(reference))
            self.referenceGenes = reference = uniqueRefGenes.keys()
        else:
            self.referenceGenes = reference = self.org.get_genes()
        self.uniqueGenesDict = uniqueGenes
        self.genes = uniqueGenes.keys()
        self.revUniqueGenesDict = dict([(val, key) for key, val in self.uniqueGenesDict.items()])
        self.progressBarInit()
        self.pathways = self.org.get_enriched_pathways_by_genes(self.genes, reference, callback=self.progressBarSet)
        self.progressBarFinished()
        self.UpdateListView()
        self.UpdatePathwayView()

    def SelectObjects(self, objs):
        if (not self.selectedObjects or self.ctrlPressed) and not objs:
            return
        if self.ctrlPressed:
            for id, rect in objs:
                if id in self.selectedObjects[rect]:
                    self.selectedObjects[rect].pop(self.selectedObjects[rect].index(id))
                    if not self.selectedObjects[rect]:
                        del self.selectedObjects[rect]
                else:
                    self.selectedObjects[rect].append(id)
        else:
            self.selectedObjects.clear()
            for id, rect in objs:
                self.selectedObjects[rect].append(id)
        if self.autoCommit:
            self.Commit()

    def Commit(self):
        if self.useAttrNames:
            selectedGenes = reduce(lambda s,b:s.union(b), self.selectedObjects.values(), set())
            selectedVars = [self.data.domain[self.uniqueGenesDict[gene]] for gene in selectedGenes]
            newDomain = orange.Domain(selectedVars ,0)
            self.send("Selected Examples", orange.ExampleTable(newDomain, self.data))
        else:
            geneAttr = self.geneAttrCandidates[min(self.geneAttrIndex, len(self.geneAttrCandidates)-1)]
            selectedExamples = []
            otherExamples = []
            selectedGenes = reduce(lambda s,b:s.union(b), self.selectedObjects.values(), set())
            for ex in self.data:
                name = self.revUniqueGenesDict.get(str(ex[geneAttr]).strip(), None)
                if name and name in selectedGenes:
                    selectedExamples.append(ex)
                else:
                    otherExamples.append(ex)
            self.send("Selected Examples", selectedExamples and orange.ExampleTable(selectedExamples) or None)
            self.send("Unselected Examples", otherExamples and orange.ExampleTable(otherExamples) or None)
        
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

if __name__=="__main__":
    app = QApplication(sys.argv)
    data = orange.ExampleTable("../../orange/doc/datasets/brown-selected.tab")
    w = OWKEGGPathwayBrowser()
    app.setMainWidget(w)
    w.show()
    w.SetData(data)
    app.exec_loop()
    w.saveSettings()