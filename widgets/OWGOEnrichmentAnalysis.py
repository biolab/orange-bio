"""
<name>GO Enrichment Analysis</name>
<description>GO Enrichment Analysis</description>
<contact>Ales Erjavec</contact>
<icon>icons/GOTermFinder.png</icon>
<priority>103</priority>
"""

import go
import sys
import OWGUI

from OWWidget import *
from qt import *
from qttable import *

class TreeNode(object):
    def __init__(self, tuple, children):
        self.tuple = tuple
        self.children = children

def paintSection(self, painter, index, fr):
    print self
    if index!=5:
        QHeader.paintSection(self, painter, index, fr)
    else:
        pass

class MyListView(QListView):
    def __init__(self, parent):
        apply(QListView.__init__,(self, parent))

    #def paintEvent(self, event):
    
class OWGOEnrichmentAnalysis(OWWidget):
    settingsList=["annotationIndex", "useReferenceDataset", "aspectIndex", "geneAttrIndex",
                    "filterByNumOfInstances", "minNumOfInstances", "filterByPValue", "maxPValue", "selectionDirectAnnotation", "selectionDisjoint",
                    "selectionAddTermAsClass", "useAttrNames"]
    def __init__(self, parent=None, signalManager=None, name="GO Enrichment Analysis"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Cluster Examples", ExampleTable, self.SetClusterDataset, Default), ("Reference Examples", ExampleTable, self.SetReferenceDataset, Single + NonDefault)] #, ("Structured Data", DataFiles, self.chipdata, Single + NonDefault)]
        self.outputs = [("Selected Examples", ExampleTable, Default), ("Unselected Examples", ExampleTable, Default), ("Example With Unknown Genes", ExampleTable, Default)] #, ("Selected Structured Data", DataFiles, Single + NonDefault)]

        self.annotationIndex = 0
        self.useReferenceDataset  = 0
        self.aspectIndex = 0
        self.geneAttrIndex = 0
        self.useAttrNames = False
        self.filterByNumOfInstances = False
        self.minNumOfInstances = 1
        self.filterByPValue = True
        self.maxPValue = 0.1
        self.selectionDirectAnnotation = 0
        self.selectionDisjoint = 0
        self.selectionAddTermAsClass = 0
        # check usage of all evidences
        for etype in go.evidenceTypesOrdered:
            varName = "useEvidence"+etype 
##            self.settingsList.append( varName)
            code = compile("self.%s = True" % (varName), ".", "single")
            exec(code)

        self.updateGODataBaseWidget = None
        try:
            from OWUpdateGODataBase import OWUpdateGODataBase
            self.updateGODataBaseWidget = OWUpdateGODataBase()
            self.connect(w, PYSIGNAL("closed()"), self.UpdateAnnotationComboBox)
        except:
            pass
        
        self.annotationCodes = go.listDownloadedOrganisms()
        #############
        ##GUI
        #############
        self.tabs = QTabWidget(self.controlArea, 'tabWidget')
        ##Input tab
        self.inputTab = QVGroupBox(self)
        box = OWGUI.widgetBox(self.inputTab, "Organism annotation", addSpace=True)
        self.annotationComboBox = OWGUI.comboBox(box, self, "annotationIndex", items = self.annotationCodes, callback=self.SetAnnotationCallback)
        #box = OWGUI.widgetBox(box, "Evidence codes in annotation", addSpace=True)
        box = OWGUI.widgetBox(self.inputTab, "Evidence codes in annotation", addSpace=True)
        box.setMaximumWidth(150)
        self.evidenceCheckBoxDict = {}
        for etype in go.evidenceTypesOrdered:
            self.evidenceCheckBoxDict[etype] = OWGUI.checkBox(box, self, "useEvidence"+etype, etype, callback=self.Update)
        OWGUI.radioButtonsInBox(self.inputTab, self, "useReferenceDataset", ["Annotation", "Signal"], box="Reference From", callback=self.Update)
        OWGUI.radioButtonsInBox(self.inputTab, self, "aspectIndex", ["Biological process", "Cellular component", "Molecular function"], box="Aspect", callback=self.Update)
        self.geneAttrIndexCombo = OWGUI.comboBox(self.inputTab, self, "geneAttrIndex", box="Gene attribute", callback=self.Update)
        OWGUI.checkBox(self.geneAttrIndexCombo.box, self, "useAttrNames", "Use attribute names", callback=self.SetUseAttrNamesCallback)
        self.tabs.insertTab(self.inputTab, "Input")
        box = OWGUI.widgetBox(self.inputTab, "GO update")
        b = OWGUI.button(box, self, "Update", callback = self.UpdateGOAndAnnotation)
        box.setMaximumWidth(150)
        
        ##Filter tab
        self.filterTab = QVGroupBox(self)
        box = OWGUI.widgetBox(self.filterTab, "Filter GO Term Nodes", addSpace=True)
        OWGUI.checkBox(box, self, "filterByNumOfInstances", "Number of instances", callback=self.FilterAndDisplayGraph)
        OWGUI.qwtHSlider(box, self, 'minNumOfInstances', label='#:', labelWidth=33, minValue=1, maxValue=1000, step=1.0, precision=1, ticks=0, maxWidth=80, callback=self.FilterAndDisplayGraph)
        OWGUI.checkBox(box, self, "filterByPValue", "p-value",callback=self.FilterAndDisplayGraph)
        OWGUI.qwtHSlider(box, self, 'maxPValue', label='p:', labelWidth=33, minValue=0, maxValue=1, step=0.001, precision=3, ticks=0, maxWidth=80, callback=self.FilterAndDisplayGraph)
        self.tabs.insertTab(self.filterTab, "Filter")
        
        ##Select tab
        self.selectTab=QVGroupBox(self)
        #box = OWGUI.widgetBox(self.selectTab, "Annotated genes", addSpace=True)
        box = OWGUI.radioButtonsInBox(self.selectTab, self, "selectionDirectAnnotation", ["Directly or Indirectly", "Directly"], box="Annotated genes", callback=self.ExampleSelection)
        box = OWGUI.widgetBox(self.selectTab, "Output", addSpace=True)
        OWGUI.checkBox(box, self, "selectionDisjoint", "Disjoint/Inclusive", callback=self.ExampleSelection)
        OWGUI.checkBox(box, self, "selectionAddTermAsClass", "Add GO Term as class", callback=self.ExampleSelection)
        self.tabs.insertTab(self.selectTab, "Select")

        # ListView for DAG, and table for significant GOIDs
        self.DAGcolumns = ['GO term', 'Cluster frequency', 'Reference frequency', 'p value', 'Genes', 'Enrichment']
        self.layout=QVBoxLayout(self.mainArea)
        self.splitter = QSplitter(QSplitter.Vertical, self.mainArea)
        self.layout.add(self.splitter)

        # list view
        self.listView = MyListView(self.splitter)
        self.listView.setMultiSelection(1)
        self.listView.setAllColumnsShowFocus(1)
        self.listView.addColumn(self.DAGcolumns[0])
        self.listView.setColumnWidth(0, 300)
        self.listView.setColumnWidthMode(0, QListView.Manual)
        self.listView.setColumnAlignment(0, QListView.AlignLeft)
        self.listView.setSorting(-1)
        for dagColumnTitle in self.DAGcolumns[1:]:
            col = self.listView.addColumn(dagColumnTitle)
            self.listView.setColumnWidth(col, 100)
            self.listView.setColumnWidthMode(col, QListView.Manual)
            self.listView.setColumnAlignment(col, QListView.AlignCenter)
        self.connect(self.listView, SIGNAL("selectionChanged()"), self.ViewSelectionChanged)
        self.listView.setRootIsDecorated (True)

        # table of significant GO terms
        self.sigTermsTable = QTable(self.splitter)
        self.sigTermsTable.setNumCols(6)
        self.sigTermsTable.setNumRows(4)
        ## hide the vertical header
        self.sigTermsTable.verticalHeader().hide()
        self.sigTermsTable.setLeftMargin(0)
        self.sigTermsTable.setSelectionMode(QTable.Multi)
        self.sigTermsTable.setColumnWidth(0, 300)
        for col in range(1, self.sigTermsTable.numCols()):
            self.sigTermsTable.setColumnWidth(col, 100)
        self.header = self.sigTermsTable.horizontalHeader()
        for i in range(len(self.DAGcolumns)):
            self.header.setLabel(i, self.DAGcolumns[i])
        self.connect(self.sigTermsTable, SIGNAL("selectionChanged()"), self.TableSelectionChanged)
        self.splitter.show()

        self.sigTableTermsSorted = []
        self.graph = {}
        #header = self.listView.header()
        #from new import instancemethod
        #print(dir(header))
        #header.paintSection = instancemethod(paintSection, header, type(header))
        self.resize(900, 800)
        
    def SetAnnotationCallback(self):
        self.LoadAnnotation()
        if self.clusterDataset:
            self.SetClusterDataset(self.clusterDataset)

    def UpdateSelectedEvidences(self):
        self.SetClusterDataset(self.clusterDataset)

    def SetReferenceCallback(self):
        self.SetClusterDataset(self.clusterDataset)

    def SetAspectCallback(self):
        self.SetClusterDataset(self.clusterDataset)

    def SetUseAttrNamesCallback(self):
        self.geneAttrIndexCombo.setDisabled(bool(self.useAttrNames))
        self.Update()

    def Update(self):
        self.SetClusterDataset(self.clusterDataset)

    def UpdateGOAndAnnotation(self):
        from OWUpdateGODataBase import OWUpdateGODataBase
        w = OWUpdateGODataBase()
        w.show()
        self.connect(w, PYSIGNAL("closed()"), self.UpdateAnnotationComboBox)

    def UpdateAnnotationComboBox(self):
        curr = self.annotationCodes[self.annotationIndex]
        self.annotationCodes = go.listDownloadedOrganisms()
        index = self.annotationCodes.index(curr)
        self.annotationComboBox.clear()
        self.annotationComboBox.insertStrList(self.annotationCodes)
        self.annotationComboBox.setCurrentItem(index)
        #print go.getDataDir()

    def SetGenesComboBox(self):
        self.candidateGeneAttrs = self.clusterDataset.domain.attributes + self.clusterDataset.domain.getmetas().values()
        self.candidateGeneAttrs = filter(lambda v: v.varType==orange.VarTypes.String or v.varType==orange.VarTypes.Other or v.varType==orange.VarTypes.Discrete, self.candidateGeneAttrs)
        self.geneAttrIndexCombo.clear()
        self.geneAttrIndexCombo.insertStrList([a.name for a in  self.candidateGeneAttrs])

    def FindBestGeneAttrAndOrganism(self):
        organismGenes = dict([(o,set(go.getCachedGeneNames(o))) for o in self.annotationCodes])
        candidateGeneAttrs = self.clusterDataset.domain.attributes + self.clusterDataset.domain.getmetas().values()
        candidateGeneAttrs = filter(lambda v: v.varType==orange.VarTypes.String or v.varType==orange.VarTypes.Other or v.varType==orange.VarTypes.Discrete, candidateGeneAttrs)
        attrNames = [v.name for v in self.clusterDataset.domain.variables]
        cn = {}
        for attr in candidateGeneAttrs:
            vals = [str(e[attr]) for e in self.clusterDataset]
            for organism, s in organismGenes.items():
                l = filter(lambda a: a in s, vals)
                cn[(attr,organism)] = len(l)
        for organism, s in organismGenes.items():
            l = filter(lambda a: a in s, attrNames)
            cn[("_var_names_", organism)] = len(l)
            
        cn = cn.items()
        cn.sort(lambda a,b:-cmp(a[1],b[1]))
        bestAttr, organizm = cn[0][0]
        self.annotationIndex = self.annotationCodes.index(organizm)
        if bestAttr=="_var_names_":
            self.useAttrNames = True
            self.geneAttrIndexCombo.setDisabled(True)
            self.geneAttrIndex = 0
        else:
            self.useAttrNames = False
            self.geneAttrIndexCombo.setDisabled(False)
            self.geneAttrIndex = candidateGeneAttrs.index(bestAttr)
    
    def SetClusterDataset(self, data=None):
        self.clusterDataset = data
        if data:
            self.SetGenesComboBox()
            if not go.loadedGO:
                self.LoadGO()
            self.FindBestGeneAttrAndOrganism()
            if not go.loadedAnnotation:
                self.LoadAnnotation()
            self.FilterUnknownGenes()
            graph = self.Enrichment()
            #print graph
            self.SetGraph(graph)
        else:
            self.ClearGraph()
            self.send("Selected Examples", None)
            self.send("Unselected Examples", None)
            self.send("Example With Unknown Genes", None)

    def SetReferenceDataset(self, data=None):
        self.referenceDataset=data
        if data and self.useReferenceDataset:
            graph = self.Enrichment(self.data)
            self.SetGraph(graph)

    def FilterUnknownGenes(self):
        if not self.useAttrNames:
            geneAttr = self.candidateGeneAttrs[self.geneAttrIndex]
            examples = []
            for ex in self.clusterDataset:
                if str(ex[geneAttr]) not in go.loadedAnnotation.aliasMapper:
                    examples.append(ex)
            self.send("Example With Unknown Genes", examples and orange.ExampleTable(examples) or None)
        else:
            self.send("Example With Unknown Genes", None)

    def LoadGO(self):
        self.progressBarInit()
        go.loadGO(progressCallback=self.progressBarSet)
        self.progressBarFinished()
        
    def LoadAnnotation(self):
        self.progressBarInit()
        go.loadAnnotation(self.annotationCodes[self.annotationIndex], progressCallback=self.progressBarSet)
        self.progressBarFinished()
        count = dict([(etype, 0) for etype in go.evidenceTypesOrdered])
        geneSets = dict([(etype, set()) for etype in go.evidenceTypesOrdered])
        for anno in go.loadedAnnotation.annotationList:
            count[anno.evidence]+=1
            geneSets[anno.evidence].add(anno.geneName)
        for etype in go.evidenceTypesOrdered:
            self.evidenceCheckBoxDict[etype].setEnabled(bool(count[etype]))
            self.evidenceCheckBoxDict[etype].setText(etype+": %i annots(%i genes)" % (count[etype], len(geneSets[etype])))
        
    def Enrichment(self):
        if self.useAttrNames:
            clusterGenes = [v.name for v in self.clusterDataset.domain.variables]
        else:
            geneAttr = self.candidateGeneAttrs[self.geneAttrIndex]
            clusterGenes = [str(ex[geneAttr]) for ex in self.clusterDataset if not ex[geneAttr].isSpecial()]
        self.clusterGenes = clusterGenes = filter(lambda g: g in go.loadedAnnotation.aliasMapper, clusterGenes)
        referenceGenes = None
        if self.useReferenceDataset:
            try:
                if self.useAttrNames:
                    referenceGenes = [v.name for v in self.referenceDataset.domain.variables]
                else:
                    referenceGenes = [str(ex[geneAttr]) for ex in self.referenceDataset if not ex[geneAttr].isSpecial()]
                referenceGenes = filter(lambda g: g in go.loadedAnnotation.aliasMapper, referenceGenes)
                self.information()
            except Exception, er:
                self.information(str(er)+" Using the annotation for reference")
        else:
            self.information()
            self.referenceGenes = go.loadedAnnotation.geneNames
        evidences = []
        for etype in go.evidenceTypesOrdered:
            if getattr(self, "useEvidence"+etype):
                evidences.append(etype)
        aspect = ["P", "F", "C"][self.aspectIndex]
        self.progressBarInit()
        if clusterGenes:
            self.terms = terms = go.GOTermFinder(clusterGenes, referenceGenes, evidences, aspect=aspect, progressCallback=self.progressBarSet)
        else:
            self.terms = terms = {}
        self.progressBarFinished()
        self.treeStructDict = {}
        ids = self.terms.keys()
        for term in self.terms:
            self.treeStructDict[term] = TreeNode(self.terms[term], filter(lambda t:term in go.loadedGO.termDict[t].parents, ids))
            if not go.loadedGO.termDict[term].parents:
                self.treeStructRootKey = term
        return terms
        
    def FilterGraph(self, graph):
        if self.filterByPValue:
            graph = go.filterByPValue(graph, self.maxPValue)
        if self.filterByNumOfInstances:
            graph = dict(filter(lambda (id,(genes, p, rc)):len(genes)>=self.minNumOfInstances, graph.items()))
        return graph

    def FilterAndDisplayGraph(self):
        self.graph = self.FilterGraph(self.originalGraph)
        self.ClearGraph()
        self.DisplayGraph()

    def SetGraph(self, graph=None):
        self.originalGraph = graph
        if graph:
            self.FilterAndDisplayGraph()
        else:
            print self.annotationCodes[self.annotationIndex], self.candidateGeneAttrs[self.geneAttrIndex]
            self.graph = {}
            self.ClearGraph()

    def ClearGraph(self):
        self.listView.clear()
        self.listViewItems=[]
        self.sigTermsTable.setNumRows(0)
        #self.sigTableItems=[]

    def DisplayGraph(self):
        fromParentDict = {}
        self.termListViewItemDict = {}
        self.listViewItems=[]
        enrichment = lambda t:float(len(t[0])) / t[2] * (float(len(self.referenceGenes))/len(self.clusterGenes))
        maxFoldEnrichment = max([enrichment(term) for term in self.graph.values()] or [1])
        def addNode(term, parent, parentDisplayNode):
            if (parent, term) in fromParentDict:
                return
            if term in self.graph:
                displayNode = MyListViewItem(parentDisplayNode)
                displayNode.setText(0, go.loadedGO.termDict[term].name)
                displayNode.setText(1, str(len(self.graph[term][0])))
                displayNode.setText(2, str(self.graph[term][2]))
                displayNode.setText(3, "%.4f" % self.graph[term][1])
                displayNode.setText(4, ", ".join(self.graph[term][0]))
                displayNode.setText(5, "%.4f" % (enrichment(self.graph[term])/maxFoldEnrichment)) #(float(len(self.graph[term][0]))/self.graph[term][2]))
                displayNode.setOpen(True)
                displayNode.term=term
                self.listViewItems.append(displayNode)
                if term in self.termListViewItemDict:
                    self.termListViewItemDict[term].append(displayNode)
                else:
                    self.termListViewItemDict[term] = [displayNode]
                fromParentDict[(parent, term)] = True
                parent = term
            else:
                displayNode = parentDisplayNode
            
            for c in self.treeStructDict[term].children:
                addNode(c, parent, displayNode)
        addNode(self.treeStructRootKey, None, self.listView)

        terms = self.graph.items()
        terms.sort(lambda a,b:cmp(a[1][1],b[1][1]))
        self.sigTableTermsSorted = [t[0] for t in terms]
        self.sigTermsTable.setNumRows(len(terms))
        for i, (id, (genes, p_value, refCount)) in enumerate(terms):
            text = [go.loadedGO.termDict[id].name, str(len(genes)), str(refCount), "%.4f" % p_value, " ,".join(genes), "%.2f" % enrichment((genes, p_value, refCount))]
            for j,t in enumerate(text):
                self.sigTermsTable.setText(i, j, t)

    def ViewSelectionChanged(self):
        selected = filter(lambda lvi: lvi.isSelected(), self.listViewItems)
        self.selectedTerms = set([lvi.term for lvi in selected])
        self.ExampleSelection()
        
    def TableSelectionChanged(self):
        self.selectedTerms=[]
        for i, term in enumerate(self.sigTableTermsSorted):
            selected = self.sigTermsTable.isRowSelected(i, False)
            if selected:
                self.selectedTerms.append(term)
            for lvi in self.termListViewItemDict[term]:
                try:
                    lvi.setSelected(selected)
                    self.listView.repaintItem(lvi)
                    if selected: lvi.setOpen(True)
                except RuntimeError:    ##Underlying C/C++ object deleted (why??)
                    pass
        self.listView.triggerUpdate()
        self.ExampleSelection()
            
    
    def ExampleSelection(self):
        selectedExamples = []
        unselectedExamples = []
        selectedGenes = []
        if self.selectionDirectAnnotation:
            s = filter(lambda anno: anno.GOId in self.selectedTerms, go.loadedAnnotation.annotationList)
            selectedGenes = [anno.geneName for anno in s]
        else:        
            map(selectedGenes.extend, [v[0] for id, v in self.graph.items() if id in self.selectedTerms])
            
        if self.selectionDisjoint:
            count = dict([(g, 0) for g in self.clusterGenes])
            for term in self.selectedTerms:
                for g in self.graph[term][0]:
                    count[g]+=1
            selectedGenes = [gene for gene, c in count.items() if c==1 and gene in selectedGenes]

        if self.useAttrNames:
            vars = [self.clusterDataset.domain[gene] for gene in set(selectedGenes)]
            newDomain = orange.Domain(vars, 0)
            self.send("Selected Examples", orange.ExampleTable(newDomain, self.clusterDataset))
            self.send("Unselected Examples", None)
        else:
            geneAttr = self.candidateGeneAttrs[self.geneAttrIndex]            
            newClass = orange.EnumVariable("GO Term", values=list(self.selectedTerms))
            newDomain = orange.Domain(self.clusterDataset.domain.variables, newClass)
            for ex in self.clusterDataset:
                if not ex[geneAttr].isSpecial() and str(ex[geneAttr]) in selectedGenes:
                    if self.selectionDisjoint and self.selectionAddTermAsClass:
                        c = filter(lambda term: str(ex[geneAttr]) in self.graph[term][0], self.selectedTerms)[0]
                        ex =  orange.Example(newDomain, ex)
                        ex.setclass(newClass(c))
                    selectedExamples.append(ex)
                else:
                    unselectedExamples.append(ex)
            self.send("Selected Examples", selectedExamples and orange.ExampleTable(selectedExamples) or None)
            self.send("Unselected Examples", unselectedExamples and orange.ExampleTable(unselectedExamples) or None)
            
class MyListViewItem(QListViewItem):
    enrichmentColumn = 5
    def paintCell(self, painter, colorgroup, column, width, align):
        if column!=self.enrichmentColumn:
            QListViewItem.paintCell(self, painter, colorgroup, column, width, align)
        else:
            f = float(str(self.text(self.enrichmentColumn)))
            painter.setBrush(QBrush(Qt.white, QBrush.SolidPattern))
            painter.drawRect(0, 0, width-1, self.height()-1)
            painter.setBrush(QBrush(Qt.blue, QBrush.SolidPattern))
            painter.drawRect((1-f)*(width-1), 0, width-1, self.height()-1)

    def width(self):
        return 100
    
if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    w=OWGOEnrichmentAnalysis()
    data = orange.ExampleTable("../../doc/datasets/brown-selected.tab")
    w.SetClusterDataset(data)
    app.setMainWidget(w)
    w.show()
    app.exec_loop()
    w.saveSettings()
        
