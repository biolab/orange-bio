"""
<name>GO Enrichment Analysis</name>
<description>Enrichment analysis for Gene Ontology terms.</description>
<contact>Ales Erjavec</contact>
<icon>icons/GOTermFinder.png</icon>
<priority>103</priority>
"""

import obiGO
import obiKEGG
import obiProb
import obiTaxonomy
import sys, os, tarfile, math
import OWGUI
import orngServerFiles

from os.path import join as p_join
from OWWidget import *
from collections import defaultdict
from obiGeneMatch import GeneMatchMk2
from functools import partial

dataDir = orngServerFiles.localpath("GO")

def listAvailable():
    import orngServerFiles
    files = orngServerFiles.listfiles("GO")
    ret = {}
    for file in files:
        tags = orngServerFiles.info("GO", file)["tags"]
        td = dict([tuple(tag.split(":")) for tag in tags if tag.startswith("#") and ":" in tag])
        if "association" in file.lower():
            ret[td.get("#organism", file)] = file
    orgMap = {"352472":"44689"}
    essential = [obiGO.from_taxid(orgMap.get(id, id)).pop() for id in obiTaxonomy.essential_taxids()]
    essentialNames = [obiTaxonomy.name(orgMap.get(id, id)) for id in obiTaxonomy.essential_taxids()]
    ret.update(zip(essentialNames, essential))
    return ret

def getOrgFileName(org):
    import orngServerFiles
    files = orngServerFiles.listfiles("go")
    return [f for f in files if org in f].pop()

class TreeNode(object):
    def __init__(self, tuple, children):
        self.tuple = tuple
        self.children = children

class GOTreeWidget(QTreeWidget):
    def contextMenuEvent(self, event):
        QTreeWidget.contextMenuEvent(self, event)
##        print event.x(), event.y()
        term = self.itemAt(event.pos()).term
        self._currMenu = QMenu()
        self._currAction = self._currMenu.addAction("View term on AmiGO website")
##        self.connect(self, SIGNAL("triggered(QAction*)"), partial(self.BrowserAction, term))
        self.connect(self._currAction, SIGNAL("triggered()"), lambda :self.BrowserAction(term))
        self._currMenu.popup(event.globalPos())

    def BrowserAction(self, term):
        import webbrowser
        webbrowser.open("http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term="+term)

        
    
class OWGOEnrichmentAnalysis(OWWidget):
    settingsList=["annotationIndex", "useReferenceDataset", "aspectIndex", "geneAttrIndex",
                    "filterByNumOfInstances", "minNumOfInstances", "filterByPValue", "maxPValue", "selectionDirectAnnotation", "selectionDisjoint", "selectionType",
                    "selectionAddTermAsClass", "useAttrNames", "probFunc", "useFDR"]
    contextHandlers = {"": DomainContextHandler("", ["geneAttrIndex", "useAttrNames", "annotationIndex"], matchValues=1)}
    def __init__(self, parent=None, signalManager=None, name="GO Enrichment Analysis"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Cluster Examples", ExampleTable, self.SetClusterDataset, Default), ("Reference Examples", ExampleTable, self.SetReferenceDataset, Single + NonDefault)] #, ("Structured Data", DataFiles, self.chipdata, Single + NonDefault)]
        self.outputs = [("Selected Examples", ExampleTable, Default), ("Unselected Examples", ExampleTable, Default), ("Example With Unknown Genes", ExampleTable, Default)] #, ("Selected Structured Data", DataFiles, Single + NonDefault)]

        self.annotationIndex = 0
        self.autoFindBestOrg = False
        self.useReferenceDataset  = 0
        self.aspectIndex = 0
        self.geneAttrIndex = 0
        self.useAttrNames = False
        self.filterByNumOfInstances = False
        self.minNumOfInstances = 1
        self.filterByPValue = True
        self.maxPValue = 0.1
        self.probFunc = 0
        self.useFDR = True
        self.selectionDirectAnnotation = 0
        self.selectionDisjoint = 0
        self.selectionAddTermAsClass = 0
        self.selectionChanging = 0
        
        self.loadSettings()
        
        # check usage of all evidences
        for etype in obiGO.evidenceTypesOrdered:
            varName = "useEvidence" + etype 
##            self.settingsList.append( varName)
            code = compile("self.%s = True" % (varName), ".", "single")
            exec(code)
        self.annotationFiles = listAvailable()
        self.annotationCodes = sorted(self.annotationFiles.keys())
        if not self.annotationCodes:
            self.error(0, "No downloaded annotations!!\nUse the Update Genomics Databases widget and download annotations for at least one organism!")
        else:
            self.error(0)
        #############
        ##GUI
        #############
        self.tabs = OWGUI.tabWidget(self.controlArea)
        ##Input tab
        self.inputTab = OWGUI.createTabPage(self.tabs, "Input")
        box = OWGUI.widgetBox(self.inputTab, "Info")
        self.infoLabel = OWGUI.widgetLabel(box, "No data on input\n")
        OWGUI.button(box, self, "Ontology/Annotation Info", callback=self.ShowInfo, tooltip="Show information on loaded ontology and annotations")
        box = OWGUI.widgetBox(self.inputTab, "Organism", addSpace=True)
        self.annotationComboBox = OWGUI.comboBox(box, self, "annotationIndex", items = self.annotationCodes, callback=self.SetAnnotationCallback, tooltip="Select organism")
        self.geneAttrIndexCombo = OWGUI.comboBox(self.inputTab, self, "geneAttrIndex", box="Gene names", callback=self.Update, tooltip="Use this attribute to extract gene names from input data")
        OWGUI.checkBox(self.geneAttrIndexCombo.box, self, "useAttrNames", "Use data attributes names", disables=[(-1, self.geneAttrIndexCombo)], callback=self.SetUseAttrNamesCallback, tooltip="Use attribute names for gene names")
        
        self.referenceRadioBox = OWGUI.radioButtonsInBox(self.inputTab, self, "useReferenceDataset", ["Entire genome", "Reference set (input)"], tooltips=["Use entire genome for reference", "Use genes from Referece Examples input signal as reference"], box="Reference", callback=self.Update)
        self.referenceRadioBox.buttons[1].setDisabled(True)
        OWGUI.radioButtonsInBox(self.inputTab, self, "aspectIndex", ["Biological process", "Cellular component", "Molecular function"], box="Aspect", callback=self.Update)
        self.geneAttrIndexCombo.setDisabled(bool(self.useAttrNames))
##        self.geneInfoLabel = OWGUI.label(self.geneAttrIndexCombo.box, self, "0 genes on input signal")
       
##        box = OWGUI.widgetBox(self.inputTab, "GO update")
##        b = OWGUI.button(box, self, "Update", callback = self.UpdateGOAndAnnotation)
##        box.setMaximumWidth(150)
        
        ##Filter tab
        self.filterTab = OWGUI.createTabPage(self.tabs, "Filter")
        box = OWGUI.widgetBox(self.filterTab, "Filter GO Term Nodes", addSpace=True)
        OWGUI.checkBox(box, self, "filterByNumOfInstances", "Genes", callback=self.FilterAndDisplayGraph, tooltip="Filter by number of input genes mapped to a term")
##        OWGUI.qwtHSlider(box, self, 'minNumOfInstances', label='#:', labelWidth=5, minValue=1, maxValue=100, step=1.0, precision=1, ticks=0, maxWidth=60, callback=self.FilterAndDisplayGraph)
        OWGUI.spin(OWGUI.indentedBox(box), self, 'minNumOfInstances', 1, 100, step=1, label='#:', labelWidth=15, callback=self.FilterAndDisplayGraph, callbackOnReturn=True, tooltip="Min. number of input genes mapped to a term")
        OWGUI.checkBox(box, self, "filterByPValue", "Significance",callback=self.FilterAndDisplayGraph, tooltip="Filter by term p-value")
##        OWGUI.qwtHSlider(box, self, 'maxPValue', label='p:', labelWidth=5, minValue=0.001, maxValue=1, step=0.001, precision=3, ticks=0, logarithmic=True, maxWidth=60, callback=self.FilterAndDisplayGraph)
        OWGUI.doubleSpin(OWGUI.indentedBox(box), self, 'maxPValue', 1e-8, 1, step=1e-8,  label='p:', labelWidth=15, callback=self.FilterAndDisplayGraph, callbackOnReturn=True, tooltip="Max term p-value")
        box = OWGUI.widgetBox(box, "Significance test")
        OWGUI.radioButtonsInBox(box, self, "probFunc", ["Binomial", "Hypergeometric"], tooltips=["Use binomial distribution test", "Use hypergeometric distribution test"], callback=self.Update)
        OWGUI.checkBox(box, self, "useFDR", "Use FDR (False Discovery Rate)", callback=self.Update, tooltip="Use False Discovery Rate correction")
        box = OWGUI.widgetBox(self.filterTab, "Evidence codes in annotation", addSpace=True)
##        box.setMaximumWidth(150)
        self.evidenceCheckBoxDict = {}
        for etype in obiGO.evidenceTypesOrdered:
            self.evidenceCheckBoxDict[etype] = OWGUI.checkBox(box, self, "useEvidence"+etype, etype, callback=self.UpdateSelectedEvidences, tooltip=obiGO.evidenceTypes[etype])
        
        ##Select tab
        self.selectTab = OWGUI.createTabPage(self.tabs, "Select")
        #box = OWGUI.widgetBox(self.selectTab, "Annotated genes", addSpace=True)
        box = OWGUI.radioButtonsInBox(self.selectTab, self, "selectionDirectAnnotation", ["Directly or Indirectly", "Directly"], box="Annotated genes", callback=self.ExampleSelection)
        box = OWGUI.widgetBox(self.selectTab, "Output", addSpace=True)
##        OWGUI.checkBox(box, self, "selectionDisjoint", "Disjoint/Inclusive", callback=self.ExampleSelection)
        OWGUI.radioButtonsInBox(box, self, "selectionDisjoint", btnLabels=["All selected genes", "Term-specific genes"], tooltips=["Outputs genes annotated to all selected GO terms", "Outputs genes that appear in only one of selected GO terms"], callback=self.ExampleSelection)
        OWGUI.checkBox(box, self, "selectionAddTermAsClass", "Add GO Term as class", callback=self.ExampleSelection)

        # ListView for DAG, and table for significant GOIDs
        self.DAGcolumns = ['GO term', 'Cluster', 'Reference', 'p value', 'Genes', 'Enrichment']
        #self.layout=QVBoxLayout(self.mainArea)
        self.splitter = QSplitter(Qt.Vertical, self.mainArea)
        self.mainArea.layout().addWidget(self.splitter)

        # list view
        self.listView = GOTreeWidget(self.splitter)
        self.listView.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(len(self.DAGcolumns))
        self.listView.setHeaderLabels(self.DAGcolumns)
        
        self.listView.header().setClickable(True)
        self.listView.header().setSortIndicatorShown(True)
        self.listView.setSortingEnabled(True)
        self.listView.setItemDelegateForColumn(5, EnrichmentColumnItemDelegate(self))
        self.listView.setRootIsDecorated(True)

        
        self.connect(self.listView, SIGNAL("itemSelectionChanged()"), self.ViewSelectionChanged)
        
        # table of significant GO terms
        self.sigTerms = QTreeWidget(self.splitter)
        self.sigTerms.setColumnCount(len(self.DAGcolumns))
        self.sigTerms.setHeaderLabels(self.DAGcolumns)
        self.sigTerms.setSortingEnabled(True)
        self.sigTerms.setSelectionMode(QAbstractItemView.ExtendedSelection)
        
        self.connect(self.sigTerms, SIGNAL("itemSelectionChanged()"), self.TableSelectionChanged)
        self.splitter.show()

        self.sigTableTermsSorted = []
        self.graph = {}
        
        self.loadedAnnotationCode = "---"
        
        self.inputTab.layout().addStretch(1)
        self.filterTab.layout().addStretch(1)
        self.selectTab.layout().addStretch(1)
        
        self.resize(900, 800)

        self.keggOrg = None
        self.clusterDataset = None
        self.referenceDataset = None
        self.ontology = None
        self.annotations = None
        self.probFunctions = [obiProb.Binomial(), obiProb.Hypergeometric()]
        
    def SetAnnotationCallback(self):
        self.LoadAnnotation()
        if self.clusterDataset:
            self.FilterUnknownGenes()
            graph = self.Enrichment()
            self.SetGraph(graph)

    def UpdateSelectedEvidences(self):
        if self.clusterDataset:
            self.FilterUnknownGenes()
            graph = self.Enrichment()
            self.SetGraph(graph)

    def SetReferenceCallback(self):
        if self.clusterDataset:
            self.FilterUnknownGenes()
            graph = self.Enrichment()
            self.SetGraph(graph)

    def SetAspectCallback(self):
        if self.clusterDataset:
            self.FilterUnknownGenes()
            graph = self.Enrichment()
            self.SetGraph(graph)

    def SetUseAttrNamesCallback(self):
##        self.geneAttrIndexCombo.setDisabled(bool(self.useAttrNames))
        self.Update()

    def Update(self):
        if self.clusterDataset:
            self.FilterUnknownGenes()
            graph = self.Enrichment()
            self.SetGraph(graph)

    def UpdateGOAndAnnotation(self, tags=[]):
        from OWUpdateGenomicsDatabases import OWUpdateGenomicsDatabases
        w = OWUpdateGenomicsDatabases(parent = self, searchString=" ".join(tags))
        w.setModal(True)
        w.show()
        self.UpdateAnnotationComboBox()
##        self.connect(w, SIGNAL("closed()"), self.UpdateAnnotationComboBox)

    def UpdateAnnotationComboBox(self):
        if self.annotationCodes:
            curr = self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
        else:
            curr = None
        self.annotationFiles = listAvailable()
        self.annotationCodes = self.annotationFiles.keys()
        index = curr and self.annotationCodes.index(curr) or 0
        self.annotationComboBox.clear()
        self.annotationComboBox.addItems(self.annotationCodes)
        self.annotationComboBox.setCurrentIndex(index)
##        print "updated annotations"
        if not self.annotationCodes:
            self.error(0, "No downloaded annotations!!\nClick the update button and update annotationa for at least one organism!")
        else:
            self.error(0)

    def SetGenesComboBox(self):
        self.candidateGeneAttrs = self.clusterDataset.domain.variables + self.clusterDataset.domain.getmetas().values()
        self.candidateGeneAttrs = filter(lambda v: v.varType==orange.VarTypes.String or v.varType==orange.VarTypes.Other or v.varType==orange.VarTypes.Discrete, self.candidateGeneAttrs)
        self.geneAttrIndexCombo.clear()
        self.geneAttrIndexCombo.addItems([a.name for a in  self.candidateGeneAttrs])

    def FindBestGeneAttrAndOrganism(self):
        if self.autoFindBestOrg:  
            organismGenes = dict([(o,set(go.getCachedGeneNames(o))) for o in self.annotationCodes])
        else:
            currCode = self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
            filename = p_join(dataDir, self.annotationFiles[currCode])
            try:
                f = tarfile.open(filename)
                info = [info for info in f.getmembers() if info.name.startswith("gene_names")].pop()
                geneNames = cPickle.loads(f.extractfile(info).read().replace("\r\n", "\n"))
            except Exception, ex:
                geneNames = cPickle.loads(open(p_join(filename, "gene_names.pickle")).read().replace("\r\n", "\n"))
            organismGenes = {currCode: set(geneNames)}
        candidateGeneAttrs = self.clusterDataset.domain.attributes + self.clusterDataset.domain.getmetas().values()
        candidateGeneAttrs = filter(lambda v: v.varType==orange.VarTypes.String or v.varType==orange.VarTypes.Other or v.varType==orange.VarTypes.Discrete, candidateGeneAttrs)
        attrNames = [v.name for v in self.clusterDataset.domain.variables]
        cn = {}
        for attr in candidateGeneAttrs:
            vals = [str(e[attr]) for e in self.clusterDataset]
            if any("," in val for val in vals):
                vals = reduce(list.__add__, (val.split(",") for val in vals))
            for organism, s in organismGenes.items():
                l = filter(lambda a: a in s, vals)
                cn[(attr,organism)] = len(set(l))
        for organism, s in organismGenes.items():
            l = filter(lambda a: a in s, attrNames)
            cn[("_var_names_", organism)] = len(set(l))
            
        cn = cn.items()
        cn.sort(lambda a,b:-cmp(a[1],b[1]))
        ((bestAttr, organism), count) = cn[0]
##        print "match count:", count
        if bestAttr=="_var_names_" and count<=len(attrNames)/10.0 or \
           bestAttr!="_var_names_" and count<=len(self.clusterDataset)/10.0:
            return
        
        self.annotationIndex = self.annotationCodes.index(organism)
        if bestAttr=="_var_names_":
            self.useAttrNames = True
##            self.geneAttrIndexCombo.setDisabled(True)
            self.geneAttrIndex = 0
        else:
            self.useAttrNames = False
##            self.geneAttrIndexCombo.setDisabled(False)
            self.geneAttrIndex = candidateGeneAttrs.index(bestAttr)
    
    def SetClusterDataset(self, data=None):
        self.closeContext()
        self.clusterDataset = data
        self.infoLabel.setText("\n")
        if data:
            self.SetGenesComboBox()
##            index = self.annotationIndex = -42
            self.openContext("", data)
            if not self.ontology:
                self.LoadOntology()
            if not self.annotations or self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]!= self.loadedAnnotationCode:
                self.LoadAnnotation()

##            if self.annotationIndex == -42:
##                self.FindBestGeneAttrAndOrganism()
##            else:
##                self.annotationIndex = index
            
            self.FilterUnknownGenes()
            graph = self.Enrichment()
            self.SetGraph(graph)
        else:
            self.infoLabel.setText("No data on input\n")
            self.openContext("", None)
            self.ClearGraph()
            self.send("Selected Examples", None)
            self.send("Unselected Examples", None)
            self.send("Example With Unknown Genes", None)

    def SetReferenceDataset(self, data=None):
        self.referenceDataset=data
        self.referenceRadioBox.buttons[1].setDisabled(not bool(data))
        if self.clusterDataset and self.useReferenceDataset:
            graph = self.Enrichment()
            self.SetGraph(graph)

    def FilterUnknownGenes(self):
        if not self.useAttrNames:
            geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]
            examples = []
            for ex in self.clusterDataset:
##                if not any(n in go.loadedAnnotation.aliasMapper for n in str(ex[geneAttr]).split(",")):
                if not any(n.strip() in self.annotations.aliasMapper or n.strip() in self.annotations.additionalAliases for n in str(ex[geneAttr]).split(",")):
                    examples.append(ex)
##                if str(ex[geneAttr]) not in go.loadedAnnotation.aliasMapper:
##                    examples.append(ex)
            self.send("Example With Unknown Genes", examples and orange.ExampleTable(examples) or None)
        else:
            self.send("Example With Unknown Genes", None)

    def LoadOntology(self):
        try:
            self.progressBarInit()
            self.ontology = obiGO.Ontology.Load(progressCallback=self.progressBarSet)
            self.progressBarFinished()
        except IOError, er:
            response = QMessageBox.warning(self, "GOEnrichmentAnalysis", "Unable to load the ontology.\nClik OK to download it?", "OK", "Cancel", "", 0, 1)
            if response==0:
                self.UpdateGOAndAnnotation(tags = ["ontology", "GO", "essential"])
                self.ontology = obiGO.Ontology.Load(progressCallback=self.progressBarSet)
                self.progressBarFinished()
            else:
                raise
        
    def LoadAnnotation(self):
        if not self.annotationCodes:
            response = QMessageBox.warning(self, "GOEnrichmentAnalysis", "No annotations found.\nClick OK to download it", "OK", "Cancel", "", 0, 1)
            if response==0:
                self.UpdateGOAndAnnotation(tags=["annotation", "go"])
        if self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]!= self.loadedAnnotationCode:
            self.progressBarInit()
            try:
                self.annotations = None
##                filename = p_join(dataDir, self.annotationFiles[self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]])
##                self.annotations = obiGO.Annotations(filename, ontology = self.ontology, progressCallback=self.progressBarSet)
                self.annotations = obiGO.Annotations(self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)], ontology = self.ontology, progressCallback=self.progressBarSet)
            except IOError, er:
                raise
            self.progressBarFinished()
##            count = dict([(etype, 0) for etype in go.evidenceTypesOrdered])
##            geneSets = dict([(etype, set()) for etype in go.evidenceTypesOrdered])
            count = defaultdict(int)
            geneSets = defaultdict(set)
##            for anno in go.loadedAnnotation.annotationList:
            for anno in self.annotations.annotations:
                count[anno.evidence]+=1
                geneSets[anno.evidence].add(anno.geneName)
            for etype in obiGO.evidenceTypesOrdered:
                self.evidenceCheckBoxDict[etype].setEnabled(bool(count[etype]))
                self.evidenceCheckBoxDict[etype].setText(etype+": %i annots(%i genes)" % (count[etype], len(geneSets[etype])))
            self.loadedAnnotationCode=self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
            if self.loadedAnnotationCode in GeneMatchMk2.dbOrgMap:
                self.keggOrg = obiKEGG.KEGGOrganism(GeneMatchMk2.dbOrgMap[self.loadedAnnotationCode], update=False)
                self.keggOrg.api.download_progress_callback = self.progressBarSet
            else:
                self.keggOrg = None

    def UpdateGOAliases(self, genes):
##        genes = [gene for gene in genes if gene not in go.loadedAnnotation.aliasMapper]
        genes = [gene for gene in genes if gene not in self.annotations.aliasMapper]
        if not self.keggOrg or not os.path.isfile(os.path.join(self.keggOrg.local_database_path,"genes//organisms//"+self.keggOrg.org+"//_genes.pickle")):
            print "Gene translation failed"
            return
##        old = dict(go.loadedAnnotation.__annotation.aliasMapper)
        dbNames = set([anno.DB for anno in go.loadedAnnotation.annotationList])
        dbNames = [GeneMatchMk2.dbNameMap[db] for db in dbNames if db in GeneMatchMk2.dbNameMap]
        org = GeneMatchMk2.dbOrgMap[self.loadedAnnotationCode]
        try:
            self.progressBarInit()
            unique, c, u = self.keggOrg.get_unique_gene_ids(genes)
            self.progressBarFinished()
            for k, gene in unique.items():
                links = self.keggOrg.api._genes[org][k].get_db_links()
                for db in dbNames:
                    if db in links and len(set([link for link in links[db] if link in go.loadedAnnotation.aliasMapper]))==1:
                        go.loadedAnnotation.aliasMapper[gene] = go.loadedAnnotation.aliasMapper[set([link for link in links[db] if link in go.loadedAnnotation.aliasMapper]).pop()]
                        break
                altNames = self.keggOrg.api._genes[org][k].get_alt_names()
                altNames = [name for name in altNames if name in go.loadedAnnotation.aliasMapper]
                if len(set([go.loadedAnnotation.aliasMapper[name] for name in altNames]))==1:
                    go.loadedAnnotation.aliasMapper[gene] = go.loadedAnnotation.aliasMapper[altNames[0]]
        except IOError, ex:
            print ex
        return
    
    def Enrichment(self):
        if not self.annotations.ontology:
            self.annotations.ontology = self.ontology
        if self.useAttrNames:
            clusterGenes = [v.name for v in self.clusterDataset.domain.variables]
            self.information(0)
        else:
            geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]
            clusterGenes = [str(ex[geneAttr]) for ex in self.clusterDataset if not ex[geneAttr].isSpecial()]
            if any("," in gene for gene in clusterGenes):
                self.information(0, "Separators detected in cluster gene names. Assuming multiple genes per example.")
                clusterGenes = reduce(list.__add__, (genes.split(",") for genes in clusterGenes))
            else:
                self.information(0)
##        self.UpdateGOAliases(clusterGenes)
##        self.geneInfoLabel.setText("%i genes on input" % len(clusterGenes))
        genesCount = len(clusterGenes)
##        self.clusterGenes = clusterGenes = filter(lambda g: g in go.loadedAnnotation.aliasMapper, clusterGenes)
        self.clusterGenes = clusterGenes = filter(lambda g: g in self.annotations.aliasMapper or g in self.annotations.additionalAliases, clusterGenes)
        self.infoLabel.setText("%i genes on input\n%i (%.1f%%) gene names matched" % (genesCount, len(clusterGenes), 100.0*len(clusterGenes)/genesCount if genesCount else 0.0))
##        print len(self.clusterGenes), self.clusterGenes[:5]
        referenceGenes = None
        if self.useReferenceDataset:
            try:
                if self.useAttrNames:
                    referenceGenes = [v.name for v in self.referenceDataset.domain.variables]
                    self.information(1)
                else:
                    referenceGenes = [str(ex[geneAttr]) for ex in self.referenceDataset if not ex[geneAttr].isSpecial()]
                    if any("," in gene for gene in clusterGenes):
                        self.information(1, "Separators detected in reference gene names. Assuming multiple genes per example.")
                        referenceGenes = reduce(list.__add__, (genes.split(",") for genes in referenceGenes))
                    else:
                        self.information(1)
##                self.UpdateGOAliases(referenceGenes)
##                referenceGenes = filter(lambda g: g in go.loadedAnnotation.aliasMapper, referenceGenes)
                referenceGenes = filter(lambda g: g in self.annotations.aliasMapper or g in self.annotations.additionalAliases, referenceGenes)
                self.information(2)
            except Exception, er:
                self.information(2, str(er)+" Using entire genome for reference")
                referenceGenes = self.annotations.geneNames
                self.useReferenceDataset = 0
        else:
            self.information(2)
##            referenceGenes = go.loadedAnnotation.geneNames
            referenceGenes = self.annotations.geneNames
        self.referenceGenes = referenceGenes
        evidences = []
        for etype in obiGO.evidenceTypesOrdered:
            if getattr(self, "useEvidence"+etype):
                evidences.append(etype)
        aspect = ["P", "C", "F"][self.aspectIndex]
        self.progressBarInit()
        if clusterGenes:
##            print clusterGenes[:5], referenceGenes[:5], evidences, aspect
##            self.terms = terms = go.GOTermFinder(clusterGenes, referenceGenes, evidences, aspect=aspect, progressCallback=self.progressBarSet)
            self.terms = terms = self.annotations.GetEnrichedTerms(clusterGenes, referenceGenes, evidences, aspect=aspect,
                                                                   prob=self.probFunctions[self.probFunc], progressCallback=self.progressBarSet)
            if self.useFDR:
                terms = sorted(terms.items(), key=lambda (_1, (_2, p, _3)): p)
                p_vals = obiProb.FDR([p for _, (_, p, _) in terms])
                self.terms = terms = dict([(id, (genes, p, ref)) for p, (id, (genes, _, ref)) in zip(p_vals, terms)])
##            go.loadedAnnotation.__annotation.aliasMapper = old
        else:
            self.terms = terms = {}
        if not self.terms:
            self.warning(0, "No terms found")
        else:
            self.warning(0)
        self.progressBarFinished()
        self.treeStructDict = {}
        ids = self.terms.keys()
        for term in self.terms:
##            self.treeStructDict[term] = TreeNode(self.terms[term], filter(lambda t:term in go.loadedGO.termDict[t].parents, ids))
            parents = lambda t: [term for typeId, term in  self.ontology[t].related]
            self.treeStructDict[term] = TreeNode(self.terms[term], [id for id in ids if term in parents(id)])
##            if not go.loadedGO.termDict[term].parents:
            if not self.ontology[term].related and not getattr(self.ontology[term], "is_obsolete", False):
                self.treeStructRootKey = term
        return terms
        
    def FilterGraph(self, graph):
        if self.filterByPValue:
            graph = obiGO.filterByPValue(graph, self.maxPValue)
        if self.filterByNumOfInstances:
            graph = dict(filter(lambda (id,(genes, p, rc)):len(genes)>=self.minNumOfInstances, graph.items()))
        return graph

    def FilterAndDisplayGraph(self):
        if self.clusterDataset:
            self.graph = self.FilterGraph(self.originalGraph)
            self.ClearGraph()
            self.DisplayGraph()

    def SetGraph(self, graph=None):
        self.originalGraph = graph
        if graph:
            self.FilterAndDisplayGraph()
        else:
            self.graph = {}
            self.ClearGraph()

    def ClearGraph(self):
        self.listView.clear()
        self.listViewItems=[]
        self.sigTerms.clear()
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
                displayNode = GOTreeWidgetItem(self.ontology[term], self.graph[term], len(self.clusterGenes), len(self.referenceGenes), maxFoldEnrichment, parentDisplayNode)
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
##        self.sigTermsTable.setRowCount(len(terms))
        self.sigTerms.clear()
        for i, (id, (genes, p_value, refCount)) in enumerate(terms):
##            text = [go.loadedGO.termDict[id].name, str(len(genes)), str(refCount), "%.4f" % p_value, " ,".join(genes), "%.2f" % enrichment((genes, p_value, refCount))]
            text = [self.ontology[id].name, str(len(genes)), str(refCount), "%.4f" % p_value, " ,".join(genes), "%.2f" % enrichment((genes, p_value, refCount))]
            item = GOTreeWidgetItem(self.ontology[id], (genes, p_value, refCount), len(self.clusterGenes),
                                    len(self.referenceGenes), maxFoldEnrichment, self.sigTerms)
            item.goId = id
##            for j,t in enumerate(text):
##                self.sigTermsTable.setItem(i, j, QTableWidgetItem(t))
                
        self.listView.expandAll()
        self.listView.resizeColumnToContents(0)
##        print [item.sizeHint(0).width() for item in self.listViewItems]
        width = min(self.listView.columnWidth(0), 500)
        self.listView.setColumnWidth(0, width)
        self.sigTerms.setColumnWidth(0, width)
        
    def ViewSelectionChanged(self):
        if self.selectionChanging:
            return
        
        self.selectionChanging = 1
        self.selectedTerms = []
        #selected = filter(lambda lvi: lvi.isSelected(), self.listViewItems)
        selected = self.listView.selectedItems()
        self.selectedTerms = list(set([lvi.term for lvi in selected]))
        self.ExampleSelection()
        self.selectionChanging = 0
        
        
    def TableSelectionChanged(self):
        if self.selectionChanging:
            return
        
        self.selectionChanging = 1
        self.selectedTerms = []
        selectedIds = set([self.sigTerms.itemFromIndex(index).goId for index in self.sigTerms.selectedIndexes()])
        
        for i in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(i)
            selected = item.goId in selectedIds
##            term = self.sigTableTermsSorted[row]
            term = item.goId
            
            if selected:
                self.selectedTerms.append(term)
                
            for lvi in self.termListViewItemDict[term]:
                try:
                    lvi.setSelected(selected)
                    #self.listView.repaintItem(lvi)
                    if selected: lvi.setExpanded(True)
                except RuntimeError:    ##Underlying C/C++ object deleted (why??)
##                    print "error 11"
                    pass
                
        #self.listView.triggerUpdate()
        self.ExampleSelection()
        self.selectionChanging = 0
            
    
    def ExampleSelection(self):
        selectedExamples = []
        unselectedExamples = []
        selectedGenes = []

        #change by Marko. don't do anything if there is no3 dataset 
        if not self.clusterDataset:
            return

        if self.selectionDirectAnnotation:
##            s = filter(lambda anno: anno.GOId in self.selectedTerms, go.loadedAnnotation.annotationList)
            s = filter(lambda anno: anno.GOId in self.selectedTerms, self.annotations.annotations)
            selectedGenes = set([anno.geneName for anno in s])
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
            newDomain = orange.Domain(vars, self.clusterDataset.domain.classVar)
            self.send("Selected Examples", orange.ExampleTable(newDomain, self.clusterDataset))
            self.send("Unselected Examples", None)
        else:
            geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]
            newClass = orange.EnumVariable("GO Term", values=list(self.selectedTerms))
            newDomain = orange.Domain(self.clusterDataset.domain.variables, newClass)
            for ex in self.clusterDataset:
                if not ex[geneAttr].isSpecial() and any(gene in selectedGenes for gene in str(ex[geneAttr]).split(",")):
                    if self.selectionDisjoint and self.selectionAddTermAsClass:
                        c = filter(lambda term: any(gene in self.graph[term][0] for gene in str(ex[geneAttr]).split(",")) , self.selectedTerms)[0]
                        ex =  orange.Example(newDomain, ex)
                        ex.setclass(newClass(c))
                    selectedExamples.append(ex)
                else:
                    unselectedExamples.append(ex)
            self.send("Selected Examples", selectedExamples and orange.ExampleTable(selectedExamples) or None)
            self.send("Unselected Examples", unselectedExamples and orange.ExampleTable(unselectedExamples) or None)

    def ShowInfo(self):
        dialog = QDialog(self)
        dialog.setModal(False)
        dialog.setLayout(QVBoxLayout())
        label = QLabel(dialog)
        label.setText("Ontology:\n"+self.ontology.header if self.ontology else "Ontology not loaded!")
        dialog.layout().addWidget(label)

        label = QLabel(dialog)
        label.setText("Annotations:\n"+self.annotations.header.replace("!", "") if self.annotations else "Annotations not loaded!")
        dialog.layout().addWidget(label)
        dialog.show()

    def sendReport(self):
        self.reportSettings("Settings", [("Organism", self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes) - 1)]),
                                         ("Significance test", ("Binomial" if self.probFunc == 0 else "Hypergeometric") + (" with FDR" if self.useFDR else ""))])
        self.reportSettings("Filter", ([("Min cluster size", self.minNumOfInstances)] if self.filterByNumOfInstances else []) + \
                                      ([("Max p-value", self.maxPValue)] if self.filterByPValue else []))
##        self.reportSubsection("Enriched terms:")
##        text = [["Term:", "list:", "reference:", "p-value:", "enrichment:"]] + \
##               [[self.sigTerms.topLevelItem(index).text(i) for i in (range(4) + [5])] for index in range(self.sigTerms.topLevelItemCount())]
##        widths = reduce(lambda widths, line: [max(w, len(t)) for w, t in zip(widths, line)], text, [0]*5)
##        table = " ".join(["%-" + str(w) + "s" for w in widths]) % tuple(text[0])
##        fmt = " ".join(["%" + str(w) + "s" for w in widths])
##        table += "\n" + "\n".join([fmt % tuple(line) for line in text[1:]])
##        table = "<PRE>%s\n</PRE>" % table
        def _itemText(item, level):
            text = '<tr><td>%s%s</td>' % ('&nbsp;' * level * 2 + '-', item.text(0))
            text += ''.join('<td>%s</td>' % item.text(i) for i in range(1, 4) + [5]) + '</tr>\n'
            for i in range(item.childCount()):
                text += _itemText(item.child(i), level + 1)
            return text 
        tableText = '<table border= "1"><caption>Significant terms</caption><tr>' + ''.join('<th>%s</th>' % s for s in ["Term:", "List:", "Reference:", "P-value:", "Enrichment:"]) + '</td>'
        treeText = '<table border= "1"><caption>Significant terms in ontology structure</caption><tr>' + ''.join('<th>%s</th>' % s for s in ["Term:", "List:", "Reference:", "P-value:", "Enrichment:"]) + '</td>'
        
        for index in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(index)
            tableText += _itemText(item, 0) 
##            text += '<tr>' + ''.join('<td>%s</td>' % item.text(i) for i in (range(4) + [5])) + '</tr>'
        tableText += '</table>' 
        
        for index in range(self.listView.topLevelItemCount()):
            item = self.listView.topLevelItem(index)
            treeText += _itemText(item, 0)
        
        self.reportRaw(tableText)
        self.reportRaw(treeText)
        
        
        for index in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(index) 

class GOTreeWidgetItem(QTreeWidgetItem):
    def __init__(self, term, enrichmentResult, nClusterGenes, nRefGenes, maxFoldEnrichment, parent):
        QTreeWidgetItem.__init__(self, parent)
        self.term = term
        self.enrichmentResult = enrichmentResult
        self.nClusterGenes = nClusterGenes
        self.nRefGenes = nRefGenes
        self.maxFoldEnrichment = maxFoldEnrichment
        self.enrichment = enrichment = lambda t:float(len(t[0])) / t[2] * (float(nRefGenes)/nClusterGenes)
        self.setText(0, term.name)
        fmt = "%" + str(-int(math.log(nClusterGenes))) + "i (%.2f%%)"
##        self.setText(1, "%i (%.2f%%)" % (len(enrichmentResult[0]), 100.0*len(self.enrichmentResult[0])/nClusterGenes))
        self.setText(1, fmt % (len(enrichmentResult[0]), 100.0*len(self.enrichmentResult[0])/nClusterGenes))
        fmt = "%" + str(-int(math.log(nRefGenes))) + "i (%.2f%%)"
##        self.setText(2, "%i (%.2f%%)" % (enrichmentResult[2], 100.0*enrichmentResult[2]/nRefGenes))
        self.setText(2, fmt % (enrichmentResult[2], 100.0*enrichmentResult[2]/nRefGenes))
        self.setText(3, "%.4f" % enrichmentResult[1])
        self.setText(4, ", ".join(enrichmentResult[0]))
        self.setText(5, "%.2f" % (enrichment(enrichmentResult))) #(float(len(self.graph[term][0]))/self.graph[term][2]))
        self.setToolTip(0, "<p>" + term.__repr__()[6:].strip().replace("\n", "<br>"))
        self.sortByData = [term.name, len(self.enrichmentResult[0]), enrichmentResult[2], enrichmentResult[1], ", ".join(enrichmentResult[0]), enrichment(enrichmentResult)]

    def data(self, col, role):
        if role == Qt.UserRole:
            return QVariant(self.enrichment(self.enrichmentResult) / self.maxFoldEnrichment)
        else:
            return QTreeWidgetItem.data(self, col, role)

    def __lt__(self, other):
        col = self.treeWidget().sortColumn()
        return self.sortByData[col] < other.sortByData[col]
    
class EnrichmentColumnItemDelegate(QItemDelegate):
    def paint(self, painter, option, index):
        self.drawBackground(painter, option, index)
        value, ok = index.data(Qt.UserRole).toDouble()
        if ok:
            painter.save()
            painter.setBrush(QBrush(Qt.white, Qt.SolidPattern))
            painter.drawRect(option.rect)
            painter.setBrush(QBrush(Qt.blue, Qt.SolidPattern))
            painter.drawRect(option.rect.x(), option.rect.y(), value*(option.rect.width()-1), option.rect.height()-1)
            painter.restore()
        else:
            QItemDelegate.paint(self, painter, option, index)
        
        
if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    w=OWGOEnrichmentAnalysis()
    data = orange.ExampleTable("../../orange/doc/datasets/brown-selected.tab")
    w.show()
    w.SetClusterDataset(data)
    app.exec_()
    w.saveSettings()
        
