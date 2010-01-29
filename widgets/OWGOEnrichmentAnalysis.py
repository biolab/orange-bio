"""
<name>GO Enrichment Analysis</name>
<description>Enrichment analysis for Gene Ontology terms.</description>
<contact>Ales Erjavec</contact>
<icon>icons/GOTermFinder.png</icon>
<priority>2020</priority>
"""

from __future__ import with_statement

import obiGO
import obiProb
import obiTaxonomy
import obiGene
import sys, os, tarfile, math
import gc
import OWGUI
import orngServerFiles

from os.path import join as p_join
from OWWidget import *
from collections import defaultdict
from functools import partial

dataDir = orngServerFiles.localpath("GO")

def listAvailable():
    files = orngServerFiles.listfiles("GO")
    ret = {}
    for file in files:
        tags = orngServerFiles.info("GO", file)["tags"]
        td = dict([tuple(tag.split(":")) for tag in tags if tag.startswith("#") and ":" in tag])
        if "association" in file.lower():
            ret[td.get("#organism", file)] = file
    orgMap = {"352472":"44689"}
    essential = ["gene_association.%s.tar.gz" % obiGO.from_taxid(id) for id in obiTaxonomy.essential_taxids() if obiGO.from_taxid(id)]
    essentialNames = [obiTaxonomy.name(id) for id in obiTaxonomy.essential_taxids() if obiGO.from_taxid(id)]
    ret.update(zip(essentialNames, essential))
    return ret

class _disablegc(object):
    def __enter__(self):
        gc.disable()
    def __exit__(self, *args):
        gc.enable()

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
    settingsList=["annotationIndex", "useReferenceDataset", "aspectIndex", "geneAttrIndex", "geneMatcherSettings",
                    "filterByNumOfInstances", "minNumOfInstances", "filterByPValue", "maxPValue", "selectionDirectAnnotation", "selectionDisjoint", "selectionType",
                    "selectionAddTermAsClass", "useAttrNames", "probFunc", "useFDR"]
    contextHandlers = {"": DomainContextHandler("", ["geneAttrIndex", "useAttrNames", "annotationIndex", "geneMatcherSettings"], matchValues=1)}
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
        self.geneMatcherSettings = [True, False, False, False]
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
#        self.progressBarInit()
#        with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
#            self.annotationFiles = listAvailable()
#        self.progressBarFinished()
#        self.annotationCodes = sorted(self.annotationFiles.keys())
#        if not self.annotationCodes:
#            self.error(0, "No downloaded annotations!!\nUse the Update Genomics Databases widget and download annotations for at least one organism!")
#        else:
#            self.error(0)
        self.annotationCodes = []
        
        #############
        ##GUI
        #############
        self.tabs = OWGUI.tabWidget(self.controlArea)
        ##Input tab
        self.inputTab = OWGUI.createTabPage(self.tabs, "Input")
        box = OWGUI.widgetBox(self.inputTab, "Info")
        self.infoLabel = OWGUI.widgetLabel(box, "No data on input\n")
        OWGUI.button(box, self, "Ontology/Annotation Info", callback=self.ShowInfo, tooltip="Show information on loaded ontology and annotations", debuggingEnabled=0)
        box = OWGUI.widgetBox(self.inputTab, "Organism", addSpace=True)
        self.annotationComboBox = OWGUI.comboBox(box, self, "annotationIndex", items = self.annotationCodes, callback=self.Update, tooltip="Select organism", debuggingEnabled=0)
        
        self.signalManager.setFreeze(1) ## freeze until annotation combo box is updateded with available annotations.
        QTimer.singleShot(0, self.UpdateOrganismComboBox)
        
        self.geneAttrIndexCombo = OWGUI.comboBox(self.inputTab, self, "geneAttrIndex", box="Gene names", callback=self.Update, tooltip="Use this attribute to extract gene names from input data")
        OWGUI.checkBox(self.geneAttrIndexCombo.box, self, "useAttrNames", "Use data attributes names", disables=[(-1, self.geneAttrIndexCombo)], callback=self.Update, tooltip="Use attribute names for gene names")
        OWGUI.button(self.geneAttrIndexCombo.box, self, "Gene matcher settings", callback=self.UpdateGeneMatcher, tooltip="Open gene matching settings dialog", debuggingEnabled=0)
        
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
            self.evidenceCheckBoxDict[etype] = OWGUI.checkBox(box, self, "useEvidence"+etype, etype, callback=self.Update, tooltip=obiGO.evidenceTypes[etype])
        
        ##Select tab
        self.selectTab = OWGUI.createTabPage(self.tabs, "Select")
        #box = OWGUI.widgetBox(self.selectTab, "Annotated genes", addSpace=True)
        box = OWGUI.radioButtonsInBox(self.selectTab, self, "selectionDirectAnnotation", ["Directly or Indirectly", "Directly"], box="Annotated genes", callback=self.ExampleSelection)
        box = OWGUI.widgetBox(self.selectTab, "Output", addSpace=True)
##        OWGUI.checkBox(box, self, "selectionDisjoint", "Disjoint/Inclusive", callback=self.ExampleSelection)
        OWGUI.radioButtonsInBox(box, self, "selectionDisjoint", btnLabels=["All selected genes", "Term-specific genes", "Common term genes"], tooltips=["Outputs genes annotated to all selected GO terms", "Outputs genes that appear in only one of selected GO terms", "Outputs genes common to all selected GO terms"], callback=self.ExampleSelection)
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
        
        self.resize(1000, 800)

        self.clusterDataset = None
        self.referenceDataset = None
        self.ontology = None
        self.annotations = None
        self.probFunctions = [obiProb.Binomial(), obiProb.Hypergeometric()]
        self.selectedTerms = []
        
    def UpdateOrganismComboBox(self):
        try:
            if self.annotationCodes and len(self.annotationCodes) > self.annotationIndex:
                currAnnotationCode = self.annotationCodes[self.annotationIndex]
            else:
                currAnnotationCode = None
            self.progressBarInit()
            with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
                self.annotationFiles = listAvailable()
            self.progressBarFinished()
            self.annotationCodes = sorted(self.annotationFiles.keys())
    #        if not self.annotationCodes:
    #            self.error(0, "No downloaded annotations!!\nUse the Update Genomics Databases widget and download annotations for at least one organism!")
    #        else:
    #            self.error(0)
            self.annotationComboBox.clear()
            self.annotationComboBox.addItems(self.annotationCodes)
#            self.annotationIndex = self.annotationCodes.index(currAnnotationCode) if currAnnotationCode in self.annotationCodes else 0
            self.annotationComboBox.setCurrentIndex(self.annotationIndex)
#            print "update", self.annotationIndex, currAnnotationCode
        finally:
            self.signalManager.setFreeze(0)
            
    def UpdateGeneMatcher(self):
        dialog = GeneMatcherDialog(self, defaults=self.geneMatcherSettings, modal=True)
        if dialog.exec_():
            self.geneMatcherSettings = [getattr(dialog, item[0]) for item in dialog.items]
            if self.annotations:
                self.SetGeneMatcher()
                if self.clusterDataset:
                    self.Update()
                
    def Update(self):
        if self.clusterDataset:
            pb = OWGUI.ProgressBar(self, 100)
            self.Load(pb=pb)
            self.FilterUnknownGenes()
            graph = self.Enrichment(pb=pb)
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
            self.openContext("", data)
#            if not self.ontology:
#                self.LoadOntology()
#            if not self.annotations or self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]!= self.loadedAnnotationCode:
#                self.LoadAnnotation()
#            
#            self.Load()
#            self.FilterUnknownGenes()
#            graph = self.Enrichment()
#            self.SetGraph(graph)
            self.Update()
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
        self.referenceRadioBox.buttons[1].setText("Reference set")
        if self.clusterDataset and self.useReferenceDataset:
            self.useReferenceDataset = 0 if not data else 1
            graph = self.Enrichment()
            self.SetGraph(graph)
        elif self.clusterDataset:
            self.UpdateReferenceSetButton()
            
        
    def UpdateReferenceSetButton(self):
        allgenes, refgenes = None, None
        if self.referenceDataset:
            try:
                allgenes = self.GenesFromExampleTable(self.referenceDataset)
            except Exception:
                allgenes = []
            refgenes, unknown = self.FilterAnnotatedGenes(allgenes)
        self.referenceRadioBox.buttons[1].setDisabled(not bool(allgenes))
        self.referenceRadioBox.buttons[1].setText("Reference set " + ("(%i genes, %i matched)" % (len(allgenes), len(refgenes)) if allgenes and refgenes else ""))

    def GenesFromExampleTable(self, data):
        if self.useAttrNames:
            genes = [v.name for v in data.domain.variables]
        else:
            attr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs) - 1)]
            genes = [str(ex[attr]) for ex in data if not ex[attr].isSpecial()]
            if any("," in gene for gene in genes):
                self.information(0, "Separators detected in gene names. Assuming multiple genes per example.")
                genes = reduce(list.__add__, (genes.split(",") for genes in genes))
        return genes
        
    def FilterAnnotatedGenes(self, genes):
        matchedgenes = self.annotations.GetGeneNamesTranslator(genes).values()
        return matchedgenes, [gene for gene in genes if gene not in matchedgenes]
        
    def FilterUnknownGenes(self):
        if not self.useAttrNames:
            geneAttr = self.candidateGeneAttrs[min(self.geneAttrIndex, len(self.candidateGeneAttrs)-1)]
            examples = []
            for ex in self.clusterDataset:
                if not any(self.annotations.genematcher.match(n.strip()) for n in str(ex[geneAttr]).split(",")):
                    examples.append(ex)

            self.send("Example With Unknown Genes", examples and orange.ExampleTable(examples) or None)
        else:
            self.send("Example With Unknown Genes", None)

    def Load(self, pb=None):
        go_files, tax_files = orngServerFiles.listfiles("GO"), orngServerFiles.listfiles("Taxonomy")
        calls = []
        pb, finish = (OWGUI.ProgressBar(self, 0), True) if pb is None else (pb, False)
        count = 0
        if not tax_files:
            calls.append(("Taxonomy", "ncbi_taxnomy.tar.gz"))
            count += 1
        org = self.annotationCodes[min(self.annotationIndex, len(self.annotationCodes)-1)]
        if org != self.loadedAnnotationCode:
            count += 1
            if self.annotationFiles[org] not in go_files:
                calls.append(("GO", self.annotationFiles[org]))
                count += 1
                
        if "gene_ontology_edit.obo.tar.gz" not in go_files:
            calls.append(("GO", "gene_ontology_edit.obo.tar.gz"))
            count += 1
        if not self.ontology:
            count += 1
        pb.iter += count*100
#        self.progressBarInit()
        for i, args in enumerate(calls):
#            with orngServerFiles.DownloadProgress.setredirect(lambda value: self.progressBarSet(100.0 * i / count + value/count)):
#                print args
            orngServerFiles.localpath_download(*args, **dict(callback=pb.advance))
            
        i = len(calls)
        if not self.ontology:
            self.ontology = obiGO.Ontology(progressCallback=lambda value: pb.advance()) #self.progressBarSet(100.0 * i / count + value/count))
            i+=1
        if org != self.loadedAnnotationCode:
            code = self.annotationFiles[org].split(".")[-3]
            self.annotations = obiGO.Annotations(code, genematcher=obiGene.GMDirect(), progressCallback=lambda value: pb.advance())#self.progressBarSet(100.0 * i / count + value/count))
            i+=1
            self.loadedAnnotationCode = org
            count = defaultdict(int)
            geneSets = defaultdict(set)

            for anno in self.annotations.annotations:
                count[anno.evidence]+=1
                geneSets[anno.evidence].add(anno.geneName)
            for etype in obiGO.evidenceTypesOrdered:
                self.evidenceCheckBoxDict[etype].setEnabled(bool(count[etype]))
                self.evidenceCheckBoxDict[etype].setText(etype+": %i annots(%i genes)" % (count[etype], len(geneSets[etype])))
        if finish:
            pb.finish()
#        self.progressBarFinished()
            
    def SetGeneMatcher(self):
        if self.annotations:
            taxid = self.annotations.taxid
            matchers = []
            for matcher, use in zip([obiGene.GMGO, obiGene.GMKEGG, obiGene.GMNCBI, obiGene.GMAffy], self.geneMatcherSettings):
                if use:
                    try:
                        matchers.append(matcher(taxid))
                    except Exception, ex:
                        print ex
            matchers.reverse()
#            print matchers
            self.annotations.genematcher = obiGene.matcher(matchers)
#            self.progressBarInit()
#            with orngServerFiles.DownloadProgress.setredirect(self.progressBarSet):
            self.annotations.genematcher.set_targets(self.annotations.geneNames)
#            self.progressBarFinished()
            
    def Enrichment(self, pb=None):
        pb = OWGUI.ProgressBar(self, 100) if pb is None else pb
        if not self.annotations.ontology:
            self.annotations.ontology = self.ontology
            
        if isinstance(self.annotations.genematcher, obiGene.GMDirect):
            self.SetGeneMatcher()
            
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

        genesCount = len(clusterGenes)
        
        self.clusterGenes = clusterGenes = self.annotations.GetGeneNamesTranslator(clusterGenes).values()
        
#        self.clusterGenes = clusterGenes = filter(lambda g: g in self.annotations.aliasMapper or g in self.annotations.additionalAliases, clusterGenes)
        self.infoLabel.setText("%i genes on input\n%i (%.1f%%) gene names matched" % (genesCount, len(clusterGenes), 100.0*len(clusterGenes)/genesCount if genesCount else 0.0))
        
        referenceGenes = None
        if self.referenceDataset:
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

                refc = len(referenceGenes)
#                referenceGenes = filter(lambda g: g in self.annotations.aliasMapper or g in self.annotations.additionalAliases, referenceGenes)
                referenceGenes = self.annotations.GetGeneNamesTranslator(referenceGenes).values()
                self.referenceRadioBox.buttons[1].setText("Reference set (%i genes, %i matched)" % (refc, len(referenceGenes)))
                self.referenceRadioBox.buttons[1].setDisabled(False)
                self.information(2)
            except Exception, er:
                if not self.referenceDataset:
                    self.information(2, "Unable to extract gene names from reference dataset. Using entire genome for reference")
                else:
                    self.referenceRadioBox.buttons[1].setText("Reference set")
                    self.referenceRadioBox.buttons[1].setDisabled(True)
                referenceGenes = self.annotations.geneNames
                self.useReferenceDataset = 0
        else:
            self.useReferenceDataset = 0
        if not self.useReferenceDataset:
            self.information(2)
            self.information(1)
            referenceGenes = self.annotations.geneNames
        self.referenceGenes = referenceGenes
        evidences = []
        for etype in obiGO.evidenceTypesOrdered:
            if getattr(self, "useEvidence"+etype):
                evidences.append(etype)
        aspect = ["P", "C", "F"][self.aspectIndex]
#        self.progressBarInit()
        if clusterGenes:
            self.terms = terms = self.annotations.GetEnrichedTerms(clusterGenes, referenceGenes, evidences, aspect=aspect,
                                                                   prob=self.probFunctions[self.probFunc], progressCallback=lambda value:pb.advance() )#self.progressBarSet)
            if self.useFDR:
                terms = sorted(terms.items(), key=lambda (_1, (_2, p, _3)): p)
                p_vals = obiProb.FDR([p for _, (_, p, _) in terms])
                self.terms = terms = dict([(id, (genes, p, ref)) for p, (id, (genes, _, ref)) in zip(p_vals, terms)])
        else:
            self.terms = terms = {}
        if not self.terms:
            self.warning(0, "No terms found")
        else:
            self.warning(0)
#        self.progressBarFinished()
        pb.finish()
        self.treeStructDict = {}
        ids = self.terms.keys()
        for term in self.terms:
            parents = lambda t: [term for typeId, term in  self.ontology[t].related]
            self.treeStructDict[term] = TreeNode(self.terms[term], [id for id in ids if term in parents(id)])
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
                displayNode.goId = term
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
        for i in range(4):
            self.listView.resizeColumnToContents(i)
            self.sigTerms.resizeColumnToContents(i)
        self.sigTerms.resizeColumnToContents(5)
##        print [item.sizeHint(0).width() for item in self.listViewItems]
        width = min(self.listView.columnWidth(0), 350)
        self.listView.setColumnWidth(0, width)
        self.sigTerms.setColumnWidth(0, width)
        
    def ViewSelectionChanged(self):
        if self.selectionChanging:
            return
        
        self.selectionChanging = 1
        self.selectedTerms = []
        #selected = filter(lambda lvi: lvi.isSelected(), self.listViewItems)
        selected = self.listView.selectedItems()
        self.selectedTerms = list(set([lvi.term.id for lvi in selected]))
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
        
        selectedGenes = reduce(set.union, [v[0] for id, v in self.graph.items() if id in self.selectedTerms], set())
        evidences = []
        for etype in obiGO.evidenceTypesOrdered:
            if getattr(self, "useEvidence"+etype):
                evidences.append(etype)
        allTerms = self.annotations.GetAnnotatedTerms(selectedGenes, directAnnotationOnly=self.selectionDirectAnnotation, evidenceCodes=evidences)
            
        if self.selectionDisjoint:
            ##count = dict([(g, 0) for g in self.clusterGenes])
            count = defaultdict(int)
            for term in self.selectedTerms:
                ##for g in self.graph[term][0]:
                for g in allTerms.get(term, []):
                    count[g]+=1
            ccount = 1 if self.selectionDisjoint==1 else len(self.selectedTerms)
            selectedGenes = [gene for gene, c in count.items() if c==ccount and gene in selectedGenes]
        else:
            selectedGenes = reduce(set.union, [allTerms.get(term, []) for term in self.selectedTerms], set())

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

        def treeDepth(item):
            return 1 + max([treeDepth(item.child(i)) for i in range(item.childCount())] +[0])
        
        def printTree(item, level, treeDepth):
            text = '<tr>' + '<td width=16px></td>' * level
            text += '<td colspan="%i">%s: %s</td>' % (treeDepth - level, item.term.id, item.term.name)
            text += ''.join('<td>%s</td>' % item.text(i) for i in range(1, 4) + [5]) + '</tr>\n'
            for i in range(item.childCount()):
                text += printTree(item.child(i), level + 1, treeDepth)
            return text
        
        treeDepth = max([treeDepth(self.listView.topLevelItem(i)) for i in range(self.listView.topLevelItemCount())] + [0])
        
        tableText = '<table>\n<tr>' + ''.join('<th>%s</th>' % s for s in ["Term:", "List:", "Reference:", "P-value:", "Enrichment:"]) + '</tr>'
        
        treeText = '<table>\n' +  '<th colspan="%i">%s</th>' % (treeDepth, "Term:") 
        treeText += ''.join('<th>%s</th>' % s for s in ["List:", "Reference:", "P-value:", "Enrichment:"]) + '</tr>'
        
        for index in range(self.sigTerms.topLevelItemCount()):
            item = self.sigTerms.topLevelItem(index)
            tableText += printTree(item, 0, 1) 
###            text += '<tr>' + ''.join('<td>%s</td>' % item.text(i) for i in (range(4) + [5])) + '</tr>'
        tableText += '</table>' 
        
        for index in range(self.listView.topLevelItemCount()):
            item = self.listView.topLevelItem(index)
            treeText += printTree(item, 0, treeDepth)
        
        self.reportSection("Enriched Terms")
        self.reportRaw(tableText)
        
        self.reportSection("Enriched Terms in the Ontology Tree")
        self.reportRaw(treeText)

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
        
        
class GeneMatcherDialog(OWWidget):
    items = [("useGO", "Use gene names from Gene Ontology annotations"),
             ("useKEGG", "Use gene names from KEGG Genes database"),
             ("useNCBI", "Use gene names from NCBI Gene info database"),
             ("useAffy", "Use Affymetrix platform reference ids")]
    settingsList = [item[0] for item in items]
    def __init__(self, parent=None, defaults=[True, False, False, False], enabled=[False, True, True, True], **kwargs):
        OWWidget.__init__(self, parent, **kwargs)
        for item, default in zip(self.items, defaults):
            setattr(self, item[0], default)
            
        self.loadSettings()
        for item, enable in zip(self.items, enabled):
            cb = OWGUI.checkBox(self, self, *item)
            cb.setEnabled(enable)
            
        box = OWGUI.widgetBox(self, orientation="horizontal")
        OWGUI.button(box, self, "OK", callback=self.accept)
        OWGUI.button(box, self, "Cancel", callback=self.reject)
        
        
if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    w=OWGOEnrichmentAnalysis()
    data = orange.ExampleTable("../../orange/doc/datasets/brown-selected.tab")
    w.show()
    w.SetClusterDataset(data)
    app.exec_()
    w.saveSettings()
        
