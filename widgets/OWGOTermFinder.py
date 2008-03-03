"""
<name>GO Term Finder</name>
<description>GO Term Finder</description>
<contact>Tomaz Curk</contact>
<icon>icons/GOTermFinder.png</icon>
<priority>100</priority>
"""

import orange, math, glob
#import GOlib ## function needed to handle the GO and annotation
import OWGUI
from qt import *
from qtcanvas import *
from OWWidget import *
from OWOptions import *
from qttable import *
from qwt import *
from sets import Set

from OWDataFiles import DataFiles, ExampleSelection

#import pywin.debugger

try:
    import go
except:
    QMessageBox.warning( None, "Missing GOLib", "GOLib libray for handling GO ontology not found.\nYou can get it at ..." , QMessageBox.Ok)

DEBUG = 0

class OWGOTermFinder(OWWidget):	
    settingsList = ["AnnotationFileName", "RecentAnnotations", "ReferenceType", "RecentGOaspects",
                    "FilterNumEnabled", "FilterNumValue", "FilterPvalEnabled", "FilterPvalue", "FilterDepthEnabled", "FilterDepthValue",
                    "SelectMode", "SelectDisjoint", "AddGOclass"]

    def __init__(self, parent=None, signalManager = None, name='OWGOTermFinder'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Cluster Examples", ExampleTable, self.clusterDataset, Default), ("Reference Examples", ExampleTable, self.referenceDataset, Single + NonDefault), ("Structured Data", DataFiles, self.chipdata, Single + NonDefault)]
        self.outputs = [("Examples", ExampleTable, Default), ("Classified Examples", ExampleTableWithClass, Default), ("Example Selection", ExampleSelection, Default), ("Selected Structured Data", DataFiles, Single + NonDefault)]
        #set default settings
        # annotation
        self.AnnotationFileName = self.GOaspectFileName = None # these are names of files
        self.RecentAnnotations = []
        self.BAnnotationIndx = 0
        self.genesInAnnotationFile = {}
        self.organizmCodes=go.listDownloadedOrganizms()
        # reference
        self.ReferenceType = 0 ## get the reference from the annotation
        # GO
        self.RecentGOaspects = []
        self.BGOaspectIndx = 0
        #
        self.maxDepth = 8
        self.FilterNumEnabled = False
        self.FilterNumValue = 1
        self.FilterPvalEnabled = True
        self.FilterPvalue = 0.50
        self.FilterDepthEnabled = False
        self.FilterDepthValue = 8
        self.SelectMode = 0 # sub graph
        self.SelectDisjoint = False # output inclusive
        self.AddGOclass = False
        # check usage of all evidences
        for etype in go.evidenceTypesOrdered:
            varName = "UseEvidence"+etype 
##            self.settingsList.append( varName)
            code = compile("self.%s = True" % (varName), ".", "single")
            exec(code)

        self.loadSettings()
        self.data = None
        self.chipdata = None
        # check if files exist and remove those that don't
        # check that all files in directories "Annotation" and "GO" are already included
        self.RecentAnnotations = filter(os.path.exists, self.RecentAnnotations)
        self.RecentGOaspects = filter(os.path.exists, self.RecentGOaspects)
        widgetDir = os.path.dirname(os.path.abspath(__file__)) + "/"
        ## add all annotations in "./Annotation" directory
        annotList = glob.glob(widgetDir + 'Annotation/*.annotation')
        for f in annotList:
            if f not in self.RecentAnnotations:
                self.RecentAnnotations.append( f)
        genesInAnnotationFile = {}

        ## add all GOs in "./GO" directory
        GOlist = glob.glob(widgetDir + 'GO/*.go')
        for f in GOlist:
            if f not in self.RecentGOaspects:
                self.RecentGOaspects.append( f)

        # tmp structures - loaded by user
        self.annotation = None
        self.GO = {}
        self.GO["relationTypes"]=dict([("is_a",0),("part_of",1)])

        # received by signals
        self.candidateGeneIDsFromSignal = [] ## list of discrete attributes present in clusterSet data signal
        self.BgeneIDattrIndx = -1 ## index of attribute in candidateGeneIDsFromSignal that was selected to represent the gene IDs
        self.geneIDattr = None ## self.geneIDattr is set accordingly
        # should read from 'GeneName' column in input signal "Examples"
        self.clusterGenes = [] ## ['YPD1', 'WHI4', 'SHS1', 'GCS1', 'HO', 'YDL228C', 'SSB1', 'PTP1', 'BRE4', 'OST4', 'YDL233W', 'GYP7']
        self.clusterData = None
        # should read from 'GeneName' column in input signal "Examples Reference"
        self.referenceGenes = None
        self.referenceData = None

        # calculated from tmp structures and received signals
        # updated by filters
        self.GOIDsFound = [] # sorted by p value, so we know when to stop, if filtering by p. value
        self.significantGOIDs = [] # selected by filters
        self.GOtermValues = {}
        self.dag = None
        self.goLVitem2GOID = {}

        # GUI definition
        self.tabs = QTabWidget(self.controlArea, 'tabWidget')

        # INPUT TAB
        self.inputTab = QVGroupBox(self)
        box = QVButtonGroup("Annotation", self.inputTab)
        box2 = QHButtonGroup(box)
        box2.setMaximumSize(250,50)
        # annotation
        self.annotationCombo = OWGUI.comboBox(box2, self, "BAnnotationIndx", items=[], callback=self.loadAnnotation)
        self.annotationCombo.setMaximumSize(150, 20)
        #self.setFilelist(self.annotationCombo, self.RecentAnnotations)
##        box2.hide()
##        box2.show()
        #self.annotationBrowse = OWGUI.button(box2, self, 'Browse', callback=self.browseAnnotation)
        self.annotationDownload = OWGUI.button(box2, self, "Download", callback=self.downloadAnnotation)
        self.annotationDownload.setMaximumSize(60, 40)
        self.evidencesBox = QVButtonGroup("Evidence codes in annotation", box)
        self.evidenceCheckBoxes = {}
        for etype in go.evidenceTypesOrdered:
            varName = "UseEvidence"+etype
            tmpCB = OWGUI.checkBox(self.evidencesBox, self, varName, etype, box='', tooltip=go.evidenceTypes.get(etype, '?unknown?'), callback=self.findTermsBuildDAG)
            tmpCB.setEnabled(False)
            self.evidenceCheckBoxes[etype] = tmpCB

        # reference
        OWGUI.radioButtonsInBox(self.inputTab, self, 'ReferenceType', ['Annotation', 'Signal'], box='Reference from', callback=self.findTermsBuildDAG)
        self.GOAspectRadioBoxButtons=OWGUI.radioButtonsInBox(self.inputTab, self, "BGOaspectIndx", ["Biological process", "Cellular component", "Molecular function"], box="GO aspect", callback=self.loadGOaspect)
        OWGUI.button(self.inputTab, self, "Download latest GO", callback=self.downloadGO)
        # gene name attribute
        box = QHButtonGroup("Gene ID Attribute", self.inputTab)
        box.setMaximumSize(250, 50)
        self.geneIDAttrCombo = OWGUI.comboBox(box, self, 'BgeneIDattrIndx', items=[], callback=self.geneIDchanged)
        self.geneIDAttrCombo.setMaximumSize(160, 20)
        self.setGeneIDAttributeList()
        self.tabs.insertTab(self.inputTab, "Input")

        # FILTER TAB
        self.filterTab = QVGroupBox(self)
        box = QVButtonGroup("Filter GO Term Nodes", self.filterTab)
        #
        OWGUI.checkBox(box, self, 'FilterNumEnabled', "Number of instances", callback=self.setFilterNumEnabled)
        self.sliderFilterNumValue = OWGUI.qwtHSlider(box, self, 'FilterNumValue', label='#:', labelWidth=33, minValue=1, maxValue=1000, step=1.0, precision=1, ticks=0, maxWidth=80, callback=self.runFilters)
        if not self.FilterNumEnabled:
            self.sliderFilterNumValue.box.setDisabled(1)
        #
        OWGUI.checkBox(box, self, 'FilterPvalEnabled', "p value", callback=self.setFilterPvalEnabled)
        self.sliderFilterPvalue = OWGUI.qwtHSlider(box, self, 'FilterPvalue', label='p:', labelWidth=33, minValue=0.0, maxValue=1.0, step=0.001, precision=3.0, ticks=0, maxWidth=80, callback=self.runFilters)
        if not self.FilterPvalEnabled:
            self.sliderFilterPvalue.box.setDisabled(1)
        #
        OWGUI.checkBox(box, self, 'FilterDepthEnabled', "GO depth", callback=self.setFilterDepthEnabled)
        self.sliderFilterDepthValue = OWGUI.qwtHSlider(box, self, 'FilterDepthValue', label='p:', labelWidth=33, minValue=0.0, maxValue=100, step=1.0, precision=1.0, ticks=0, maxWidth=80, callback=self.runFilters)
        if not self.FilterDepthEnabled:
            self.sliderFilterDepthValue.box.setDisabled(1)
        self.tabs.insertTab(self.filterTab, "Filter")

        # SELECT TAB
        self.selectTab = QVGroupBox(self)
        OWGUI.radioButtonsInBox(self.selectTab, self, 'SelectMode', ['Directly or Indirectly', 'Directly'], box='Annotated Genes', callback=self.viewSelectionChanged)
        box = QVButtonGroup('Output', self.selectTab)
        OWGUI.checkBox(box, self, 'SelectDisjoint', 'Disjoint/Inclusive', callback=self.viewSelectionChanged)
        OWGUI.checkBox(box, self, 'AddGOclass', 'Add GO term as new class', callback=self.viewSelectionChanged)
        self.tabs.insertTab(self.selectTab, "Select")

        # ListView for DAG, and table for significant GOIDs
        self.DAGcolumns = ['GO term', 'Cluster frequency', 'Reference frequency', 'p value', 'Genes']
        self.layout=QVBoxLayout(self.mainArea)
        self.splitter = QSplitter(QSplitter.Vertical, self.mainArea)
        self.layout.add(self.splitter)

        # list view
        self.goLV = QListView(self.splitter)
        self.goLV.setMultiSelection(1)
        self.goLV.setAllColumnsShowFocus(1)
        self.goLV.addColumn(self.DAGcolumns[0])
        self.goLV.setColumnWidth(0, 300)
        self.goLV.setColumnWidthMode(0, QListView.Manual)
        self.goLV.setColumnAlignment(0, QListView.AlignLeft)
        self.goLV.setSorting(-1)
        for dagColumnTitle in self.DAGcolumns[1:]:
            col = self.goLV.addColumn(dagColumnTitle)
            self.goLV.setColumnWidth(col, 100)
            self.goLV.setColumnWidthMode(col, QListView.Manual)
            self.goLV.setColumnAlignment(col, QListView.AlignCenter)
        self.connect(self.goLV, SIGNAL("selectionChanged()"), self.viewSelectionChanged)

        # table of significant GO terms
        self.sigTermsTable = QTable(self.splitter)
        self.sigTermsTable.setNumCols(5)
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
        self.connect(self.sigTermsTable, SIGNAL("selectionChanged()"), self.tableSelectionChanged)
        self.splitter.show()

        self.resize(1000, 800)
        self.layout.activate() # this is needed to scale the widget correctly

        self.organizmGenes={}
        self.setAnnotationCodesList()

    def geneIDchanged(self):
        if len(self.candidateGeneIDsFromSignal) > self.BgeneIDattrIndx:
            self.geneIDAttrCombo.setCurrentItem(self.BgeneIDattrIndx)
            self.geneIDattr = self.candidateGeneIDsFromSignal[self.BgeneIDattrIndx]
        else:
            self.geneIDattr = None
        if DEBUG: print "changing geneID attribute to: " + str(self.geneIDattr)
        self.clusterDatasetChanged()
        self.referenceDatasetChanged()
        self.findTermsBuildDAG()
        self.geneIDAttrCombo.setDisabled(len(self.candidateGeneIDsFromSignal) == 0)

    def setGeneIDAttributeList(self):
        ## refresh the list
        self.geneIDAttrCombo.clear()
        for f in self.candidateGeneIDsFromSignal:
            self.geneIDAttrCombo.insertItem(str(f.name))
        self.geneIDAttrCombo.setDisabled(len(self.candidateGeneIDsFromSignal) == 0)

    def setAnnotationCodesList(self):
        self.organizmCodes=go.listDownloadedOrganizms()
        self.annotationCombo.clear()
        for o in self.organizmCodes:
            self.annotationCombo.insertItem(o)
        
    def updateGeneID2annotationfile(self):
        ## geneID is key, item is list of indexes in self.RecentAnnotations that have that geneID
        self.progressBarInit()
        self.geneID2annotationfile = {}
        cn = 0
        allcn = len(self.RecentAnnotations)
        for f in self.RecentAnnotations:
            if f not in self.genesInAnnotationFile.keys():
                try:
                    loadedAnnotation = cPickle.load(open(f, 'r'))
                    self.genesInAnnotationFile[f] = loadedAnnotation['gene2GOID'].keys()
                except:
                    self.genesInAnnotationFile[f] = []
                    continue
            genesInAnnotation = self.genesInAnnotationFile[f]
            for geneID in genesInAnnotation:
                tmpl = self.geneID2annotationfile.get(geneID, [])
                if cn not in tmpl:
                    tmpl.append(cn)
                    self.geneID2annotationfile[geneID] = tmpl
            cn += 1
            self.progressBarSet(int(round(100.0*cn/allcn)))
        self.progressBarFinished()

    ## this is called only when new data token is received
    def findMostAppropriateGeneIDandAnnotation(self):
        if self.clusterData == None:
            self.candidateGeneIDsFromSignal = []
            self.BgeneIDattrIndx = -1
            self.geneIDattr = None
            self.setGeneIDAttributeList()
            return

        ## all discrete and string type attributes are good candidates
        self.candidateGeneIDsFromSignal = [a for a in self.clusterData.domain.attributes + self.clusterData.domain.getmetas().values() if a.varType == orange.VarTypes.Discrete or a.varType == orange.VarTypes.Other or a.varType == orange.VarTypes.String]
        self.setGeneIDAttributeList()
        self.geneIDAttrCombo.setDisabled(1)

        if not self.organizmGenes:        
            self.organizmGenes=dict([(o,Set(go.getCachedGeneNames(o))) for o in self.organizmCodes])
        cn={}
        for attr in self.candidateGeneIDsFromSignal:
            vals=[str(e[attr]) for e in self.clusterData]
            for organizm, set in self.organizmGenes.items():
                l=filter(lambda a: a in set, vals)
                cn[(attr,organizm)]=len(l)
        cn=cn.items()
        cn.sort(lambda a,b:-cmp(a[1],b[1]))
        bestAttr, organizm=cn[0][0]
        self.BAnnotationIndx=self.organizmCodes.index(organizm)
        self.BgeneIDattrIndx=self.candidateGeneIDsFromSignal.index(bestAttr)
        self.geneIDchanged()

    ##########################################################################
    # handling of input/output signals
    def chipdata(self, chipdata):
        if chipdata:
            self.chipdata = chipdata
#            self.clusterData =  self.chipdata[0][1][0]
            self.clusterDataset(self.chipdata[0][1][0], 0)
        else:
            self.chipdata = None
            # chip
            self.send("Example Selection", None)
            self.send("Selected Structured Data", None)
            # "regular" data
            self.send("Examples", None)
            self.send("Classified Examples", None)

    def clusterDataset(self, data):
        self.clusterData = data
        self.findMostAppropriateGeneIDandAnnotation()
        self.clusterDatasetChanged()
        #pywin.debugger.set_trace()
        self.loadGOaspect(forced=0) ## usually only at the first run, the GO aspect is not loaded
##        self.findTermsBuildDAG()

    def clusterDatasetChanged(self):
        self.clusterGenes = []
        if self.clusterData <> None:
            if DEBUG: print "clusterDatasetChanged, self.geneIDattr: " + str(self.geneIDattr)
            if self.geneIDattr in self.candidateGeneIDsFromSignal:
                for e in self.clusterData:
                    g = str(e[self.geneIDattr])
                    if g not in self.clusterGenes:
                        self.clusterGenes.append( g)
        if DEBUG: print "input cluster genes: " + str(len(self.clusterGenes))
        ## self.findTermsBuildDAG() need to call it, if you call clusterDatasetChanged directly

    def referenceDataset(self, data):
        self.referenceGenes = None
        self.referenceData = data
        self.referenceDatasetChanged()
        self.findTermsBuildDAG()

    def referenceDatasetChanged(self):
        if DEBUG: print "reference: " + str(self.referenceData)
        if self.referenceData <> None:
            self.referenceGenes = []
            dattrs = [a for a in self.referenceData.domain.attributes + self.referenceData.domain.getmetas().values()]
            if self.geneIDattr in dattrs:
                for e in self.referenceData:
                    g = str(e[self.geneIDattr])
                    if g not in self.referenceGenes:
                        self.referenceGenes.append( g)
        else:
            self.referenceGenes = None
        ## self.findTermsBuildDAG() need to call it, if you call referenceDatasetChanged directly

    def tableSelectionChanged(self):
        self.progressBarInit()
        tot = len(self.significantGOIDs)
        totcn = 0
        for i in range(len(self.significantGOIDs)):
            b = self.sigTermsTable.isRowSelected(i, False)
            GOID = self.significantGOIDs[i]
            for (li, liGOID) in self.goLVitem2GOID.items():
                if liGOID == GOID:
                    self.goLV.setSelected(li, b)
                    if b: self.goLV.ensureItemVisible(li)
            self.progressBarSet(int(round(totcn*100.0/tot)))
            totcn += 1
        self.progressBarFinished()

    def viewSelectionChanged(self):
        geneToGOterm = {}
        allGOterms = []
        for li in self.goLVitem2GOID.keys():
            if li.isSelected():
                GOID = self.goLVitem2GOID.get(li, None)
                GOterm, x, G, pval, genesInGOID, genesInGOIDdirect = self.GOtermValues.get(GOID, (GOID+'?', '', '', '', [], []))
                if GOID == 'root': ## put real aspect instead of 'root'
                    GOterm = self.GO.get('aspect', GOID+'?')
                
                if GOterm not in allGOterms:
                    allGOterms.append( GOterm)

                ## make gene -> GOterm annotations only for some genes; depending on the selection type
                if self.SelectMode == 1: 
                    geneList = genesInGOIDdirect # node specific: use just the direct annotations
                else:
                    geneList = genesInGOID # subgraph: use both directly and indirectly annotated genes

                for gene in geneList:
                    tmpl = geneToGOterm.get(gene, [])
                    if GOterm not in tmpl:
                        tmpl.append(GOterm)
                        geneToGOterm[gene] = tmpl
        self.sendSelectedGenesData(geneToGOterm, allGOterms)

    def sendSelectedGenesData(self, gtg, agot):
        if self.clusterData:
            # class value; GO terms
            # new domain
            if self.AddGOclass:
                newclass = orange.EnumVariable("GO class", values=agot)
                newdomain = orange.Domain( self.clusterData.domain.attributes, newclass)
            else:
                newdomain = orange.Domain( self.clusterData.domain)
            metas = self.clusterData.domain.getmetas()
            for key in metas:
                newdomain.addmeta(key, metas[key])
            # new exampletable into where to put the filtered examples
            newdata = orange.ExampleTable(newdomain)
            sel = []
            selDescription = []
            for e in self.clusterData:
                g = str(e[self.geneIDattr])
                geneTermList = gtg.get(g, [])
                if self.SelectDisjoint and len(geneTermList) > 1: ## this gene should be omitted, because belongs to many GOterms
                    sel.append( 0)
                    continue
                sel.append(int(g in gtg.keys()))
                for goterm in geneTermList:
                    if len(selDescription) < 6 and goterm not in selDescription:
                        selDescription.append( goterm)
                    nex = orange.Example(newdomain, e)
                    if self.AddGOclass:
                        nex.setclass(goterm)
                    newdata.append( nex)

            if self.chipdata:
                P = []
                for (strainname, tmpdata) in self.chipdata:
                    dataP = [d.select(sel) for d in tmpdata]
                    for i in range(len(tmpdata)):
                        dataP[i].name = tmpdata[i].name
                    P.append((strainname, dataP))
                self.send("Selected Structured Data", P)
            if newdata.domain.classVar:
                self.send("Classified Examples", newdata)
            else:
                self.send("Classified Examples", None)
            self.send("Examples", newdata)
            if self.chipdata:
                if len(selDescription) == 6:
                    selDescription[-1] = "..."
                self.send("Example Selection", ("GO: " + ", ".join(selDescription), sel))
            else:
                self.send("Selected Structured Data", None)
                self.send("Example Selection", None) 
        else:
            self.send("Example Selection", None)
            self.send("Selected Structured Data", None)
            self.send("Classified Examples", None)
            self.send("Examples", None)

    ##########################################################################
    # callback functions
    
    def loadAnnotation(self):
        if DEBUG: print "loadAnnotation"
        """
        self.annotation, self.BAnnotationIndx = self.loadRemember(self.RecentAnnotations, self.annotationCombo, self.BAnnotationIndx)
        fn = str(self.RecentAnnotations[0]) ## the loaded one is (becomes) always moved to 0 position
        self.genesInAnnotationFile[fn] = self.annotation['gene2GOID'].keys() ## update, in case the file content changed
        self.updateEvidences()"""
        self.progressBarInit()
        go.loadAnnotation(self.organizmCodes[self.BAnnotationIndx], progressCallback=self.progressBarSet)
        self.progressBarFinished()
        self.updateEvidences()
        self.findTermsBuildDAG()        

    def loadGOaspect(self, forced=1):
        if DEBUG: print "loadGOaspect"
        ## load if forced, or if index has changed
        ## if forced = 0 and index has not changed (still is 0) then don't reload the annotation data
        """if forced == 1 or self.BGOaspectIndx <> 0 or self.GO == None:
            if DEBUG: print "1:", str(self.RecentGOaspects) + "," + str(self.BGOaspectIndx)
            self.GO, self.BGOaspectIndx = self.loadRemember(self.RecentGOaspects, self.GOaspectCombo, self.BGOaspectIndx)
            if DEBUG: print "2:", str(self.RecentGOaspects) + "," + str(self.BGOaspectIndx)
            self.updateEvidences()
            """
        self.progressBarInit()
        go.loadGO(progressCallback=self.progressBarSet)
        go.loadAnnotation(self.organizmCodes[self.BAnnotationIndx], progressCallback=self.progressBarSet)
        self.progressBarFinished()
        self.updateEvidences()
        self.findTermsBuildDAG()

    def updateEvidences(self):
        """
        if not(self.annotation) or not(self.GO): ## if data missing, just disable everything
            for (etype, tmpCB) in self.evidenceCheckBoxes.items():
                tmpCB.setText(etype)
                tmpCB.setEnabled(False)
            return
        """
        if not go.loadedAnnotation:
            for (etype, tmpCB) in self.evidenceCheckBoxes.items():
                tmpCB.setText(etype)
                tmpCB.setEnabled(False)
            return

        evidenceCount=dict([(etype,0) for etype in go.evidenceDict.keys()])
        evidenceGenes=dict([(etype,Set()) for etype in go.evidenceDict.keys()])
        for ann in go.loadedAnnotation.annotationList:
            evidenceCount[ann.evidence]+=1
            evidenceGenes[ann.evidence].add(ann.geneName)
        for (etype, tmpCB) in self.evidenceCheckBoxes.items():
            if evidenceCount[etype]:
                tmpCB.setEnabled(True)
                tmpCB.setText('%s: %d annots (%d genes)' % (etype, evidenceCount[etype], len(evidenceGenes[etype])))
            else:
                tmpCB.setEnabled(False)
                tmpCB.setText(etype)
            
    def setFilterNumEnabled(self):
        self.sliderFilterNumValue.box.setDisabled(not self.FilterNumEnabled)
        self.runFilters()

    def setFilterPvalEnabled(self):
        self.sliderFilterPvalue.box.setDisabled(not self.FilterPvalEnabled)
        self.runFilters()

    def setFilterDepthEnabled(self):
        self.sliderFilterDepthValue.box.setDisabled(not self.FilterDepthEnabled)
        self.runFilters()

    def downloadAnnotation(self):
        organizmList=go.listOrganizms()
        w=DownloadDialog(self, None, "", organizmList)

    def downloadGO(self):
        self.progressBarInit()
        go.downloadGO(self.progressBarSet)
        self.progressBarFinished()

    ##########################################################################
    # GO DAG calculations and filtering
    def runFilters(self):
        self.significantGOIDs = [] ## significant GOID to display in GO
        for (p, x, GOID) in self.GOIDsFound:
            if self.FilterPvalEnabled and p > self.FilterPvalue:
                break ## end of significant GO terms reached
            if self.FilterNumEnabled and x < self.FilterNumValue:
                continue ## not worth mentioning
            self.significantGOIDs.append(GOID)

        ##self.dag = GOlib.createGODAGtoDisplay(self.GO, self.GOtermValues.keys(), self.significantGOIDs)
        daglist=go.extractGODAG(self.GOTermFinderResult.keys())
        if self.FilterDepthEnabled:
            daglist=go.DAGFilterForDepth(daglist, self.FilterDepthValue)
            daglist=go.extractGODAG([t.id for t in daglist])
        self.dag={"root":[("","")]}
        if not daglist:
            self.updateDAG()
            return
        self.dag=dict([(term.id,[]) for term in daglist])
        for term in daglist:
            if term.parents:
                for parent in term.parents:
                    rType=term.rType[parent]
                    if parent in self.dag:
                        self.dag[parent].append((term.id, rType))
                    else:
                        print "er 1!", parent
                        self.dag[parent]=[(term.id, rType)]
            else:
                root=term
        self.dag["root"]=self.dag[root.id]
        for term in daglist:
            if term.id not in self.dag:
                print "er 2",term.id
                self.dag[term.id]=[]
        del self.dag[root.id]
            
        """if self.FilterDepthEnabled:
            self.dag = GOlib.DAGfilterForDepth(self.dag, 'root', self.FilterDepthValue)"""
        self.updateDAG()

    def findTermsBuildDAG(self):
        self.dag = {}
        if DEBUG: print "findTermsBuildDAG, self.annotation: " + str(self.annotation <> None)
        if DEBUG: print "findTermsBuildDAG, self.GO: " + str(self.GO <> None)
        if not self.clusterData:
            self.significantGOIDs=[]
            self.updateDAG()
            return
        if go.loadedAnnotation and go.loadedGO:
            self.progressBarInit()
            evidences = [etype for (etype, tmpCB) in self.evidenceCheckBoxes.items() if tmpCB.isChecked()]
            
            aspect=["biological_process","cellular_component","molecular_function"][self.BGOaspectIndx]
            self.GO["aspect"]=aspect
            #pywin.debugger.set_trace()
            if self.ReferenceType==0:
                referenceGenesSet=[]    #will use all the genes in the loaded annotation
            else:
                referenceGenesSet=self.referenceGenes
            for g in self.clusterGenes:
                print g
            self.GOTermFinderResult=result=go.GOTermFinder(self.clusterGenes, referenceGenesSet, evidences, False, aspect, self.progressBarSet)
            directAnnotations=go.findTerms(self.clusterGenes, aspect=[aspect], directAnnotationOnly=True, evidenceCodes=evidences, reportEvidence=False, progressCallback=self.progressBarSet)
            if not self.GOTermFinderResult:
                self.significantGOIDs=[]
                self.dag={}
                print "no GO terms"
                self.updateDAG()
                return 
            self.GOtermValues=dict([(GOID,(go.loadedGO.termDict[GOID].name, len(genes), numRef, pvalue, genes, directAnnotations.get(GOID,[]))) for GOID,(genes, pvalue, numRef) in self.GOTermFinderResult.items()])
            root=filter(lambda t:not go.loadedGO.termDict[t].parents, self.GOtermValues.keys())[0]
            #print root, go.loadedGO.termDict[root]
            self.GOtermValues["root"]=self.GOtermValues[root]
            #del self.GOtermValues[root]
            ##fill in sorted GOIDsFound (pValue, i, ID)
            self.GOIDsFound=[(value[1],len(value[0]), GOID) for GOID,value in result.items()]
            self.GOIDsFound.sort()
            
            n=len(self.clusterGenes); N=self.referenceGenes and len(self.referenceGenes) or len(go.loadedAnnotation.geneNames)
##            print n, N

            ## find the max number of cluster gene istances in a GO term
            maxNumIstances=max([len(t[0]) for t in result.values()])
            ## create a DAG with all the possible nodes
            ## and find the max depth of the DAG
            sigGOIDs = [goid for (_, _, goid) in self.GOIDsFound]
            dag=go.extractGODAG(sigGOIDs)
            
            #pywin.debugger.set_trace()
            maxDepth=go.DAGDepth(dag)
                              
            ## update the filter controls
            self.sliderFilterNumValue.setRange(1, maxNumIstances, 1)
            self.sliderFilterDepthValue.setRange(0, maxDepth, 1)

            self.runFilters()
            self.progressBarSet(95)
            self.updateDAG()
            self.progressBarFinished()

    ##########################################################################
    # drawing
    def updateDAG(self, updateonly=0):
        def walkupdate(listviewitem):
            GOID = self.goLVitem2GOID[listviewitem]
            GOterm, x, G, pval, genesInGOID, genesInGOIDdirect = self.GOtermValues.get(GOID, (GOID+'?', '', '', '', [], []))
            if GOID == 'root': ## put real aspect instead of 'root'
                GOterm = self.GO.get('aspect', GOID+'?')
            if len(genesInGOID):
                genesInGOIDstr = str(genesInGOID[0])
            else:
                genesInGOIDstr = ''
            for gene in genesInGOID[1:]:
                genesInGOIDstr += ", " + str(gene)
        
            if pval: pval = "%.4g" % pval
            vals = [GOterm, len(genesInGOID), G, pval, genesInGOIDstr] 
            for i in range(len(vals)):
                listviewitem.setText(i, str(vals[i]))

            child = listviewitem.firstChild()
            while child:
                walkupdate(child)
                child = child.nextSibling()

        def walkcreate(node, parent):
            for (childNode, rtype) in self.dag.get(node, []):
                bd = str(childNode)
                li = QListViewItem(parent, bd)
                li.setOpen(1)
                self.goLVitem2GOID[li] = childNode
                walkcreate(childNode, li)

        if not(self.dag):           
            self.goLV.clear()
            return
        
        self.goLV.setRootIsDecorated(1)
        if not updateonly:
            self.goLV.clear()
            self.goLVitem2GOID = {}
            self.GOid2LVitem = {}
            li = QListViewItem(self.goLV, 'root')
            li.setOpen(1)
            self.goLVitem2GOID[li] = 'root'
            walkcreate('root', li)
        walkupdate(self.goLV.firstChild())
        self.goLV.show()

        # update table of significant/filtered Terms
        self.sigTermsTable.setNumRows(len(self.significantGOIDs))
        #print self.significantGOIDs
        for i in range(len(self.significantGOIDs)): ## sorted by the p value
            GOID = self.significantGOIDs[i]
            GOterm, x, G, pval, genesInGOID, genesInGOIDdirect = self.GOtermValues.get(GOID, (GOID+'?', '', '', '', [], []))
            if GOID == 'root': ## put real aspect instead of 'root'
                GOterm = self.GO.get('aspect', GOID+'?')

            if len(genesInGOID):
                genesInGOIDstr = str(genesInGOID[0])
            else:
                genesInGOIDstr = ''
            for gene in genesInGOID[1:]:
                genesInGOIDstr += ", " + str(gene)
            
            if pval: pval = "%.4g" % pval
            vals = [GOterm, len(genesInGOID), G, pval, genesInGOIDstr]
            for j in range(len(vals)):
                self.sigTermsTable.setText(i, j, str(vals[j]))

class DownloadDialog(OWWidget):
    def __init__(self, parent=None, signalManager=None, name="Download annotations", codes=[]):
        OWWidget.__init__(self, parent, signalManager, name)
        self.codesInd=0
        self.codes=codes
        self.master=parent
        self.rBox=OWGUI.comboBox(self.controlArea, self, "codesInd", items=codes, box="Organizm Codes")
        OWGUI.button(self.controlArea, self, "OK", callback=self.download)
        self.resize(100,100)
        self.show()

    def download(self):
        self.progressBarInit()
        go.downloadAnnotation(self.codes[self.codesInd], self.progressBarSet)
        self.progressBarFinished()
        self.master.setAnnotationCodesList()
        
if __name__=="__main__":
    import orange
    a = QApplication(sys.argv)
    ow = OWGOTermFinder()
    a.setMainWidget(ow)

##    d = orange.ExampleTable('testClusterSet.tab', dontCheckStored=1)
##    d = orange.ExampleTable('hjSmall.tab', dontCheckStored=1)
    d = orange.ExampleTable('dicty2.tab', dontCheckStored=1)
    ow.show()
    #pywin.debugger.set_trace()
    ow.clusterDataset(d)
    #ow.clusterDataset(None)
    a.exec_loop()
    ow.saveSettings()
