"""
<name>GO Term Finder</name>
<description>GO Term Finder</description>
<icon>icons/GOTermFinder.png</icon>
<priority>100</priority>
"""

import orange, math, glob
import GOlib ## function needed to handle the GO and annotation
import OWGUI
from qt import *
from qtcanvas import *
from OWWidget import *
from OWOptions import *
from qttable import *
from qwt import *

DEBUG = 0

class OWGOTermFinder(OWWidget):	
    settingsList = ["AnnotationFileName", "RecentAnnotations", "ReferenceType", "RecentGOaspects",
                    "FilterNumEnabled", "FilterNumValue", "FilterPvalEnabled", "FilterPvalue", "FilterDepthEnabled", "FilterDepthValue",
                    "SelectMode", "SelectDisjoint", "AddGOclass"]

    def __init__(self, parent=None, name='OWGoTermFinder'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, name, 'GO Term Finder', FALSE, FALSE) 

        self.inputs = [("Cluster Examples", ExampleTable, self.clusterDataset, 0), ("Reference Examples", ExampleTable, self.referenceDataset, 0)]
        self.outputs = [("Examples", ExampleTable), ("Classified Examples", ExampleTableWithClass)]

        #set default settings
        # annotation
        self.AnnotationFileName = self.GOaspectFileName = None # these are names of files
        self.RecentAnnotations = []
        self.BAnnotationIndx = 0
        self.genesInAnnotationFile = {}
        # reference
        self.ReferenceType = 0 ## get the reference from the annotation
        # GO
        self.RecentGOaspects = []
        self.BGOaspectIndx = 0
        #
        self.FilterNumEnabled = False
        self.FilterNumValue = 1
        self.FilterPvalEnabled = True
        self.FilterPvalue = 0.05
        self.FilterDepthEnabled = False
        self.FilterDepthValue = 8
        self.SelectMode = 0 # sub graph
        self.SelectDisjoint = False # output inclusive
        self.AddGOclass = False
        # check usage of all evidences
        for etype in GOlib.evidenceTypesOrdered:
            varName = "UseEvidence"+etype 
##            self.settingsList.append( varName)
            code = compile("self.%s = True" % (varName), ".", "single")
            exec(code)

        self.loadSettings()
        self.data = None
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
        self.GO = None

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
        self.setFilelist(self.annotationCombo, self.RecentAnnotations)
##        box2.hide()
##        box2.show()
        self.annotationBrowse = OWGUI.button(box2, self, 'Browse', callback=self.browseAnnotation)
        self.annotationBrowse.setMaximumSize(50, 30)
        self.evidencesBox = QVButtonGroup("Evidence codes in annotation", box)
        self.evidenceCheckBoxes = {}
        for etype in GOlib.evidenceTypesOrdered:
            varName = "UseEvidence"+etype
            tmpCB = OWGUI.checkBox(self.evidencesBox, self, varName, etype, box='', tooltip=GOlib.evidenceTypes.get(etype, '?unknown?'), callback=self.findTermsBuildDAG)
            tmpCB.setEnabled(False)
            self.evidenceCheckBoxes[etype] = tmpCB

        # reference
        OWGUI.radioButtonsInBox(self.inputTab, self, 'ReferenceType', ['From Annotation', 'From Signal'], box='Reference', callback=self.findTermsBuildDAG)
        # GO aspects
        box = QHButtonGroup("GO Aspect", self.inputTab)
        box.setMaximumSize(250, 50)
        self.GOaspectCombo = OWGUI.comboBox(box, self, 'BGOaspectIndx', items=[], callback=self.loadGOaspect)
        self.GOaspectCombo.setMaximumSize(160, 20)
        self.setFilelist(self.GOaspectCombo, self.RecentGOaspects)
        self.GOaspectBrowse = OWGUI.button(box, self, 'Browse', callback=self.browseGOaspect)
        self.GOaspectBrowse.setMaximumSize(50, 30)
        # gene name attribute
        box = QHButtonGroup("Gene ID attribute", self.inputTab)
        box.setMaximumSize(250, 50)
        self.geneIDAttrCombo = OWGUI.comboBox(box, self, 'BgeneIDattrIndx', items=[], callback=self.geneIDchanged)
        self.geneIDAttrCombo.setMaximumSize(160, 20)
        self.setGeneIDAttributeList()
        self.tabs.insertTab(self.inputTab, "Input")

        # FILTER TAB
        filterTab = QVGroupBox(self)
        box = QVButtonGroup("Filter GO Term Nodes", filterTab)
        #
        OWGUI.checkBox(box, self, 'FilterNumEnabled', "Number of instances", callback=self.setFilterNumEnabled)
        self.sliderFilterNumValue = OWGUI.qwtHSlider(box, self, 'FilterNumValue', label='#:', labelWidth=33, minValue=1, maxValue=1000, step=1.0, precision=1, ticks=0, maxWidth=80, callback=self.runFilters)
        if not self.FilterNumEnabled:
            self.sliderFilterNumValue.box.setDisabled(1)
        #
        OWGUI.checkBox(box, self, 'FilterPvalEnabled', "p. value", callback=self.setFilterPvalEnabled)
        self.sliderFilterPvalue = OWGUI.qwtHSlider(box, self, 'FilterPvalue', label='p:', labelWidth=33, minValue=0.0, maxValue=1.0, step=0.001, precision=3.0, ticks=0, maxWidth=80, callback=self.runFilters)
        if not self.FilterPvalEnabled:
            self.sliderFilterPvalue.box.setDisabled(1)
        #
        OWGUI.checkBox(box, self, 'FilterDepthEnabled', "GO depth", callback=self.setFilterDepthEnabled)
        self.sliderFilterDepthValue = OWGUI.qwtHSlider(box, self, 'FilterDepthValue', label='p:', labelWidth=33, minValue=0.0, maxValue=100, step=1.0, precision=1.0, ticks=0, maxWidth=80, callback=self.runFilters)
        if not self.FilterDepthEnabled:
            self.sliderFilterDepthValue.box.setDisabled(1)
        self.tabs.insertTab(filterTab, "Filter")

        # SELECT TAB
        selectTab = QVGroupBox(self)
        OWGUI.radioButtonsInBox(selectTab, self, 'SelectMode', ['Subgraph', 'Node specific'], box='Mode', callback=self.viewSelectionChanged)
        box = QVButtonGroup('Output', selectTab)
        OWGUI.checkBox(box, self, 'SelectDisjoint', 'Disjoint/Inclusive', callback=self.viewSelectionChanged)
        OWGUI.checkBox(box, self, 'AddGOclass', 'Add GO term as new class', callback=self.viewSelectionChanged)
        self.tabs.insertTab(selectTab, "Select")

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

    def updateGeneID2annotationfile(self):
        ## geneID is key, item is list of indexes in self.RecentAnnotations that have that geneID
        self.geneID2annotationfile = {}
        cn = 0
        for f in self.RecentAnnotations:
            if f not in self.genesInAnnotationFile.keys():
                loadedAnnotation = cPickle.load(open(f, 'r'))
                self.genesInAnnotationFile[f] = loadedAnnotation['gene2GOID'].keys()
            genesInAnnotation = self.genesInAnnotationFile[f]
            for geneID in genesInAnnotation:
                tmpl = self.geneID2annotationfile.get(geneID, [])
                if cn not in tmpl:
                    tmpl.append(cn)
                    self.geneID2annotationfile[geneID] = tmpl
            cn += 1

    ## this is called only when new data token is received
    def findMostAppropriateGeneIDandAnnotation(self):
        if self.clusterData == None:
            self.candidateGeneIDsFromSignal = []
            self.BgeneIDattrIndx = -1
            self.geneIDattr = None
            self.setGeneIDAttributeList()
            return

        ## all discrete and string type attributes are good candidates
        self.candidateGeneIDsFromSignal = [a for a in self.clusterData.domain.attributes + self.clusterData.domain.getmetas().values() if a.varType == orange.VarTypes.Discrete or a.varType == orange.VarTypes.Other]
        self.setGeneIDAttributeList()
        self.geneIDAttrCombo.setDisabled(1)

        ## check if there are new annotation files present
        ## remove from self.geneID2annotationfile those not present in the RecentAnnotations list
        self.updateGeneID2annotationfile() 

        ## for each attribute look how many genesID are there, that are also present in geneID2annotationfile
        ## if current self.geneIDattr has count 0
        ## then select attribute with highest count
        ## else keep self.geneIDattr

        ## when best attribute selected, check if the loaded annotation is ok
        ## otherwise suggest the most appropriate annotation
        bestAttr = '' ## key is attribute, item is number of recognized geneIDs
        bestCn = 0
        bestAnnotation = 0
        lst = self.candidateGeneIDsFromSignal
        if self.geneIDattr <> None and self.geneIDattr in self.candidateGeneIDsFromSignal: lst = [self.geneIDattr] + lst

        for attr in lst:
            vals = [ex[attr] for ex in self.clusterData]

            ## calculate the frequency of each annotation file to which this geneID belongs to
            annotationFrequency = {}
            cn = 0
            for v in vals:
                v = str(v)
                i = self.geneID2annotationfile.get(v, -1) ## -1, not present
                if i <> -1:
                    for ai in i:
                        af = annotationFrequency.get(ai, 0)
                        annotationFrequency[ai] = af + 1
                    cn += 1
            if cn > bestCn or (cn > 0 and attr == self.geneIDattr):
                bestAttr = attr
                bestCn = cn
                afs = [(f, anindex) for (anindex, f) in annotationFrequency.items()]
                if len(afs) > 0:
                    afs.sort()
                    afs.reverse() ## most frequent first
                    bestAnnotation = afs[0][1]
                else:
                    bestAnnotation = 0 ## keep current
        if DEBUG: print "best attribute: " + str(bestAttr) + " with " + str(bestCn) + " gene IDs from annotations"
        if DEBUG: print "bestAnnotation: " + str(self.RecentAnnotations[bestAnnotation])

        self.geneIDattr = bestAttr
        try:
            self.BgeneIDattrIndx = self.candidateGeneIDsFromSignal.index(self.geneIDattr)
        except:
            self.BgeneIDattrIndx = 0

        ## load annotation if a better one found
        if bestAnnotation <> 0 or self.annotation == None:
            self.BAnnotationIndx = bestAnnotation
            self.annotation, self.BAnnotationIndx = self.loadRemember(self.RecentAnnotations, self.annotationCombo, self.BAnnotationIndx)
            fn = self.RecentAnnotations[0] ## the loaded one is (becomes) always moved to 0 position
            self.genesInAnnotationFile[fn] = self.annotation['gene2GOID'].keys() ## update, in case the file content changed
##            self.annotationCombo.setCurrentItem(self.BAnnotationIndx)

        ## select the geneID, and rerun the GO term finding

        self.geneIDchanged()

    def setFilelist(self, filecombo, fileList):
        filecombo.clear()
        if fileList != []:
            for file in fileList:
                (dir, filename) = os.path.split(file)
                #leave out the path
                fnToDisp = filename
                filecombo.insertItem(fnToDisp)
            filecombo.setDisabled(False)
        else:
            filecombo.insertItem("(none)")
            filecombo.setDisabled(True)

    ##########################################################################
    # handling of input/output signals
    def clusterDataset(self, data, id):
        self.clusterData = data
        self.findMostAppropriateGeneIDandAnnotation()
        self.clusterDatasetChanged()
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

    def referenceDataset(self, data, id):
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
                for e in data:
                    g = str(e[self.geneIDattr])
                    if g not in self.referenceGenes:
                        self.referenceGenes.append( g)
        else:
            self.referenceGenes = None
        ## self.findTermsBuildDAG() need to call it, if you call referenceDatasetChanged directly

    def tableSelectionChanged(self):
        for i in range(len(self.significantGOIDs)):
            b = self.sigTermsTable.isRowSelected(i, False)
            GOID = self.significantGOIDs[i]
            for (li, liGOID) in self.goLVitem2GOID.items():
                if liGOID == GOID:
                    self.goLV.setSelected(li, b)

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
        if self.clusterData:
            # class value; GO terms
            # new domain
            if self.AddGOclass:
                newclass = orange.EnumVariable("GO class", values=allGOterms)
                newdomain = orange.Domain( self.clusterData.domain, newclass)
            else:
                newdomain = orange.Domain( self.clusterData.domain)
            metas = self.clusterData.domain.getmetas()
            for key in metas:
                newdomain.addmeta(key, metas[key])
            # new exampletable into where to put the filtered examples
            newdata = orange.ExampleTable(newdomain)
            for e in self.clusterData:
                g = str(e[self.geneIDattr])
                geneTermList = geneToGOterm.get(g, [])
                if self.SelectDisjoint and len(geneTermList) > 1: ## this gene should be omitted, because belongs to many GOterms
                    continue
                for goterm in geneTermList:
                    ne = [str(e[a]) for a in self.clusterData.domain]
                    if self.AddGOclass:
                        ne += [goterm]
                    nex = orange.Example(newdomain, ne)
                    for (id, var) in metas.items():
                        nex.setmeta(id, e.getmeta(id))
                    newdata.append( nex)

            if newdata.domain.classVar:
                self.send("Classified Examples", newdata)
            else:
                self.send("Classified Examples", None)
            self.send("Examples", newdata)
        else:
            self.send("Classified Examples", None)
            self.send("Examples", None)

    ##########################################################################
    # callback functions
    def browseRemember(self, lst, indx, loadMethod, dialogText, dialogTitle):
        if lst == []:
            startfile = "."
        else:
            startfile = lst[0]
        filename = QFileDialog.getOpenFileName(startfile, dialogText, None, dialogTitle)
        fn = str(filename)
        if fn in lst: # if already in list, remove it
            lst.remove(fn)
        lst.insert(0, fn)
        indx = 0
        loadMethod()

    def loadRemember(self, lst, filecombo, indx):
        loadedData = None
        if indx < len(lst):
            fn = lst[indx]
            if fn != "(none)":
                # remember the recent file list
                if fn in lst: # if already in list, remove it
                    lst.remove(fn)
                lst.insert(0, fn) # add to beginning of list
                self.setFilelist(filecombo, lst) # update combo
                loadedData = cPickle.load(open(fn, 'r'))
                indx = 0
        return loadedData, indx

    def browseAnnotation(self):
        self.browseRemember(self.RecentAnnotations, self.BAnnotationIndx, self.loadAnnotation, 'Annotation files (*.annotation)\nAll files(*.*)', 'Annotation Pickle File')
        self.BAnnotationIndx = 0

    def loadAnnotation(self):
        if DEBUG: print "loadAnnotation"
        self.annotation, self.BAnnotationIndx = self.loadRemember(self.RecentAnnotations, self.annotationCombo, self.BAnnotationIndx)
        fn = self.RecentAnnotations[0] ## the loaded one is (becomes) always moved to 0 position
        self.genesInAnnotationFile[fn] = self.annotation['gene2GOID'].keys() ## update, in case the file content changed
        self.updateEvidences()
        self.findTermsBuildDAG()

    def browseGOaspect(self):
        self.browseRemember(self.RecentGOaspects, self.BGOaspectIndx, self.loadGOaspect, 'GO files (*.go)\nAll files(*.*)', 'Gene Ontology Pickle File')
        self.BGOaspectIndx = 0

    def loadGOaspect(self, forced=1):
        if DEBUG: print "loadGOaspect"
        ## load if forced, or if index has changed
        ## if forced = 0 and index has not changed (still is 0) then don't reload the annotation data
        if forced == 1 or self.BGOaspectIndx <> 0 or self.GO == None:
            print "1:", str(self.RecentGOaspects) + "," + str(self.BGOaspectIndx)
            self.GO, self.BGOaspectIndx = self.loadRemember(self.RecentGOaspects, self.GOaspectCombo, self.BGOaspectIndx)
            print "2:", str(self.RecentGOaspects) + "," + str(self.BGOaspectIndx)
            self.updateEvidences()
            self.findTermsBuildDAG()

    def updateEvidences(self):
        if not(self.annotation) or not(self.GO): ## if data missing, just disable everything
            for (etype, tmpCB) in self.evidenceCheckBoxes.items():
                tmpCB.setText(etype)
                tmpCB.setEnabled(False)
            return

        # count the number of evidence in each type and number of genes with evidence; update the checkboxes
        evidenceTypeCn = {}
        for (gene, geneAnns) in self.annotation['gene2GOID'].items():
            for (daGOID, daNOT, daEvidence, daAspect, daDB_Object_Type) in geneAnns:
                if daAspect <> self.GO['aspect']: continue # skip annotations that are not for the loaded aspect
                (cn, lst) = evidenceTypeCn.get(daEvidence, (0, []))
                if gene not in lst:
                    lst = lst + [gene]
                evidenceTypeCn[daEvidence] = (cn + 1, lst)

        for (etype, tmpCB) in self.evidenceCheckBoxes.items():
            eCnLst = evidenceTypeCn.get(etype, None)
            if eCnLst:
                cn, lst = eCnLst
                tmpCB.setEnabled(True)
                tmpCB.setText('%s: %d annots (%d genes)' % (etype, cn, len(lst)))
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

    ##########################################################################
    # GO DAG calculations and filtering
    def runFilters(self):
        self.significantGOIDs = [] ## significant GOID to display in GO
        for (p, x, GOID) in self.GOIDsFound:
            if self.FilterPvalEnabled and p > self.FilterPvalue:
                break ## end of significant GO terms reached
            if self.FilterNumEnabled and x < self.FilterNumValue:
                continue ## not worth mentioning
            self.significantGOIDs.append( GOID)

        self.dag = GOlib.createGODAGtoDisplay(self.GO, self.GOtermValues.keys(), self.significantGOIDs)
        if self.FilterDepthEnabled:
            self.dag = GOlib.DAGfilterForDepth(self.dag, 'root', self.FilterDepthValue)
        self.updateDAG()

    def findTermsBuildDAG(self):
        self.dag = {}
        if DEBUG: print "findTermsBuildDAG, self.annotation: " + str(self.annotation <> None)
        if DEBUG: print "findTermsBuildDAG, self.GO: " + str(self.GO <> None)
        if self.annotation <> None and self.GO <> None:
            self.progressBarInit()
            evidences = [etype for (etype, tmpCB) in self.evidenceCheckBoxes.items() if tmpCB.isChecked()]
            if self.ReferenceType == 0: # from annotation
                ## for reference use the whole genome
                self.GOIDsFound, self.GOtermValues, clusterSet, referenceSet = GOlib.findTerms(self.annotation, self.GO, self.clusterGenes, None, evidences, self.progressBarSet, 0.0, 75.0)
            else: # from the given set of genes - received by signal
                ## for reference use genes in the reference list
                self.GOIDsFound, self.GOtermValues, clusterSet, referenceSet = GOlib.findTerms(self.annotation, self.GO, self.clusterGenes, self.referenceGenes, evidences, self.progressBarSet, 0.0, 75.0)
            n = len(clusterSet); N = len(referenceSet) # needed if relative frequencies need to be displayed
##            print n, N

            ## find the max number of cluster gene istances in a GO term
            maxNumIstances = max( [1] + [x for (GOterm, x, G, pval, genesInGOID, genesInGOIDdirect) in self.GOtermValues.values()])
            ## create a DAG with all the possible nodes
            ## and find the max depth of the DAG
            sigGOIDs = [goid for (_, _, goid) in self.GOIDsFound]
            tmpdag = GOlib.createGODAGtoDisplay(self.GO, self.GOtermValues.keys(), sigGOIDs)
            maxDepth = GOlib.DAGdepth(tmpdag)
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

            if pval: pval = "%1.4f" % pval
            vals = [GOterm, x, G, pval, genesInGOIDstr]
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

            if pval: pval = "%1.4f" % pval
            vals = [GOterm, x, G, pval, genesInGOIDstr]
            for j in range(len(vals)):
                self.sigTermsTable.setText(i, j, str(vals[j]))

if __name__=="__main__":
    import orange
    a = QApplication(sys.argv)
    ow = OWGOTermFinder()
    a.setMainWidget(ow)

##    d = orange.ExampleTable('testClusterSet.tab', dontCheckStored=1)
##    d = orange.ExampleTable('hjSmall.tab', dontCheckStored=1)
    d = orange.ExampleTable('hj.tab', dontCheckStored=1)
    ow.clusterDataset(d, 0)
    ow.show()
    a.exec_loop()
    ow.saveSettings()
