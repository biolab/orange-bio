"""
<name>GSEA</name>
<description>Gene set enrichment analysis.</description>
<contact>Marko Toplak (marko.toplak(@at@)gmail.com)</contact>
<priority>2025</priority>
<icon>icons/GSEA.svg</icon>
"""

from __future__ import absolute_import, with_statement
import os

from exceptions import Exception
import cPickle as pickle
from collections import defaultdict

from Orange.orng import orngServerFiles
from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import obiGene, obiGeneSets, obiGsea, obiKEGG

NAME = "GSEA"
DESCRIPTION = "Gene set enrichment analysis."
ICON = "icons/GSEA.svg"
PRIORITY = 2025

INPUTS = [("Examples", Orange.data.Table, "setData")]
OUTPUTS = [("Examples with selected genes only", Orange.data.Table),
           ("Results", Orange.data.Table),
           ("Distance Matrix", Orange.core.SymMatrix)]

REPLACES = ["_bioinformatics.widgets.OWGsea.OWGsea"]


def nth(l, n):
    return [ a[n] for a in l ]

def clearListView(lw):
    lw.clear()

def ne(a):
    return a if a != None else ""

def selectGenes(data, positions, how):
    """ Select genes on given positions.
    Parameter how specifies whether
    examples or attributes should be selected. """
    if how == "attributes":
        newatts = [ data.domain.attributes[i] for i in positions ]
        if data.domain.classVar:
            domain = orange.Domain(newatts, data.domain.classVar)
        else:
            domain = orange.Domain(newatts, False)
        domain.addmetas(data.domain.getmetas()) 
        return orange.ExampleTable(domain, data)
    elif how == "examples":
        return orange.ExampleTable(data.domain, [data[i] for i in positions ])

def comboboxItems(combobox, newitems):
    combobox.clear()
    if newitems:
        combobox.insertItems(0, newitems)
        #combobox.setCurrentItem(i)

def exportDistanceMatrix(resl):
    """
    Input: results as a list of tuples
    """
    dm = orange.SymMatrix(len(resl))

    for i in range(len(resl)-1):
        for j in range(i+1, len(resl)):
            gen1 = set(resl[i][1]['genes'])
            gen2 = set(resl[j][1]['genes'])
            dm[i,j] = float(len(gen1 & gen2)) / len(gen1 | gen2)

    return dm

def exportET(resl):
    #do not sort them inside
    
    if len(resl) <= 0:
        return None

    def collectionn(gs):
        return ",".join(gs.hierarchy) if gs.hierarchy else ""

    allCollections = sorted(set(collectionn(gs) for gs,_ in resl))

    vars = []
    vars.append(orange.StringVariable("Name"))
    vars.append(orange.EnumVariable("Collection", values=map(str,allCollections )))
    vars.append(orange.FloatVariable("NES"))
    vars.append(orange.FloatVariable("ES"))
    vars.append(orange.FloatVariable("P-value"))
    vars.append(orange.FloatVariable("FDR"))
    vars.append(orange.StringVariable("Geneset size"))
    vars.append(orange.StringVariable("Matched size"))
    vars.append(orange.StringVariable("Genes"))

    domain = orange.Domain(vars, False)

    examples = []
    for gs, dicr in resl:
        examples.append([str(ne(gs.id) + " " + ne(gs.name)), str(ne(collectionn(gs))), dicr['nes'], 
        dicr['es'], dicr['p'], min(dicr['fdr'],1.0), str(dicr['size']), str(dicr['matched_size']),  ", ".join(dicr['genes'])])

    return orange.ExampleTable(domain, examples)




class PhenotypesSelection(QGroupBox):
    """
    Window indices:
    0 - left chooser
    1 - right chooser

    wishedState: 0 not choosen anywhere, 1 choosen in left, 2 in right
    """

    def __init__(self, parent):
        QObject.__init__(self)
        grid = OWGUI.widgetBox(parent, "", orientation = "horizontal")
        grid.setMinimumWidth(250)
        grid.setMinimumHeight(100)

        self.boxes = [ OWGUI.listBox(grid, self) for a in range(2) ]

        for box in self.boxes: 
            #box.setSelectionMode(QListWidget.SingleSelection)
            box.setSelectionMode(QListWidget.MultiSelection)

        self.connect(self.boxes[0], SIGNAL("itemSelectionChanged ()"), self.highlighted1)
        self.connect(self.boxes[1], SIGNAL("itemSelectionChanged ()"), self.highlighted2)

        self.classes = []

        def createSquarePixmap(color = Qt.black):
            return OWGUI.createAttributePixmap("", color)

        self.whiteSq = createSquarePixmap(Qt.white)
        self.marked = [ createSquarePixmap(Qt.red), createSquarePixmap(Qt.blue) ]

        self.classVals = []

    def selectWanted(self):

        #prevent selection events when chenging here
        self.disableNot = True

        """
        Changes have to be calculated. Apply only changes because of potential
        troubles with flickering.
        """

        def disable(n, i):
            self.boxes[n].item(i).setIcon(self.whiteSq)
            if self.boxes[n].item(i) in self.boxes[n].selectedItems():
                modind = self.boxes[n].model().index(i, 0)
                self.boxes[n].selectionModel().select(modind, QItemSelectionModel.Deselect)

        def enable(n, i):
            self.boxes[n].item(i).setIcon(self.marked[n])
            if self.boxes[n].item(i) not in self.boxes[n].selectedItems():
                modind = self.boxes[n].model().index(i, 0)
                self.boxes[n].selectionModel().select(modind, QItemSelectionModel.Select)

        for boxi in range(2):

            toDisable = [ i for i in range(len(self.classVals)) \
                if self.wishedState[i] != boxi+1 ]

            for i in toDisable:
                disable(boxi, i)

            #enable every not choosen one that is wished
            toEnable = [ i for i in range(len(self.classVals)) \
                if self.wishedState[i] == boxi+1 ]

            for i in toEnable:
                enable(boxi, i)

        #allow selection events
        self.disableNot = False

        #print self.getSelection(), self.wishedState

    def highlighted(self, n):
        """
        Clicked on a i-th item of box n
        """

        selected = [ self.boxes[n].row(a) for a in self.boxes[n].selectedItems() ]

        if self.disableNot:
            return

        for i in range(len(self.classVals)):
            #print i, selected
            if i in selected:
                self.wishedState[i] = n+1 
            elif self.wishedState[i] == n+1:
                self.wishedState[i] = 0

        self.selectWanted()

    def highlighted1(self): return self.highlighted(0)
    def highlighted2(self): return self.highlighted(1)

    def setClasses(self, input, s1=0, s2=1):

        if input is not None:
            self.classVals = sorted(input)
            self.wishedState = [ 0 ] * len(self.classVals)

            self.wishedState[s1] = 1
            self.wishedState[s2] = 2

            self.setupBoxes()
            self.selectWanted()
        else:
            self.classVals = []
            self.setupBoxes()
            self.selectWanted()
 
    def getSelection(self):
        sels = [ [ self.classVals[i] for i,a in enumerate(self.wishedState) if a == n+1 ]
            for n in range(2) ]
        return sels

    def setupBoxes(self):
        for box in self.boxes:
            self.setupBox(box)

    def setupBox(self, box):
        # clear and fill box

        box.clear()
        for i,cv in enumerate(self.classVals):
            box.insertItem(i, cv)
            box.item(i).setIcon(self.whiteSq)

        if not self.classVals:
            box.setDisabled(True)
        else:
            box.setDisabled(False)

class OWGsea(OWWidget):
    settingsList = ["name", 
                    "perms", 
                    "minSubsetSize", 
                    "minSubsetSizeC", 
                    "maxSubsetSize", 
                    "maxSubsetSizeC", 
                    "minSubsetPart", 
                    "minSubsetPartC", 
                    "ptype", 
                    "loadFileName",
                    "gridSels",
                    "csgm",
                    "gsgo",
                    "gskegg",
                    "buildDistances",
                    "organismIndex",
                    "atLeast"]

    def UpdateOrganismComboBox(self):
        try:
            genome = obiKEGG.KEGGGenome()

            self.allOrganismCodes = genome

            essential = genome.essential_organisms()

            local = [name.split(".")[0].split("_")[-1]
                     for name in orngServerFiles.listfiles("KEGG")
                     if "kegg_genes" in name]

            entry_keys = map(genome.org_code_to_entry_key,
                             essential + local)

            entries = genome.batch_get(entry_keys)

            items = [entry.definition for entry in entries]

            self.organismTaxids = [entry.taxid for entry in entries]

            self.organismComboBox.clear()
            self.organismComboBox.addItems(items)
        finally:
            self.setBlocking(False)


    def __init__(self, parent=None, signalManager = None, name='GSEA'):
        OWWidget.__init__(self, parent, signalManager, name)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples with selected genes only", ExampleTable), ("Results", ExampleTable), ("Distance Matrix", orange.SymMatrix) ]

        self.res = None
        self.dm = None
        
        self.name = 'GSEA'
        self.minSubsetSize = 3
        self.minSubsetSizeC = True
        self.maxSubsetSize = 1000
        self.maxSubsetSizeC = True
        self.minSubsetPart = 10
        self.minSubsetPartC = True
        self.perms = 100
        self.csgm = False
        self.gsgo = False
        self.gskegg = False
        self.buildDistances = False
        self.selectedPhenVar = 0
        self.organismIndex = 0
        self.atLeast = 3

        self.permutationTypes =  [("Phenotype", "p"),("Gene", "g") ]
        self.ptype = 0

        self.correlationTypes = [ ("Signal2Noise", "s2n") ]
        self.ctype = 0
        
        self.data = None
        self.geneSets = {}

        self.tabs = OWGUI.tabWidget(self.controlArea)

        ca = OWGUI.createTabPage(self.tabs, "Basic")

        box = OWGUI.widgetBox(ca, 'Organism')
        #FROM KEGG WIDGET - organism selection
        self.allOrganismCodes = {}
        self.organismTaxids = []
        self.organismComboBox = cb = OWGUI.comboBox(box, self, "organismIndex", items=[], debuggingEnabled=0) #changed
        cb.setMaximumWidth(200)

        #OWGUI.checkBox(box, self, "csgm", "Case sensitive gene matching")

        box2 = OWGUI.widgetBox(ca, "Descriptors")
        self.phenCombo = OWGUI.comboBox(box2, self, "selectedPhenVar", items=[], callback=self.phenComboChange, label="Phenotype:")
        self.geneCombo = OWGUI.comboBox(box2, self, "selectedGeneVar", items=[], label = "Gene:")

        self.allowComboChangeCallback = False

        ma = self.mainArea

        self.listView = QTreeWidget(ma)
        ma.layout().addWidget(self.listView)
        self.listView.setAllColumnsShowFocus(1)
        self.listView.setColumnCount(9)
        self.listView.setHeaderLabels(["Collection", "Geneset", "NES", "ES", "P-value", "FDR", "Size", "Matched Size", "Genes"])
        
        self.listView.header().setStretchLastSection(True)
        self.listView.header().setClickable(True)
        self.listView.header().setSortIndicatorShown(True)
        self.listView.setSortingEnabled(True)
        #self.listView.header().setResizeMode(0, QHeaderView.Stretch)
        
        self.listView.setSelectionMode(QAbstractItemView.NoSelection)
        self.connect(self.listView, SIGNAL("itemSelectionChanged()"), self.newPathwaySelected)

        OWGUI.separator(ca)

        OWGUI.widgetLabel(ca, "Phenotype selection:")
        self.psel = PhenotypesSelection(ca)
        
        self.resize(600,50)
 
        OWGUI.separator(ca)

        OWGUI.checkBox(ca, self, "buildDistances", "Compute geneset distances")

        self.btnApply = OWGUI.button(ca, self, "&Compute", callback = self.compute, disabled=0)
        
        fileBox = OWGUI.widgetBox(ca, orientation='horizontal')
        OWGUI.button(fileBox, self, "Load", callback = self.loadData, disabled=0, debuggingEnabled=0)
        OWGUI.button(fileBox, self, "Save", callback = self.saveData, disabled=0, debuggingEnabled=0)
 
        #ca.layout().addStretch(1)

        ca = OWGUI.createTabPage(self.tabs, "Gene sets")
        
        box = OWGUI.widgetBox(ca)

        self.gridSel = []
        self.geneSel = []  #FIXME temporary disabled - use the same as in new "David" widget
        self.lbgs = OWGUI.listBox(box, self, "gridSel", "geneSel", selectionMode = QListWidget.MultiSelection)
        OWGUI.button(box, self, "From &File", callback = self.addCollection, disabled=0, debuggingEnabled=0)

        box = OWGUI.widgetBox(box, "Additional sources:")
        OWGUI.checkBox(box, self, "gskegg", "KEGG pathways")
        OWGUI.checkBox(box, self, "gsgo", "GO terms")
 
        #ca.layout().addStretch(1)

        ca = OWGUI.createTabPage(self.tabs, "Settings")
        box = OWGUI.widgetBox(ca, 'Properties')

        self.permTypeF = OWGUI.comboBoxWithCaption(box, self, "ptype", items=nth(self.permutationTypes, 0), \
            tooltip="Permutation type.", label="Permute")
        _ = OWGUI.spin(box, self, "perms", 50, 1000, orientation="horizontal", label="Times")
        self.corTypeF = OWGUI.comboBoxWithCaption(box, self, "ctype", items=nth(self.correlationTypes, 0), \
            tooltip="Correlation type.", label="Correlation")


        box = OWGUI.widgetBox(ca, 'Subset Filtering')

        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Size", 1, 10000, "minSubsetSizeC", "minSubsetSize", "") #TODO check sizes
        _,_ = OWGUI.checkWithSpin(box, self, "Max. Subset Size", 1, 10000, "maxSubsetSizeC", "maxSubsetSize", "")
        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Part (%)", 1, 100, "minSubsetPartC", "minSubsetPart", "")

        box = OWGUI.widgetBox(ca, 'Gene Filtering')

        _ = OWGUI.spin(box, self, "atLeast", 2, 10, label="Min. Values in Group")

        ca.layout().addStretch(1)

        self.addComment("Computation was not started.")

        if sys.platform == "darwin":
            self.loadFileName = os.path.expanduser("~/")
        else:
            self.loadFileName = "."

        self.gridSels = []
        self.loadSettings()

        self.setBlocking(True)
        QTimer.singleShot(0, self.UpdateOrganismComboBox)

        def cleanInvalid(maxn):
            """
            Removes invalid gene set selection
            """
            notAllOk = True

            while notAllOk:
                self.gridSels = getattr(self, "gridSels")
                notAllOk = False
                for i,a in enumerate(self.gridSels):
                    if a >= maxn:
                        self.gridSels.pop(i)
                        notAllOk = True
                        break
        
        cleanInvalid(len(self.geneSel))

        self.gridSel = self.gridSels
        self.gridSels = self.gridSel

    def addCollection(self):
        fname = self.chooseGeneSetsFile()

        if fname:
            if fname not in self.geneSel:
        
                #add it to the list, choose it and keep
                #all old choosen
                gridSelc = list(self.gridSel)

                self.geneSel.append(fname)
                self.geneSel = self.geneSel

                gridSelc.append(len(self.geneSel)-1)
            
                self.gridSel = gridSelc
                self.gridSels = self.gridSel #for saving


    def saveData(self):
        self.warning('')
        
        if self.res != None:
            filename = QFileDialog.getSaveFileName(self, 'Save GSEA data', '', 'GSEA files (*.gsea)')
            if filename:
                fn = ""
                head, tail = os.path.splitext(str(filename))
                if not tail:
                    fn = head + ".gsea"
                else:
                    fn = str(filename)
                    
                fp = open(fn, "wb" )
                pickle.dump(self.res, fp, -1)
                pickle.dump(self.dm, fp, -1)
                fp.close()
        else:
            self.warning('No internal data to save.')
    
    def loadData(self):                
        self.loadFileName = str(QFileDialog.getOpenFileName(self, 'Open GSEA data', self.loadFileName, "GSEA files (*.gsea)"))
        if self.loadFileName == "": 
            if sys.platform == "darwin":
                self.loadFileName = os.path.expanduser("~/")
            else:
                self.loadFileName = "."
            return
        
        fp = open(self.loadFileName, "rb")
        res = pickle.load(fp)
        
        try:
            dm = pickle.load(fp)
        except:
            dm = None
        
        fp.close()
        
        self.compute(res, dm)

    def newPathwaySelected(self):
        #print "newPathwaySelected"
        qApp.processEvents()

        if not self.selectable:
            return

        if self.res == None:
            return

        outat = set([])
        for item in self.listView.selectedItems():
            iname = self.lwiToGeneset[item]
            outat.update(self.res[iname]['genes'])

        #print "OUTGENES",  outat

        select = sorted(set(reduce(lambda x,y: x | y, 
            [ set(self.find_genes_dic[name]) for name in outat ],
            set())))

        #print "SELECT", select
            
        dataOut = selectGenes(self.data, select, self.find_genes_loc)
        self.send("Examples with selected genes only", dataOut)

    def resultsOut(self, data):
        self.send("Results", data)

    def genesetDistOut(self, dm):
        self.send("Distance Matrix", dm)


    def fillResults(self, res):
        clearListView(self.listView)

        self.lwiToGeneset = {}

        def writeGenes(g):
            return ", ".join(g)

        for gs, rdic in res.items():
            collection = ",".join(gs.hierarchy) if gs.hierarchy else ""
            name = ne(gs.id) + " " + ne(gs.name)
            item = QTreeWidgetItem(self.listView)
            item.setText(0, collection)
            item.setText(1, name)
            item.setText(2, "%0.3f" % rdic['nes'])
            item.setText(3, "%0.3f" % rdic['es'])
            item.setText(4, "%0.3f" % rdic['p'])
            item.setText(5, "%0.3f" % min(rdic['fdr'],1.0))
            item.setText(6, str(rdic['size']))
            item.setText(7, str(rdic['matched_size']))
            item.setText(8, writeGenes(rdic['genes']))

            self.lwiToGeneset[item] = gs

    def addComment(self, comm):
        item = QTreeWidgetItem(self.listView)
        item.setText(0, comm)   

    def setSelMode(self, bool):
        if bool:
            self.selectable = True
            self.listView.setSelectionMode(QAbstractItemView.ExtendedSelection)
        else:
            self.selectable = False
            self.listView.setSelectionMode(QListView.NoSelection)

    def compute(self, res=None, dm=None):

        collectionNames = [ self.geneSel[a] for a in self.gridSel ]

        organism = self.organismTaxids[self.organismIndex]

        if self.gsgo:
            collectionNames.append((("GO",),organism))
        if self.gskegg:
            collectionNames.append((("KEGG",),organism))

        self.geneSets = obiGeneSets.collections(*collectionNames)

        self.resultsOut(None)

        qApp.processEvents()
        self.res = res
        self.dm = dm

        clearListView(self.listView)
        self.addComment("Computing...")
        qApp.processEvents()

        self.phenVar = self.phenCands[self.selectedPhenVar][0]
        self.geneVar = self.geneCands[self.selectedGeneVar]

        if self.res == None and self.data:
            self.setSelMode(False)

            pb = OWGUI.ProgressBar(self, iterations=self.perms+2)

            if hasattr(self, "btnApply"):
                self.btnApply.setFocus()

            kwargs = {}
            dkwargs = {}

            dkwargs["phenVar"] = self.phenVar
            dkwargs["geneVar"] = self.geneVar

            if not obiGsea.already_have_correlations(self.data):

                selectedClasses = self.psel.getSelection()
                fc = "Phenotype group empty. Stopped."
                if len(selectedClasses[0]) == 0:
                    self.addComment(fc)
                    return
                elif len(selectedClasses[1]) == 0:
                    self.addComment(fc)
                    return

                dkwargs["classValues"] = selectedClasses

                dkwargs["atLeast"] = self.atLeast

                permtype = self.permutationTypes[self.ptype][1]
                kwargs["permutation"] = "class" if permtype == "p" else "genes"

            def ifr(case, t, f):
                if case: return t
                else: return f

            kwargs["minSize"] = \
                ifr(self.minSubsetSizeC, self.minSubsetSize, 1)
            kwargs["maxSize"] = \
                ifr(self.maxSubsetSizeC, self.maxSubsetSize, 1000000)
            kwargs["minPart"] = \
                ifr(self.minSubsetPartC, self.minSubsetPart/100.0, 0.0)


            #create gene matcher
            genematcher = obiGene.matcher([[obiGene.GMKEGG(organism)] + ([obiGene.GMDicty()] if organism == "352472"  else [])])

            #dkwargs["caseSensitive"] = self.csgm

            gso = obiGsea.GSEA(self.data, matcher=genematcher, **dkwargs)

            
            for gs in self.geneSets:
                gso.addGenesets([gs])
                qApp.processEvents()

            self.res = gso.compute(n=self.perms, callback=pb.advance, **kwargs)
            
            pb.finish()
            
        if self.res != None:
            if len(self.res) > 0:
                self.fillResults(self.res)
                self.setSelMode(True)
                resl = self.res.items()

                etres = exportET(resl)
                self.resultsOut(etres)

                self.find_genes_dic, self.find_genes_loc = \
                    find_genes_dic(self.res, self.data, self.geneVar)

                if self.buildDistances:
                    if self.dm == None:
                        self.dm = exportDistanceMatrix(resl)
                        
                        for ex in etres:
                            ex.name = str(ex[0])
                        self.dm.setattr("items", etres)
                else:
                    self.dm = None

            else:
                self.setSelMode(False)
                clearListView(self.listView)
                self.dm = None
                self.addComment("No genesets found.")

            self.genesetDistOut(self.dm)


    def setData(self, data):

        self.allowComboChangeCallback = False

        self.data = data

        if data:
            taxid = data_hints.get_hint(data, "taxid", None)
            try:
                self.organismIndex = self.organismTaxids.index(taxid)
            except Exception, ex:
                pass

            if obiGsea.already_have_correlations(data):

                #disable correlation type
                comboboxItems(self.corTypeF, [])
                self.corTypeF.setDisabled(True)
                #set permutation type to fixed
                self.ptype = 1
                self.permTypeF.setDisabled(True)

            else:
                #enable correlation type
                comboboxItems(self.corTypeF, nth(self.correlationTypes, 0))
                self.corTypeF.setDisabled(False)
                #allow change of permutation type
                self.permTypeF.setDisabled(False)
                #print "set classes"

            self.phenCombo.clear()
            self.phenCands = obiGsea.phenotype_cands(data)
            self.phenCombo.addItems(map(lambda x: str(x[0]), self.phenCands))

            self.phenCombo.setDisabled(len(self.phenCands) <= 1)

            self.selectedPhenVar = 0

            self.allowComboChangeCallback = True

            self.phenComboChange()



    def phenComboChange(self):

        if self.allowComboChangeCallback == False:
            return

        pv = self.phenCands[self.selectedPhenVar]

        self.psel.setClasses(pv[1] if len(pv[1]) > 0 else None)

        def conv(x):
            if x == True:
                return "Attribute names"
            else:
                return str(x)

        self.geneCombo.clear()

        look_cols = obiGsea.is_variable(pv[0])

        if obiGsea.already_have_correlations(self.data):
            look_cols = not obiGsea.need_to_transpose_single(self.data)

        self.geneCands = obiGsea.gene_cands(self.data, look_cols)
        self.geneCombo.addItems(map(lambda x: conv(x), self.geneCands))
        self.selectedGeneVar = 0

        self.geneCombo.setDisabled(len(self.geneCands) <= 1)

    def chooseGeneSetsFile(self):
        """
        Return choosen gene sets file name or None, if no file
        was choosen.
        """
        filename = str(QFileDialog.getOpenFileName(self,  "Choose gene set collection", './', "Gene Collections (*.gmt)"))
        return filename

def find_genes_dic(res, data, geneVar):
    """
    Builds a dictionary of where to find the genes to
    avoid delay with choosing them.

    If geneVar is True or it is a group parameter, then
    use columns as genes: out of genes select only
    those contained in selected pathways.
    """

    def rorq(a, name):
        """ Group annatation or question mark. """
        try: 
            return a.attributes[name]
        except: 
            return '?'

    d = defaultdict(list)
    if geneVar == True or geneVar in obiGsea.allgroups(data):
        for i,a in enumerate(data.domain.attributes):
            if geneVar == True:
                d[a.name].append(i)
            else:
                d[rorq(a, geneVar)].append(i)
        return d, "attributes"
    else:
        for i,ex in enumerate(data):
            d[ex[geneVar].value].append(i)
        return d, "examples"

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWGsea()
    ow.show()

    #d = orange.ExampleTable('/home/marko/testData.tab')
    d = orange.ExampleTable('/home/marko/orange/add-ons/Bioinformatics/sterolTalkHepa.tab')
    #d = orange.ExampleTable('/home/marko/ddd.tab')
    #d = orange.ExampleTable('tmp.tab')
    #d = orange.ExampleTable('../gene_three_lines_log.tab')

    QTimer.singleShot(1000, lambda : ow.setData(d))

    a.exec_()
    ow.saveSettings()
