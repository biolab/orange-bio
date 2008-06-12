"""
<name>GSEA</name>
<description>Gene Set Enrichment Analysis</description>
<contact>Marko Toplak (marko.toplak(@at@)gmail.com)</contact>
<priority>210</priority>
<icon>icons/GSEA.png</icon>
"""

from OWWidget import *
import OWGUI
import obiGsea
from exceptions import Exception

def nth(l, n):
    return [ a[n] for a in l ]

def clearListView(lw):
    it = lw.firstChild()
    while it:
        lw.takeItem(it)
        it = lw.firstChild()

def dataWithAttrs(data, attributes):
    attributes = dict([(a,1) for a in attributes])
    newatts = [ a for a in data.domain.attributes if a.name in attributes ]
    if data.domain.classVar:
        domain = orange.Domain(newatts, data.domain.classVar)
    else:
        domain = orange.Domain(newatts, False)
    return orange.ExampleTable(domain, data)

def comboboxItems(combobox, newitems):
    combobox.clear()
    if newitems:
        combobox.insertStrList(newitems)
        #combobox.setCurrentItem(i)

def getClasses(data):
    return [ str(a) for a in data.domain.classVar ]

class PhenotypesSelection():
    """
    Window indices:
    0 - left chooser
    1 - right chooser

    wishedState: 0 not choosen anywhere, 1 choosen in left, 2 in right
    """

    def __init__(self, parent):
        grid = QHBox(parent)
        grid.setMinimumWidth(250)
        grid.setMinimumHeight(100)

        self.boxes = [ QListBox(grid), QListBox(grid) ]

        for box in self.boxes:
            box.setSelectionMode(QListBox.Single)

        QObject.connect(self.boxes[0], SIGNAL("selected ( int )"), self.highlighted1)
        QObject.connect(self.boxes[1], SIGNAL("selected ( int )"), self.highlighted2)

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
            self.boxes[n].changeItem(self.whiteSq, self.classVals[i], i)
            self.choosen[n][i] = False

        def enable(n, i):
            self.boxes[n].changeItem(self.marked[n], self.classVals[i], i)
            self.choosen[n][i] = True

        self.boxes[0].setCurrentItem(self.lastSel[0])
        self.boxes[1].setCurrentItem(self.lastSel[1])

        #print self.wishedState, self.choosen

        for boxi in range(2):

            #disable every choosen one which is not wished any more
            toDisable = [ i for i,e in enumerate(self.choosen[boxi]) \
                    if e == True and self.wishedState[i] != boxi+1 ]

            for i in toDisable:
                disable(boxi, i)

            #enable every not choosen one that is wished
            toEnable = [ i for i,e in enumerate(self.choosen[boxi]) \
                    if e == False and self.wishedState[i] == boxi+1 ]

            for i in toEnable:
                enable(boxi, i)

        #allow selection events
        self.disableNot = False

        #print self.getSelection(), self.wishedState

    def highlighted(self, n, i):
        """
        Clicked on a i-th item of box n
        """

        if self.disableNot:
            return

        self.lastSel[n] = i

        if self.wishedState[i] == n+1:
            self.wishedState[i] = 0
        else:
            self.wishedState[i] = n+1

        self.selectWanted()

    def highlighted1(self, i): return self.highlighted(0, i)
    def highlighted2(self, i): return self.highlighted(1, i)

    def setClasses(self, input, s1=0, s2=1):

        self.classVals = sorted(input)
        self.wishedState = [ 0 ] * len(self.classVals)

        self.choosen = [ [ False ] * len(self.classVals), [ False ] * len(self.classVals) ]

        self.wishedState[s1] = 1
        self.wishedState[s2] = 2

        self.lastSel = [s1, s2]

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
        for cv in self.classVals:
            box.insertItem(self.whiteSq, cv)

        if not self.classVals:
            box.setDisabled(True)
        else:
            box.setDisabled(False)

class OWGsea(OWWidget):

    settingsList = [ "name", "perms", "minSubsetSize", "minSubsetSizeC", "maxSubsetSize", "maxSubsetSizeC", "minSubsetPart", "minSubsetPartC", "ptype", "gridSel" ]

    def __init__(self, parent=None, signalManager = None, name='GSEA'):
        OWWidget.__init__(self, parent, signalManager, name)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples with selected genes only", ExampleTable), ("Results", ExampleTable), ("Distance Matrix", orange.SymMatrix) ]

        self.res = None

        self.name = 'GSEA'
        self.minSubsetSize = 3
        self.minSubsetSizeC = True
        self.maxSubsetSize = 1000
        self.maxSubsetSizeC = True
        self.minSubsetPart = 10
        self.minSubsetPartC = True
        self.perms = 100

        self.permutationTypes =  [("Phenotype", "p"),("Gene", "g") ]
        self.ptype = 0

        self.organisms = [ ("hsa", "hsa"), ("ddi", "ddi") ]
        self.otype = 0

        self.correlationTypes = [ ("Signal2Noise", "s2n") ]
        self.ctype = 0

        self.loadSettings()

        self.data = None
        self.geneSets = {}

        ca = self.controlArea
        ca.setMaximumWidth(500)

        box = QVGroupBox(ca)
        box.setTitle('Organism')

        OWGUI.comboBox(box, self, "otype", \
            items=nth(self.organisms, 0), tooltip="Organism")

        OWGUI.separator(ca)

        box = QVGroupBox(ca)
        box.setTitle('Properties')

        self.permTypeF = OWGUI.comboBoxWithCaption(box, self, "ptype", items=nth(self.permutationTypes, 0), \
            tooltip="Permutation type.", label="Permutate")

        _ = OWGUI.spin(box, self, "perms", 50, 1000, orientation="horizontal", label="Times")

        self.corTypeF = OWGUI.comboBoxWithCaption(box, self, "ctype", items=nth(self.correlationTypes, 0), \
            tooltip="Correlation type.", label="Correlation")

        OWGUI.separator(ca)

        box = QVGroupBox(ca)
        box.setTitle('Subset Filtering')

        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Size", 1, 10000, "minSubsetSizeC", "minSubsetSize", "") #TODO check sizes
        _,_ = OWGUI.checkWithSpin(box, self, "Max. Subset Size", 1, 10000, "maxSubsetSizeC", "maxSubsetSize", "")
        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Part (%)", 1, 100, "minSubsetPartC", "minSubsetPart", "")

        OWGUI.separator(ca)

        box = QVGroupBox(ca)
        box.setTitle("Gene Sets")

        self.gridSel = []
        self.geneSel = [ a[0] for a in obiGsea.getCollectionFiles() ]
        self.lbgs = OWGUI.listBox(box, self, "gridSel", "geneSel", selectionMode = QListBox.Multi)
        #OWGUI.button(box, self, "From &File", callback = self.addCollection, disabled=0)

        ma = self.mainArea
        boxL = QVBoxLayout(ma, QVBoxLayout.TopToBottom)
        #box.setTitle("Results")

        self.listView = QListView(ma)
        for header in ["Geneset", "NES", "ES", "P-value", "FDR", "Size", "Matched Size", "Genes"]:
            self.listView.addColumn(header)
        self.listView.setSelectionMode(QListView.NoSelection)
        self.connect(self.listView, SIGNAL("selectionChanged ( QListViewItem * )"), self.newPathwaySelected)
        boxL.addWidget(self.listView)

        OWGUI.separator(ca)

        box = QVGroupBox(ca)
        box.setTitle("Phenotypes")

        self.psel = PhenotypesSelection(box)

        self.resize(600,50)
 
        OWGUI.separator(ca)
        self.btnApply = OWGUI.button(ca, self, "&Compute", callback = self.compute, disabled=0)

        #gen1 = getGenesets()
        #for name,genes in gen1.items():
        #    self.addGeneset(name, genes)

        self.addComment("Computation was not started.")

    def addCollection(self):
        fname = self.chooseGeneSetsFile()
        if fname:
            if fname not in self.geneSel:
                self.geneSel.append(fname)

    def newPathwaySelected(self, item):

        qApp.processEvents()

        if not self.selectable:
            return

        iname = self.lwiToGeneset[item]
        outat = self.res[iname][6]

        dataOut =  dataWithAttrs(self.data, outat)
        self.send("Examples with selected genes only", dataOut)

    def resultsOut(self, data):
        self.send("Results", data)

    def genesetDistOut(self, dm):
        self.send("Distance Matrix", dm)

    def exportET(self, resl):
        #do not sort them inside
        
        if len(resl) <= 0:
            return None

        vars = []
        vars.append(orange.StringVariable("Name"))
        vars.append(orange.FloatVariable("NES"))
        vars.append(orange.FloatVariable("ES"))
        vars.append(orange.FloatVariable("P-value"))
        vars.append(orange.FloatVariable("FDR"))
        vars.append(orange.StringVariable("Geneset size"))
        vars.append(orange.StringVariable("Matched size"))
        vars.append(orange.StringVariable("Genes"))
    
        domain = orange.Domain(vars, False)

        examples = []
        for name, (es, nes, pval, fdr, os, ts, genes) in resl:
            examples.append([name, nes, es, pval, min(fdr,1.0), str(os), str(ts),  ", ".join(genes)])

        return orange.ExampleTable(domain, examples)


    def exportDistanceMatrix(self, resl):
        """
        Input: results as a list of tuples
        """

        dm = orange.SymMatrix(len(resl))
    
        for i in range(len(resl)-1):
            for j in range(i+1, len(resl)):
                gen1 = set(resl[i][1][6])
                gen2 = set(resl[j][1][6])
                dm[i,j] = float(len(gen1 & gen2)) / len(gen1 | gen2)

        return dm


    def fillResults(self, res):

        clearListView(self.listView)

        self.lwiToGeneset = {}

        def writeGenes(g):
            return ", ".join(genes)

        for name, (es, nes, pval, fdr, os, ts, genes) in res.items():
            item = QListViewItem(self.listView)
            item.setText(0, name)
            item.setText(1, "%0.3f" % nes)
            item.setText(2, "%0.3f" % es)
            item.setText(3, "%0.3f" % pval)
            item.setText(4, "%0.3f" % min(fdr,1.0))
            item.setText(5, str(os))
            item.setText(6, str(ts))
            item.setText(7, writeGenes(genes))

            self.lwiToGeneset[item] = name

    def addComment(self, comm):
        item = QListViewItem(self.listView)
        item.setText(0, comm)

    def setSelMode(self, bool):
        if bool:
            self.selectable = True
            self.listView.setSelectionMode(QListView.Single)
        else:
            self.selectable = False
            self.listView.setSelectionMode(QListView.NoSelection)

    def compute(self):

        #self.lbgs.clear()

        #LOAD GENE SETS
        collectionNames = [ self.geneSel[a] for a in self.gridSel ]
        self.geneSets = obiGsea.collections(collectionNames, default=False)

        #self.geneSel.append("fffafda")

        clearListView(self.listView)
        self.addComment("Computing...")

        self.resultsOut(None)

        qApp.processEvents()

        if self.data:

            self.setSelMode(False)

            pb = OWGUI.ProgressBar(self, iterations=self.perms+2)

            if hasattr(self, "btnApply"):
                self.btnApply.setFocus()

            selectedClasses = self.psel.getSelection()
            fc = "Phenotype group empty. Stopped."
            if len(selectedClasses[0]) == 0:
                self.addComment(fc)
                return
            elif len(selectedClasses[1]) == 0:
                self.addComment(fc)
                return 

            kwargs = {}

            def ifr(case, t, f):
                if case: return t
                else: return f

            kwargs["minSize"] = \
                ifr(self.minSubsetSizeC, self.minSubsetSize, 1)
            kwargs["maxSize"] = \
                ifr(self.maxSubsetSizeC, self.maxSubsetSize, 1000000)
            kwargs["minPart"] = \
                ifr(self.minSubsetPartC, self.minSubsetPart/100.0, 0.0)

            if len(self.data) > 1:
                permtype = self.permutationTypes[self.ptype][1]
                kwargs["permutation"] = ifr(permtype == "p", "class", "genes") 

            dkwargs = {}
            if len(self.data) > 1:
                dkwargs["classValues"] = selectedClasses
 

            organism = self.organisms[self.otype][1]

            gso = obiGsea.GSEA(organism="ddi")  # ORIGINALLY HSA
            gso.setData(self.data, **dkwargs)

            for name,genes in self.geneSets.items():
                gso.addGeneset(name, genes)
                qApp.processEvents()

            self.res = gso.compute(n=self.perms, callback=pb.advance, **kwargs)
            
            pb.finish()

            if len(self.res) > 0:
                self.fillResults(self.res)
                self.setSelMode(True)
                resl = self.res.items()

                etres = self.exportET(resl)

                self.resultsOut(etres)
                dm = self.exportDistanceMatrix(resl)

                dm.setattr("items", etres)
                for ex in etres:
                    ex.name = str(ex[0])

                self.genesetDistOut(dm)

            else:
                self.setSelMode(False)
                clearListView(self.listView)
                self.addComment("No genesets found.")


    def setData(self, data):
        self.data = data

        if data:
            if len(data) == 1:
                #disable correlation type
                comboboxItems(self.corTypeF, [])
                self.corTypeF.setDisabled(True)
                #set permutation type to fixed
                self.permTypeF.setCurrentItem(1)
                self.permTypeF.setDisabled(True)
                
                self.psel.setClasses([])
            else:
                #enable correlation type
                comboboxItems(self.corTypeF, nth(self.correlationTypes, 0))
                self.corTypeF.setDisabled(False)
                #allow change of permutation type
                self.permTypeF.setDisabled(False)

                self.psel.setClasses(getClasses(data))

    def chooseGeneSetsFile(self):
        """
        Return choosen gene sets file name or None, if no file
        was choosen.
        """
        filename = str(QFileDialog.getOpenFileName("./","Gene Collections (*.gmt *.pck)", self, "open", "Choose gene set collection"))
        return filename


def unpckGS(filename):
    import pickle
    f = open(filename,'rb')
    return pickle.load(f)

def getGenesets():
    import orngRegistry
    #return unpckGS("../abc/precompPathways.pck")
    return unpckGS(orngRegistry.bufferDir + "/gsea/geneSets_Gsea_KEGGhsa.pck")
    

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWGsea()
    a.setMainWidget(ow)
    ow.show()

    #d = orange.ExampleTable('DLBCL_200a.tab')
    #d = orange.ExampleTable('brown-selected.tab')

    #d = orange.ExampleTable('testCorrelated.tab')
    #ow.setData(d)

    #d = orange.ExampleTable("sterolTalkHepa.tab")
    #ow.setData(d)

    #d = orange.ExampleTable("demo.tab")
    #ow.setData(d)

    #d = orange.ExampleTable("tmp.tab")
    #ow.setData(d)

    d = orange.ExampleTable("../abc/abc_gsea_1.tab")
    ow.setData(d)

    a.exec_loop()
    ow.saveSettings()
