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
import cPickle as pickle

def nth(l, n):
    return [ a[n] for a in l ]

def clearListView(lw):
    lw.clear()
    #it = lw.firstChild()
    #while it:
    #    lw.takeItem(it)
    #    it = lw.firstChild()

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
        combobox.insertItems(0, newitems)
        #combobox.setCurrentItem(i)

def getClasses(data):
    return [ a.value for a in data.domain.classVar ]

class PhenotypesSelection(QGroupBox):

    def __init__(self, parent, s1=0, s2=1):
        QObject.__init__(self)
        grid = OWGUI.widgetBox(parent, "", orientation = "horizontal")
        grid.setMinimumWidth(250)
        grid.setMinimumHeight(100)
        
        self.p1b = OWGUI.listBox(grid, self)
        self.p2b = OWGUI.listBox(grid, self)

        self.connect(self.p1b,  SIGNAL("currentRowChanged(int)"), self.highlighted1)
        self.connect(self.p2b,  SIGNAL("currentRowChanged(int)"), self.highlighted2)

        self.classes = []

        def createSquarePixmap(color = Qt.black):
            pixmap = QPixmap(13, 13)
            painter = QPainter()
            painter.begin(pixmap)
            painter.setPen(color);
            painter.setBrush(color);
            painter.drawRect(0, 0, 13, 13);
            painter.end()
            return pixmap

        self.whiteSq = QIcon(createSquarePixmap(Qt.white))
        self.redSq = QIcon(createSquarePixmap(Qt.red))
        self.blueSq = QIcon(createSquarePixmap(Qt.blue))

        self.classVals = []

        self.setStates(s1, s2)

    def setStates(self, s1 = 0, s2 = 1):
        self.state1 = self.ls1 = s1
        self.state2 = self.ls2 = s2

        if self.state1 == self.state2:
            if self.state1 == 0: 
                self.state2 = 1
            else: 
                self.state2 = 0

        self.selectWanted()

    def selectWanted(self):
        self.disableNot = True

        try:
            self.p1b.item(self.ls1).setIcon(self.whiteSq)
            self.p2b.item(self.ls2).setIcon(self.whiteSq)
        except:
            #except can happen only if both are illegal
            pass

        try:
            self.p1b.setCurrentRow(self.state1)
            self.p2b.setCurrentRow(self.state2)
            self.p1b.currentItem().setIcon(self.redSq)
            self.p2b.currentItem().setIcon(self.blueSq)
            self.ls1 = self.state1
            self.ls2 = self.state2
        except:
            pass

        self.disableNot = False

    def highlighted1(self, i):
        if self.disableNot:
            return
        if i == self.state2:
            self.state2 = self.state1
        self.state1 = i
        self.selectWanted()

    def highlighted2(self, i):
        if self.disableNot:
            return
        if i == self.state1:
            self.state1 = self.state2
        self.state2 = i
        self.selectWanted()

    def setClasses(self, input, s1=0, s2=1):
        self.classVals = sorted(input)
        self.setupBoxes()
        self.setStates(s1, s2)

    def getSelection(self):
        return (self.classVals[self.state1], self.classVals[self.state2])

    def setupBoxes(self):
        self.setupBox(self.p1b)
        self.setupBox(self.p2b)

    def setupBox(self, box):
        box.clear()
        for cv in self.classVals:
            box.addItem(QListWidgetItem(self.whiteSq, cv))
            
        if not self.classVals:
            box.setDisabled(True)
        else:
            box.setDisabled(False)

class OWGsea(OWWidget):

    settingsList = [ "name", "perms", "minSubsetSize", "minSubsetSizeC", "maxSubsetSize", "maxSubsetSizeC", \
        "minSubsetPart", "minSubsetPartC", "ptype" ]

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

        self.permutationTypes =  [("Phenotype", "p"),("Gene", "g") ]
        self.ptype = 0

        self.correlationTypes = [ ("Signal2Noise", "s2n") ]
        self.ctype = 0

        #self.loadSettings()
        self.data = None
        self.geneSets = {}

        ca = self.controlArea
        ca.setMaximumWidth(500)

        box = OWGUI.widgetBox(ca, 'Permutate')

        self.permTypeF = OWGUI.comboBox(box, self, "ptype", items=nth(self.permutationTypes, 0), \
            tooltip="Permutation type.")

        _ = OWGUI.spin(box, self, "perms", 50, 1000, orientation="horizontal", label="Times")

        OWGUI.separator(ca)

        box = OWGUI.widgetBox(ca, 'Correlation Calculation')

        self.corTypeF = OWGUI.comboBox(box, self, "ctype", items=nth(self.correlationTypes, 0), \
            tooltip="Correlation type.")

        OWGUI.separator(ca)

        box = OWGUI.widgetBox(ca, 'Subset Filtering')

        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Size", 1, 10000, "minSubsetSizeC", "minSubsetSize", "") #TODO check sizes
        _,_ = OWGUI.checkWithSpin(box, self, "Max. Subset Size", 1, 10000, "maxSubsetSizeC", "maxSubsetSize", "")
        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Part (%)", 1, 100, "minSubsetPartC", "minSubsetPart", "")

        ma = self.mainArea
        #boxL = QVBoxLayout(ma, QVBoxLayout.TopToBottom)
        #box.setTitle("Results")

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
        
        #for header in ["Geneset", "NES", "ES", "P-value", "FDR", "Size", "Matched Size", "Genes"]:
            #self.listView.addColumn(header)
        self.listView.setSelectionMode(QAbstractItemView.NoSelection)
        #self.connect(self.listView, SIGNAL("selectionChanged ( QListViewItem * )"), self.ja)
        self.connect(self.listView, SIGNAL("itemSelectionChanged()"), self.newPathwaySelected)

        OWGUI.separator(ca)

        box = OWGUI.widgetBox(ca, 'Phenotypes')

        self.psel = PhenotypesSelection(box)
        
        self.resize(600,50)
 
        OWGUI.separator(ca)
        self.btnApply = OWGUI.button(ca, self, "&Compute", callback = self.compute, disabled=0)
        
        fileBox = OWGUI.widgetBox(ca, orientation='horizontal')
        OWGUI.button(fileBox, self, "Load", callback = self.loadData, disabled=0)
        OWGUI.button(fileBox, self, "Save", callback = self.saveData, disabled=0)
        
        gen1 = getGenesets()

        for name,genes in gen1.items():
            self.addGeneset(name, genes)

        self.addComment("Computation was not started.")
        
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
        if sys.platform == "darwin":
            startfile = user.home
        else:
            startfile = "."
                
        filename = str(QFileDialog.getOpenFileName(self, 'Open GSEA data', startfile, "GSEA files (*.gsea)"))
        if filename == "": return
        
        fp = open(filename, "rb")
        res = pickle.load(fp)
        
        try:
            dm = pickle.load(fp)
        except:
            dm = None
        
        fp.close()
        
        self.compute(res, dm)

    def newPathwaySelected(self):
        print "newPathwaySelected"
        qApp.processEvents()

        if not self.selectable:
            return

        outat = set([])
        for item in self.listView.selectedItems():
            iname = self.lwiToGeneset[item]
            outat.update(self.res[iname][6])
            
            
        dataOut =  dataWithAttrs(self.data,list(outat))
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
        vars.append(orange.EnumVariable("Collection", values = ["KEGG", "C5", "C4", "C3", "C2", "C1"]))
        vars.append(orange.EnumVariable("Subcollection", values = ["KEGG", "C5 bp", "C5 cc", "C5 mf", "C4 cm", "C4 cgn", "C3 mir", "C3 tft", "C2 cgp", "C2 cp", "C1"]))
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
            splitndx = name.find("]")
            subcollection = name[1:splitndx]
            collection = subcollection.split(" ")[0]
            name = name[splitndx + 2:]
            examples.append([name, collection, subcollection, nes, es, pval, min(fdr,1.0), str(os), str(ts),  ", ".join(genes)])

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
            splitndx = name.find("]")
            collection = name[1:splitndx]
            name = name[splitndx + 2:]
            item = QTreeWidgetItem(self.listView)
            item.setText(0, collection)
            item.setText(1, name)
            item.setText(2, "%0.3f" % nes)
            item.setText(3, "%0.3f" % es)
            item.setText(4, "%0.3f" % pval)
            item.setText(5, "%0.3f" % min(fdr,1.0))
            item.setText(6, str(os))
            item.setText(7, str(ts))
            item.setText(8, writeGenes(genes))

            self.lwiToGeneset[item] = name

    def addComment(self, comm):
        item = QTreeWidgetItem(self.listView)
        item.setText(0, comm)   

    def setSelMode(self, bool):
        if bool:
            self.selectable = True
            self.listView.setSelectionMode(QAbstractItemView.MultiSelection)
        else:
            self.selectable = False
            self.listView.setSelectionMode(QListView.NoSelection)

    def compute(self, res=None, dm=None):
        clearListView(self.listView)
        self.addComment("Computing...")

        self.resultsOut(None)

        qApp.processEvents()
        self.res = res
        self.dm = dm
        
        if self.res == None and self.data:
            self.setSelMode(False)

            pb = OWGUI.ProgressBar(self, iterations=self.perms+2)

            if hasattr(self, "btnApply"):
                self.btnApply.setFocus()

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
                dkwargs["classValues"] = self.psel.getSelection()            
 
            gso = obiGsea.GSEA(organism="hsa")
            gso.setData(self.data, **dkwargs)

            for name,genes in self.geneSets.items():
                gso.addGeneset(name, genes)
                qApp.processEvents()

            self.res = gso.compute(n=self.perms, callback=pb.advance, **kwargs)
            
            pb.finish()
            
        if self.res != None:
            if len(self.res) > 0:
                self.fillResults(self.res)
                self.setSelMode(True)
                resl = self.res.items()

                etres = self.exportET(resl)

                self.resultsOut(etres)
                if self.dm == None:
                    self.dm = self.exportDistanceMatrix(resl)
                    
                    for ex in etres:
                        ex.name = str(ex[0])
                    
                    self.dm.setattr("items", etres)

                self.genesetDistOut(self.dm)

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
                print "set classes"
                self.psel.setClasses(getClasses(data))

    def addGeneset(self, name, genes):
        self.geneSets[name] = genes


def unpckGS(filename):
    import pickle
    f = open(filename,'rb')
    return pickle.load(f)

def getGenesets():
    import orngOrangeFoldersQt4
    return unpckGS(orngOrangeFoldersQt4.__getDirectoryNames()["bufferDir"] + "/gsea/geneSets_Gsea_KEGGhsa.pck")

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

    d = orange.ExampleTable("demo.tab")
    ow.setData(d)

    #d = orange.ExampleTable("tmp.tab")
    #ow.setData(d)



    a.exec_loop()
    ow.saveSettings()
