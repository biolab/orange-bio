"""
<name>OWGsea</name>
<description>Gene Set Enrichment Analysis</description>
<contact>Marko Toplak (marko.toplak(@at@)gmail.com)</contact>
<priority>10000</priority>
"""

from OWWidget import *
import OWGUI
import orngGsea
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


class OWGsea(OWWidget):

    settingsList = [ "name", "perms", "minSubsetSize", "minSubsetSizeC", "maxSubsetSize", "maxSubsetSizeC", \
        "minSubsetPart", "minSubsetPartC", "ptype" ]

    def __init__(self, parent=None, signalManager = None, name='GSEA'):
        OWWidget.__init__(self, parent, signalManager, name)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples only selected genes", ExampleTable) ]

        self.res = None

        self.name = 'GSEA'
        self.minSubsetSize = 3
        self.minSubsetSizeC = True
        self.maxSubsetSize = 1000
        self.maxSubsetSizeC = True
        self.minSubsetPart = 10
        self.minSubsetPartC = True
        self.perms = 100

        self.permutationTypes =  [ ("Phenotype", "p"), ("Gene","g") ]
        self.ptype = 0

        self.correlationTypes = [ ("Signal2Noise", "s2n") ]
        self.ctype = 0

        #self.loadSettings()
        self.data = None
        self.geneSets = {}

        ca = self.controlArea
        ca.setMaximumWidth(250)

        box = QVGroupBox(ca)
        box.setTitle('Permutate')

        self.permTypeF = OWGUI.comboBox(box, self, "ptype", items=nth(self.permutationTypes, 0), \
            tooltip="Permutation type.")

        _ = OWGUI.spin(box, self, "perms", 50, 1000, orientation="horizontal", label="Times")

        OWGUI.separator(ca)

        box = QVGroupBox(ca)
        box.setTitle('Correlation Calculation')

        self.corTypeF = OWGUI.comboBox(box, self, "ctype", items=nth(self.correlationTypes, 0), \
            tooltip="Correlation type.")

        OWGUI.separator(ca)

        box = QVGroupBox(ca)
        box.setTitle('Subset Filtering')

        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Size", 1, 10000, "minSubsetSizeC", "minSubsetSize", "") #TODO check sizes
        _,_ = OWGUI.checkWithSpin(box, self, "Max. Subset Size", 1, 10000, "maxSubsetSizeC", "maxSubsetSize", "")
        _,_ = OWGUI.checkWithSpin(box, self, "Min. Subset Part (%)", 1, 100, "minSubsetPartC", "minSubsetPart", "")

        OWGUI.separator(ca)
        
        self.btnApply = OWGUI.button(ca, self, "&Apply Changes", callback = self.compute, disabled=0)

        ma = self.mainArea
        boxL = QVBoxLayout(ma, QVBoxLayout.TopToBottom)
        #box.setTitle("Results")

        self.listView = QListView(ma)
        for header in ["Geneset", "NES", "ES", "P-value", "Size", "Matched Size", "Genes"]:
            self.listView.addColumn(header)

        self.listView.setSelectionMode(QListView.NoSelection)

        self.connect(self.listView, SIGNAL("selectionChanged ( QListViewItem * )"), self.newPathwaySelected)
        boxL.addWidget(self.listView)
        self.resize(600,50)

        gen1 = getGenesets()

        for name,genes in gen1.items():
            self.addGeneset(name, genes)

        self.addComment("Computation was not started.")

    def newPathwaySelected(self, item):

        qApp.processEvents()

        if not self.selectable:
            return

        iname = self.lwiToGeneset[item]
        outat = self.res[iname][6]

        dataOut =  dataWithAttrs(self.data, outat)
        self.send("Examples only selected genes", dataOut)

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
            #item.setText(4, "%0.3f" % min(fdr,1.0))
            item.setText(4, str(os))
            item.setText(5, str(ts))
            item.setText(6, writeGenes(genes))

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

        clearListView(self.listView)
        self.addComment("Computing...")

        qApp.processEvents()

        if self.data:

            self.setSelMode(False)

            pb = OWGUI.ProgressBar(self, iterations=self.perms+2)

            if hasattr(self, "btnApply"):
                self.btnApply.setFocus()

            gso = orngGsea.GSEA(organism="hsa")
            gso.setData(self.data)

            for name,genes in self.geneSets.items():
                gso.addGeneset(name, genes)
                qApp.processEvents()

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
                kwargs["permutation"] =  \
                    ifr(self.permutationTypes[self.ptype][1] == "p", "class", "a")
                
            res = gso.compute(n=self.perms, callback=pb.advance, **kwargs)
            self.res = res

            pb.finish()

            if len(res) > 0:
                self.fillResults(res)
                self.setSelMode(True)
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
            else:
                #enable correlation type
                comboboxItems(self.corTypeF, nth(self.correlationTypes, 0))
                self.corTypeF.setDisabled(False)
                #allow change of permutation type
                self.permTypeF.setDisabled(False)

    def addGeneset(self, name, genes):
        self.geneSets[name] = genes


def unpckGS(filename):
    import pickle
    f = open(filename,'rb')
    return pickle.load(f)

genesetFile = "geneSets3.pck"

def getGenesets():
    import orngRegistry
    return unpckGS(orngRegistry.outputDir + "/geneSets3.pck")

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWGsea()
    a.setMainWidget(ow)

    d = orange.ExampleTable('DLBCL_200a.tab')
    #d = orange.ExampleTable('brown-selected.tab')
    #d = orange.ExampleTable('testCorrelated.tab')
    ow.setData(d)

    gen1 = getGenesets()

    for name,genes in gen1.items():
        ow.addGeneset(name, genes)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
