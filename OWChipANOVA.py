"""
<name>ANOVA on Chip Data</name>
<description>ANOVA on chip data (strains, replicas).</description>
<category>Genomics</category>
<icon>icons/ChipANOVA.png</icon>
<priority>1120</priority>
"""

from OWWidget import *
import OWGUI
from OWChipDataFiles import ChipData
import chipstat

class ANOVAResults(orange.Orange):
    pass

class GeneSelection(orange.Orange):
    pass

class OWChipANOVA(OWWidget):
    settingsList  = ['p1', 'p2', 'p3', 'filter1', 'filter2', 'filter3', 'commitOnChange', 'selectorName', 'updateSelectorName']

    def __init__(self, parent=None, name='ANOVA on Chip Data'):
        OWWidget.__init__(self, parent, name, "ANOVA on chip data (strains, replicas).")
        self.callbackDeposit = []

        self.inputs = [("Structured Chip Data", ChipData, self.chipdata, 1)]
        self.outputs = [("Results of ANOVA", ANOVAResults), ("Gene Selection", GeneSelection), ("Selected Data", ChipData), ("Other Data", ChipData)]

        self.commitOnChange = 0
        self.chipdata = None
        self.p1, self.p2, self.p3 = (3, 3, 3)
        self.filter1, self.filter2, self.filter3 = (0, 1, 1)
        self.selectorName = 'ANOVA'
        self.updateSelectorName = 1
        # Settings
        self.loadSettings()
        self.ps = None
        self.selection = [None]*3

        # GUI
        # info
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('', box)
        self.infoc = QLabel('Total of 0 genes match criteria', box)
        OWGUI.separator(self.controlArea)

        # gene selection
        self.selectionBox = QVGroupBox("Gene Selection", self.controlArea)
        self.factors = [('First factor (time)', 'p1', 'filter1', 'H0: The expression of a gene does not depend on\ntime through the course of development.', 'Time (P<%3.2f)'),
                   ('Second factor (condition)', 'p2', 'filter2', 'H0: The expression of a gene does not depend on\nthe texperimental condition (mutation).', 'Cond (P<%3.2f)'),
                   ('Interaction factor (time*condition)', 'p3', 'filter3', 'H0: The effect of experimental conditions on\nthe expression does not depend on time. ', 'Time x Cond (P<%3.2f)')]
        self.pvals = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]

        self.numgenes = []
        for i in range(3):
            OWGUI.checkBox(self.selectionBox, self, self.factors[i][2], self.factors[i][0], callback=self.finalselection, tooltip=self.factors[i][3])
            hbox = QHBox(self.selectionBox)
            lbl = QLabel('', hbox)
            lbl.setFixedSize(20, lbl.sizeHint().height())
            lbl = QLabel('p <', hbox)
            lbl.setFixedSize(25, lbl.sizeHint().height())
            OWGUI.comboBox(hbox, self, self.factors[i][1], items = self.pvals, callback=lambda x=i: self.geneselection(x))
            self.numgenes.append(QLabel('  (0 genes)', hbox))
        self.selectionBox.setDisabled(1)

        OWGUI.separator(self.selectionBox)
        self.selectorNameLE = OWGUI.lineEdit(self.selectionBox, self, 'selectorName', label='Selector Name: ')
        OWGUI.checkBox(self.selectionBox, self, 'updateSelectorName', 'Automatically update selector name')

        # output
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

##        OWGUI.button(self.controlArea, self, 'TEST', self.setSelectorName)        

        self.resize(200,100)
        
    def chipdata(self, data):
        if data:
            self.data = data
            self.selectionBox.setEnabled(1)
            nfiles = 0
            for (n, d) in data:
                nfiles += len(d)
            self.infoa.setText("Microarray data, %d strains, total of %d data files" % (len(data), nfiles))
            d = data[0][1][0]
            self.infob.setText("Each file contains %d measurements of %d genes" % (len(d.domain.attributes), len(d)))

            self.analysis()
        else:
            self.send("Results of ANOVA", None)
        self.commitBtn.setEnabled(data <> None)

    def geneselection(self, indx):            
        margin = self.pvals[getattr(self, self.factors[indx][1])]
        p = [x[indx] for x in self.ps]
        n = len(filter(lambda x: x<margin, p))
        self.numgenes[indx].setText('  (%d %s)' % (n, ['genes', 'gene'][n==1]))
        self.selection[indx] = map(lambda x: x<margin, p)
        for x in self.selection:
            if x==None: return
        self.finalselection()

    def finalselection(self):
        if not self.ps:
            return
        self.match = [0]*len(self.ps)
        for indx in range(3):
            if getattr(self, self.factors[indx][2]) and self.selection[indx]:
                self.match = map(lambda a,b: a or b, self.match, self.selection[indx])
        n = self.match.count(1)
        self.infoc.setText('Total of %d %s match criteria' % (n, ['genes', 'gene'][n==1]))
        self.setSelectorName()
        if self.commitOnChange:
            self.senddata()

    def analysis(self):
        self.progressBarInit()
        pbStep = 90./len(self.data[0][1][0])
        self.ps, anova_results = chipstat.anova_on_genes(self.data, callback=lambda: self.progressBarAdvance(pbStep))
        self.send("Results of ANOVA", anova_results)
        self.progressBarSet(90)
        self.selection = [None]*3
        for indx in range(3):
            self.geneselection(indx)
        self.progressBarFinished()

    def senddata(self):
        self.send("Gene Selection", (self.selectorName, self.match))
        n = self.match.count(1)
        if n==len(self.data[0][1][0]): # all genes match
            self.send("Selected Data", self.data)
            self.send("Other Data", None)
        elif n==0:
            self.send("Selected Data", None)
            self.send("Other Data", self.data)
        else:
            P = []; N = []
            for (strainname, data) in self.data:
                dataP = [d.select(self.match) for d in data]
                dataN = [d.select(self.match, negate=1) for d in data]
                for i in range(len(data)):
                    dataP[i].name = dataN[i].name = data[i].name
                P.append((strainname, dataP))
                N.append((strainname, dataN))
            self.send("Selected Data", P)
            self.send("Other Data", N)

    def setSelectorName(self):
        if self.updateSelectorName:
            s = 'ANOVA '
            if len(self.data): s += reduce(lambda x,y:x+', '+y, [x[0] for x in self.data][:2]) + ['. ', '... '][len(self.data)>2]
            ss = []
            for indx in range(3):
                if getattr(self, self.factors[indx][2]):
                    ss.append(self.factors[indx][4] % self.pvals[getattr(self, self.factors[indx][1])])
            if len(ss): s += reduce(lambda x,y:x+', '+y, ss)
            self.selectorName = s

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWChipANOVA()
    a.setMainWidget(ow)
    ow.show()
    a.exec_loop()
    ow.saveSettings()
