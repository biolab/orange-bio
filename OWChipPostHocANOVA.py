"""
<name>Post Hoc ANOVA on Chip Data</name>
<description>Post Hoc ANOVA on chip data (strains, replicas).</description>
<category>Genomics</category>
<icon>icons/ChipPostHocANOVA.png</icon>
<priority>1130</priority>
"""

from OWWidget import *
import OWGUI
from OWChipDataFiles import ChipData
from OWChipANOVA import ANOVAResults, GeneSelection
import chipstat

class OWChipPostHocANOVA(OWWidget):
    settingsList  = ['commitOnChange', 'f1', 'f2', 'p1', 'p2', 'p3', 'filter1', 'filter2', 'filter3', 'selectorName', 'updateSelectorName']

    def __init__(self, parent=None, name='ANOVA on Chip Data'):
        OWWidget.__init__(self, parent, name, "ANOVA on chip data (strains, replicas).")
        self.callbackDeposit = []

        self.inputs = [("Results of ANOVA", ANOVAResults, self.anovaresults, 1)]
        self.outputs = [("PostHoc Gene Selection", GeneSelection)]

        self.commitOnChange = 0
        self.anova = None

        self.p1, self.p2, self.p3 = (2, 2, 2)
        self.filter1, self.filter2, self.filter3 = (1, 0, 0)
        self.f1 = None; self.f2 = None
        self.selectorName = "PostHoc Test"; self.updateSelectorName = 1
        # Settings
        self.loadSettings()
        self.ps = None
        self.selection = [None]*3

        # GUI
        # info
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('Total of 0 genes selected', box)
        OWGUI.separator(self.controlArea)

        # factors
        box = QVGroupBox("Factors", self.controlArea)
        self.f1Combo = OWGUI.comboBox(box, self, 'f1', sendSelectedValue=1, callback=[self.setFactorCombos, self.analysis], label="First:", labelWidth=50, orientation="horizontal")
        self.f1Combo.setDisabled(1)
        self.f2Combo = OWGUI.comboBox(box, self, 'f2', sendSelectedValue=1, callback=[self.setFactorCombos, self.analysis], label="Second:", labelWidth=50, orientation="horizontal")
        self.f2Combo.setDisabled(1)

        # gene selection
        self.selectionBox = QVGroupBox("Gene Selection", self.controlArea)
        self.factors = [('First factor (time)', 'p1', 'filter1'),
                   ('Second factor (strain)', 'p2', 'filter2'),
                   ('Interaction factor (time*strain)', 'p3', 'filter3')]
        self.pvals = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]

        self.numgenes = []
        for i in range(3):
            OWGUI.checkBox(self.selectionBox, self, self.factors[i][2], self.factors[i][0], callback=self.finalselection)
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

        self.resize(200,100)

    def setFactorCombos(self):

        def setOneCombo(combo, exclude):
            combo.clear()
            labels = [x for x in self.anova.classNames]
            if exclude in self.anova.classNames:
                labels.remove(exclude)
            for l in labels:
                combo.insertItem(l)

        items = self.anova.classNames
        if not self.f1 in items:
            self.f1 = items[items[0]==self.f2]
        if not self.f2 in items:
            self.f2 = items[items[0]==self.f2]

        setOneCombo(self.f1Combo, self.f2)
        setOneCombo(self.f2Combo, self.f1)
        self.f1 = self.f1 # this looks stupid, but it's not:
        self.f2 = self.f2 # it sets the right item to the combo boxes

    def anovaresults(self, anova):
        self.anova = anova
        self.selectionBox.setEnabled(anova <> None)
        self.commitBtn.setEnabled(anova <> None)
        if anova:
            self.setFactorCombos()
            self.f1Combo.setEnabled(1)
            self.f2Combo.setEnabled(1)
            self.analysis(commit=1)
            self.infoa.setText('ANOVA results on %d strains' % len(self.anova.classNames))
        else:
            self.send("PostHoc Gene Selection", None)
            self.f1Combo.clear()
            self.f1Combo.setDisabled(1)
        

    def geneselection(self, indx, commit=0):            
        margin = self.pvals[getattr(self, self.factors[indx][1])]
        p = [x[indx] for x in self.ps]
        n = len(filter(lambda x: x<margin, p))
        self.numgenes[indx].setText('  (%d %s)' % (n, ['genes', 'gene'][n==1]))
        self.selection[indx] = map(lambda x: x<margin, p)
        for x in self.selection:
            if x==None: return
        self.finalselection(commit)

    def finalselection(self, commit=0):
        if not self.ps:
            return
        self.match = [0]*len(self.ps)
        for indx in range(3):
            if getattr(self, self.factors[indx][2]) and self.selection[indx]:
                self.match = map(lambda a,b: a or b, self.match, self.selection[indx])
        n = self.match.count(1)
        self.infob.setText('Total of %d %s match criteria' % (n, ['genes', 'gene'][n==1]))
        if self.commitOnChange or commit:
            self.senddata()

    def analysis(self, commit=0):
        self.setSelectorName()
        self.progressBarInit()
        if not self.anova:
            return
        pbStep = 90./len(self.anova.anovaList)
        self.ps = chipstat.posthoc_anova_on_genes(self.f1, self.f2, self.anova, callback=lambda: self.progressBarAdvance(pbStep))
        self.progressBarSet(90)
        self.selection = [None]*3
        for indx in range(3):
            self.geneselection(indx, commit=1)
        self.progressBarFinished()

    def senddata(self):
        self.send("PostHoc Gene Selection", (self.selectorName, self.match))

    def setSelectorName(self):
        if self.updateSelectorName:
            s = 'PostHoc %s, %s' % (self.f1, self.f2)
##            ss = []
##            for indx in range(3):
##                if getattr(self, self.factors[indx][2]):
##                    ss.append(self.factors[indx][4] % self.pvals[getattr(self, self.factors[indx][1])])
##            if len(ss): s += reduce(lambda x,y:x+', '+y, ss)
            self.selectorName = s

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWChipPostHocANOVA()
    a.setMainWidget(ow)
    ow.show()
    a.exec_loop()
    ow.saveSettings()
