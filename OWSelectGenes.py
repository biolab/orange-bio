"""
<name>Select Genes</name>
<description>Select genes from chip based on input selectors.</description>
<icon>icons/SelectGenes.png</icon>
<priority>1150</priority>
"""

from OWWidget import *
from qttable import *
import OWGUI
from OWChipDataFiles import ChipData
from OWChipANOVA import GeneSelection
import chipstat

class OWSelectGenes(OWWidget):
    settingsList  = ['negate', 'commitOnChange']

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Select Genes')
        self.callbackDeposit = []

        self.inputs = [("Gene Selection", GeneSelection, self.loadselection, 0), ("Structured Chip Data", ChipData, self.chipdata)]
        self.outputs = [("Gene Selection", GeneSelection), ("Selected Data", ChipData), ("Other Data", ChipData)]

        self.negate = 0
        self.selectors = {}
        self.data = None
        self.commitOnChange = 1
        # Settings
        self.loadSettings()

        # GUI
        # info
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('', box)
        OWGUI.separator(self.controlArea)
        box.setMinimumWidth(170)

        # gene selection
        self.layout=QVBoxLayout(self.mainArea)
        box = QVGroupBox("Gene Selection", self.mainArea)
        self.table=QTable(box)
        self.table.setSelectionMode(QTable.NoSelection)
        self.layout.add(box)
        self.table.hide()
        self.drawtable()

        OWGUI.checkBox(box, self, 'negate', 'Negate', callback = self.selectionChange)

        # output
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

        self.resize(700,100)
        
    def chipdata(self, data):
        if data:
            self.data = data
            if self.commitOnChange and len(self.selectors)>0 and len(self.data[0][1][0]) == len(self.selectors.values()[0][1]):
                self.senddata()
        else:
            self.send("Selected Data", None)
            self.send("Other Data", None)
            self.send("Gene Selection", None)
        self.commitBtn.setEnabled(data <> None)

    def loadselection(self, selector, id):
        if selector:
            self.selectors[id] = list(selector + (1,0))    # the last two items for use and negation
        else:
            if id in self.selectors.keys(): del self.selectors[id]
        self.infoa.setText('%d selectors on input.' % len(self.selectors))
        self.drawtable()
##        print "debug OWSelectGenes.loadselection: self.commitOnChange: %s, len(self.selectors)>0: %s, self.data: %s, self.selectors.values(): %s" % (str(self.commitOnChange), str(len(self.selectors)>0), str(self.data), str(self.selectors.values()))
        if self.commitOnChange and len(self.selectors)>0 and self.data and len(self.data[0][1][0]) == len(self.selectors.values()[0][1]):
            self.senddata()
        self.commitBtn.setEnabled(self.selectors <> None and len(self.selectors))

    def drawtable(self):
        header = ['', 'Selector Name', 'Match #', 'Neg']
        self.table.setNumCols(len(header))
        self.header=self.table.horizontalHeader()
        for i in range(len(header)):
            self.header.setLabel(i, header[i])

        self.table.setNumRows(len(self.selectors))
        self.usesel_callback = [None] * len(self.selectors)
        self.tmpWidgets = []
        for (i, sel) in enumerate(self.selectors.values()):
            for (pos, k) in [(0,2), (3,3)]:
                cb = QCheckBox('', None)
                cb.setChecked(sel[k])
                self.usesel_callback[i] = lambda x, i=sel, k=k: self.usesel(x, i, k)
                self.connect(cb, SIGNAL("toggled(bool)"), self.usesel_callback[i])
                self.table.setCellWidget(i, pos, cb)
                self.tmpWidgets.append(cb)
            
            self.table.setText(i, 1, sel[0])
            self.table.setText(i, 2, '%d (%4.1f%s)' % (sel[1].count(1), 100.*sel[1].count(1)/len(sel[1]), '%') )
        self.table.setLeftMargin(0)
        self.table.verticalHeader().hide()

        for i in range(len(header)):
            self.table.adjustColumn(i)

        self.table.show()

    def usesel(self, x, sel, k):
        sel[k] = int(x)
        if self.commitOnChange:
            self.senddata()

    def senddata(self):
        if len(self.selectors):
            self.progressBarInit()
##            s = self.selectors.values()

            sel = []
            for s in self.selectors.values():
                if s[2]: # used?
                    if s[3]: # negated?
                        sel.append([not x for x in s[1]])
                    else:
                        sel.append(s[1])
            if len(sel) == 0:
                match = [0]*len(self.selectors.values()[0][1])
            elif len(sel)>1:
##                match = apply(map, (max, ) + tuple(sel))
                match = apply(map, (min, ) + tuple(sel))
            else:
                match = sel[0]
            if self.negate:
                match = [not x for x in match]

            nmatch = match.count(1)
            self.infob.setText("%d genes (%4.1f%s) match criteria" % (nmatch, 100.*nmatch/len(self.selectors.values()[0][1]), '%'))

            self.match = match
            self.progressBarSet(20)

            self.send("Gene Selection", self.match)
            if self.data:
                n = self.match.count(1)
                if n==len(self.data[0][1][0]): # all genes match
                    self.progressBarFinished()
                    self.send("Selected Data", self.data)
                    self.send("Other Data", None)
                elif n==0:
                    self.progressBarFinished()
                    self.send("Selected Data", None)
                    self.send("Other Data", self.data)
                else:
                    try:
                        P = []; N = []
                        pbStep = 80. / len(self.data)
                        for (strainname, data) in self.data:
                            dataP = [d.select(self.match) for d in data]
                            dataN = [d.select(self.match, negate=1) for d in data]
                            for i in range(len(data)):
                                dataP[i].name = dataN[i].name = data[i].name
                            P.append((strainname, dataP))
                            N.append((strainname, dataN))
                            self.progressBarAdvance(pbStep)
                        self.send("Selected Data", P)
                        self.send("Other Data", N)
                    except:
                        print "debug OWSelectGenes.senddata\n\tlen(d): %i\n\tlen(self.match): %i\n\t" % (len(d), len(self.match))
                        raise
            self.progressBarFinished()

    def selectionChange(self):
        if self.commitOnChange:
            self.senddata()

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWSelectGenes()
    a.setMainWidget(ow)
    ow.show()
    a.exec_loop()
    ow.saveSettings()
