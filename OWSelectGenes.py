"""
<name>Select Genes</name>
<description>Select genes from chip based on input selectors.</description>
<category>Genomics</category>
<icon>icons/SelectGenes.png</icon>
<priority>1150</priority>
"""

from OWWidget import *
import OWGUI
from OWChipDataFiles import ChipData
from OWChipANOVA import GeneSelection
import chipstat

class OWSelectGenes(OWWidget):
    settingsList  = ['negate', 'commitOnChange']

    def __init__(self, parent=None, name='Select Genes'):
        OWWidget.__init__(self, parent, name, "Select genes from chip based on input selectors")
        self.callbackDeposit = []

        self.inputs = [("Gene Selection", GeneSelection, self.loadselection, 0), ("Structured Chip Data", ChipData, self.chipdata, 1)]
        self.outputs = [("Selected Data", ChipData), ("Other Data", ChipData), ("Gene Selection", GeneSelection)]

        self.negate = 1
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

        # gene selection
        box = QVGroupBox("Gene Selection", self.controlArea)
        OWGUI.checkBox(box, self, 'negate', 'Negate', callback = self.selectionChange)

        # output
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

        self.resize(200,100)
        
    def chipdata(self, data):
        if data:
            self.data = data
##            nfiles = 0
##            for (n, d) in data:
##                nfiles += len(d)
##            self.infoa.setText("Microarray data, %d strains, total of %d data files" % (len(data), nfiles))
##            d = data[0][1][0]
##            self.infob.setText("Each data file contains %d measurements of %d genes" % (len(d.domain.attributes), len(d)))
            self.senddata()
        else:
            self.send("Selected Data", None)
            self.send("Other Data", None)
            self.send("Gene Selection", None)
        self.commitBtn.setEnabled(data <> None)

    def loadselection(self, selector, id):
        if selector:
            self.selectors[id] = selector
        else:
            del self.selectors[id]
        print 'SELN',
        for s in self.selectors.values():
            print len(s),
        print
##        print self.selectors
        self.infoa.setText('%d selectors on input' % len(self.selectors))
        self.senddata()
        self.commitBtn.setEnabled(self.selectors <> None and len(self.selectors))

    def senddata(self):
        if len(self.selectors):
            self.progressBarInit()
            s = self.selectors.values()
##            match = map(lambda *args: min(args), self.selectors.values())
            # this needs to be changed!!!
            match = []
            n = len(s)
            for i in range(len(s[0])):
                mm = []
                for j in range(n):
                    mm.append(s[j][i])
                match.append(min(mm))
                
            if self.negate:
                match = [(1,0)[x] for x in match]
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
            self.progressBarFinished()

    def selectionChange(self):
        pass

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWSelectGenes()
    a.setMainWidget(ow)
    ow.show()
    a.exec_loop()
    ow.saveSettings()
