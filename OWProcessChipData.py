"""
<name>Process Chip Data</name>
<description>Pre- and post-processing of chip data.</description>
<category>Genomics</category>
<icon>icons/ProcessChipData.png</icon>
<priority>1110</priority>
"""

from OWWidget import *
import OWGUI
from OWChipDataFiles import ChipData
import chipstat

class OWProcessChipData(OWWidget):
    settingsList  = ["preStdMethod", "postStdMethod", "preStdRob", "postStdRob", "mergeType", "commitOnChange"]

    def __init__(self, parent=None, name='Chip Data Files'):
        OWWidget.__init__(self, parent, name, "Process Chip Data pre- and post-processes chip data.")
        self.callbackDeposit = []

        self.inputs = [("Structured Chip Data", ChipData, self.chipdata)]
        self.outputs = [("Structured Chip Data", ChipData)]

        self.chipdata = None; self.datasets = None
        self.std = [("No preprocessing", None),
                    ("Array-based standardization", chipstat.standardize_arrays),
                    ("Gene-based standardization", chipstat.standardize_genes),
                    ("First array-, then gene-based standardization", lambda e,r: chipstat.standardize_genes(chipstat.standardize_arrays(e,r),r)),
                    ("First gene-, then array-based standardization", lambda e,r: chipstat.standardize_arrays(chipstat.standardize_genes(e,r),r))]
        # Settings
        self.data = None
        self.preStdMethod = 0; self.preStdRob = 1
        self.postStdMethod = 0; self.postStdRob = 1
        
        self.mergeType = 0
        self.commitOnChange = 0
        self.loadSettings()

        # GUI
        # info
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('', box)
        
        # preprocessing
        OWGUI.separator(self.controlArea)
        box = QVGroupBox("Preprocessing", self.controlArea)
        labels = [x[0] for x in self.std]
        OWGUI.comboBox(box, self, 'preStdMethod', label=None, labelWidth=None, orientation='vertical', items=labels, callback=self.selectionChange)
        self.preRobBtn = OWGUI.checkBox(box, self, "preStdRob", "Robust standardization", callback=self.selectionChange)
        
        # merge
        OWGUI.separator(self.controlArea)
        self.mergeTypes = [('mean', 'Mean'), ('median', 'Median'), ('min', 'Minimum expression'), ('max', 'Maximum expression')]
        labels = [x[1] for x in self.mergeTypes]
        OWGUI.radioButtonsInBox(self.controlArea, self, 'mergeType', labels, box='Merge Replicas', tooltips=None, callback=self.selectionChange)

        # postprocessing
        OWGUI.separator(self.controlArea)
        box = QVGroupBox("Postprocessing", self.controlArea)
        labels = [x[0] for x in self.std]
        OWGUI.comboBox(box, self, 'postStdMethod', label=None, labelWidth=None, orientation='vertical', items=labels, callback=self.selectionChange)
        self.postRobBtn = OWGUI.checkBox(box, self, "postStdRob", "Robust standardization", callback=self.selectionChange)

        # output
        OWGUI.separator(self.controlArea)
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.selectionChange, disabled=1)

        self.setBtnsState()
        self.resize(100,100)

    def selectionChange(self):
        self.setBtnsState()
        if self.commitOnChange:
            self.sendData()

    def chipdata(self, data):
        self.commitBtn.setEnabled(data <> None)
        if data:
            self.data = data
            nfiles = 0
            for (n, d) in data:
                nfiles += len(d)
            self.infoa.setText("Microarray data, %d strains, total of %d data files" % (len(data), nfiles))
            print data
            d = data[0][1][0]
            self.infob.setText("Each data file contains %d measurements of %d genes" % (len(d.domain.attributes), len(d)))

            self.sendData()
        else:
            self.send("Structured Chip Data", None)

    # process arrays in the structure, returns new structure
    def processArrays(self, datastructure, method, *arg):
        pbStep = 30. / sum([len(data) for (dummy, data) in datastructure])
        newdata = []
        for (strainname, data) in datastructure:
            if data:
                new = []
                for e in data:
                    new.append(apply(method, (e,)+arg))
                    self.progressBarAdvance(pbStep)
                newdata.append([strainname, new])
        return newdata

    def sendData(self):
        if not self.data:
            return
        self.send('Structured Chip Data', None) # this is required for other widgets not to mess up with two different datasets
        self.progressBarInit()
        data = self.data
        # preprocessing
        if self.preStdMethod: # 0 means no preprocessing
            data = self.processArrays(data, self.std[self.preStdMethod][1], self.preStdRob)
        self.progressBarSet(30)
        # merging
        merged = chipstat.merge_replicas(data, self.mergeTypes[self.mergeType][0])
        self.progressBarSet(70)
        # postprocessing
        if self.postStdMethod: # 0 means no preprocessing
            merged = self.processArrays(merged, self.std[self.postStdMethod][1], self.preStdRob)
        for (i,d) in enumerate(merged):
            d[1][0].name = self.data[i][0]
        self.progressBarFinished()
        self.send('Structured Chip Data', merged)
            
    def setBtnsState(self):
        self.preRobBtn.setEnabled(self.preStdMethod > 0)
        self.postRobBtn.setEnabled(self.postStdMethod > 0)

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWProcessChipData()
    a.setMainWidget(ow)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
