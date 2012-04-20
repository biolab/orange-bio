"""
<name>Process Chip Data</name>
<description>Pre- and post-processing of chip data.</description>
<icon>icons/ProcessChipData.png</icon>
<priority>1060</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

from __future__ import absolute_import

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import chipstat
from .OWDataFiles import DataFiles

class OWProcessChipData(OWWidget):
    settingsList  = ["preStdMethod", "postStdMethod", "preStdRob", "postStdRob", "mergeType", "commitOnChange"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Process Chip Data')
        self.callbackDeposit = []

        self.inputs = [("Structured Data", DataFiles, self.chipdata)]
        self.outputs = [("Structured Data", DataFiles)]

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
        self.mergeTypes = [(0, "No merging"), ('mean', 'Mean'), ('median', 'Median'), ('min', 'Minimum expression'), ('max', 'Maximum expression')]
        labels = [x[1] for x in self.mergeTypes]
        OWGUI.radioButtonsInBox(self.controlArea, self, 'mergeType', labels, box='Merge Replicas', tooltips=None, callback=self.selectionChange)

        # postprocessing
        OWGUI.separator(self.controlArea)
        self.boxPostproc = QVGroupBox("Postprocessing", self.controlArea)
        labels = [x[0] for x in self.std]
        OWGUI.comboBox(self.boxPostproc, self, 'postStdMethod', label=None, labelWidth=None, orientation='vertical', items=labels, callback=self.selectionChange)
        self.postRobBtn = OWGUI.checkBox(self.boxPostproc, self, "postStdRob", "Robust standardization", callback=self.selectionChange)

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
        print data
        self.commitBtn.setEnabled(data <> None)
        if data:
            self.data = data
            nfiles = 0
            for (n, d) in data:
                nfiles += len(d)
            self.infoa.setText("Structured data, %d sets, total of %d data files" % (len(data), nfiles))
            d = data[0][1][0]
            self.infob.setText("Each file contains %d attributes and %d examples" % (len(d.domain.attributes), len(d)))

            self.sendData()
        else:
            self.send("Structured Data", None)

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
        self.send('Structured Data', None) # this is required for other widgets not to mess up with two different datasets
        self.progressBarInit()
        data = self.data
        # preprocessing
        if self.preStdMethod: # 0 means no preprocessing
            data = self.processArrays(data, self.std[self.preStdMethod][1], self.preStdRob)
        self.progressBarSet(30)
        # merging
        if self.mergeType: data = chipstat.merge_replicas(data, self.mergeTypes[self.mergeType][0])
        self.progressBarSet(70)
        # postprocessing
        if self.postStdMethod: # 0 means no preprocessing
            data = self.processArrays(data, self.std[self.postStdMethod][1], self.preStdRob)
        for i,(strain, etList) in enumerate(data):
            if len(etList) == len(self.data[i][1]):
                for j,et in enumerate(etList):
                    et.name = self.data[i][1][j].name
            else:
                for et in etList:
                    et.name = strain
        self.progressBarFinished()
        self.send('Structured Data', data)
            
    def setBtnsState(self):
        self.preRobBtn.setEnabled(self.preStdMethod > 0)
        self.postRobBtn.setEnabled(self.postStdMethod > 0)
        self.boxPostproc.setEnabled(self.mergeType > 0)

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWProcessChipData()
    a.setMainWidget(ow)

    ow.show()
    a.exec_loop()
    ow.saveSettings()
