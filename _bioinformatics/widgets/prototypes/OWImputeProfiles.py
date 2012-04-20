## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

"""
<name>Impute Profiles</name>
<description>Imputation and non-parametric smoothing of expression profiles.</description>
<icon>icons/ImputeProfiles.png</icon>
<priority>1065</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

import math

import numpy.oldnumeric as Numeric, numpy.oldnumeric.ma as MA

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import chipimpute, chipstat
from .OWDataFiles import DataFiles

class OWImputeProfiles(OWWidget):
    settingsList  = ['impute', 'imputeK', 'smooth', 'windowSize', 'commitOnChange']

    def __init__(self, parent=None, signalManager = None):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, signalManager, 'Impute & Loess Profiles')

        self._data = None       # exampleTable
        self._dataMA = None     # input 2d masked array
        self._chipdata = None   # [(dirname0, [et0, et1, ...]), ...]
        self._chipdataMA = []   # [(dirname0, [m2d0, m2d1, ...]), ...]
        self.impute = 1
        self.imputeK = 20
        self.smooth = 1
        self.windowSize = 3
        self.commitOnChange = 1
               
        # Settings
        self.loadSettings()

        # GUI
        self.mainArea.setFixedWidth(0)
        ca=QFrame(self.controlArea)
        gl=QGridLayout(ca,4,1,5)
        # info
        box = QVGroupBox("Info", ca)
        gl.addWidget(box, 0,0)
        self.infoa = QLabel("No examples on input", box)
        self.infob = QLabel("", box)
        QLabel("", box)
        self.infoc = QLabel("No structured data on input", box)
        self.infod = QLabel("", box)
        # KNN impute
        self.boxImpute = QVGroupBox("Impute missing values", ca)
        gl.addWidget(self.boxImpute, 1,0)
        OWGUI.checkBox(self.boxImpute, self, "impute", "KNN impute", tooltip="Impute missing values from values of K nearest neighbours.", callback=self.change)
##        self.sliderK = OWGUI.hSlider(self.boxImpute, self, "imputeK", box=None, minValue=1, maxValue=7744, step=1, callback=self.imputeChange, labelFormat=" K = %i", ticks=0)
##        self.sliderK = OWGUI.qwtHSlider(self.boxImpute, self, "imputeK", box=None, label="K", labelWidth=12, minValue=1, maxValue=7744, step=0.02, precision=0, callback=self.imputeChange, logarithmic=1, ticks=0, maxWidth=200)
##        self.sliderK = OWGUI.qwtHSlider(self.boxImpute, self, "imputeK", box=None, label="K", labelWidth=12, minValue=1, maxValue=100, step=1, precision=0, callback=self.imputeKChange, logarithmic=0, ticks=0, maxWidth=200)
        self.sliderK = OWGUI.qwtHSlider(self.boxImpute, self, "imputeK", box=None, label="K", labelWidth=15, minValue=1, maxValue=999, step=1, precision=0, callback=self.imputeKChange, logarithmic=0, ticks=0, maxWidth=None)
        self.boxImpute.setDisabled(1)
        # loess
        self.boxLoess = QVGroupBox("Smoothing", ca)
        gl.addWidget(self.boxLoess, 2,0)
        OWGUI.checkBox(self.boxLoess, self, "smooth", "Loess smoothing", tooltip="Loess profiles, impute missing columns.", callback=self.change)
        lbl = QLabel("Window size (number of points)", self.boxLoess)
        lbl.setAlignment(Qt.AlignHCenter)
        self.sliderW = OWGUI.qwtHSlider(self.boxLoess, self, "windowSize", box=None, label="W", labelWidth=15, minValue=1, maxValue=999, step=1, precision=0, callback=self.smoothWChange, logarithmic=0, ticks=0, maxWidth=None)
        self.boxLoess.setDisabled(1)
        # output
        box = QVGroupBox("Output", ca)
        gl.addWidget(box, 3,0)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

        self.inputs = [("Examples", ExampleTable, self.data), ("Structured Data", DataFiles, self.chipdata)]
        self.outputs = [("Examples", ExampleTable), ("Structured Data", DataFiles)]

        # data dependent variables
        self.numRowsMissing = 0
        self.numRowsMissingChipData = 0
        self.resize(100,100)

    def change(self):
        if self.commitOnChange:
            self.senddata(0)

    def imputeKChange(self):
        if self.commitOnChange and self.impute:
            self.senddata(0)

    def smoothWChange(self):
        if self.commitOnChange and self.smooth:
            self.senddata(0)

    def setGuiCommonExpChip(self):
        """Sets GUI features that depend on both inputs, e.g. max. values of sliders, commit button.
        """
        # calculate maxK, maxW
        maxK = 1e9
        maxW = 1e9
        if self._dataMA != None:
            maxK = min(max(self._dataMA.shape[0]/50, 50), self._dataMA.shape[0], maxK)
            maxW = min(self._dataMA.shape[1], maxW)
        if self._chipdataMA:
            maxK = min(max(self._chipdataMA[0][1][0].shape[0]/50., 50), self._chipdataMA[0][1][0].shape[0], maxK)
            maxW = min(self._chipdataMA[0][1][0].shape[1], maxW)
        if maxK == 1e9: maxK = 1
        if maxW == 1e9: maxW = 1
        # enable sliders and commit button
        if self._dataMA != None or self._chipdataMA:
            # set up sliders (only if data present)
            self.sliderK.setRange(1, maxK, 1)
            if self.imputeK > maxK:
                self.sliderK.setValue(maxK)
                self.imputeK = maxK
            self.sliderW.setRange(1, maxW, 1)
            if self.windowSize > maxW:
                self.sliderW.setValue(maxW)
                self.windowSize = maxW
            self.boxImpute.setEnabled(1)
            self.boxLoess.setEnabled(1)
            self.commitBtn.setEnabled(1)
        else:
            self.boxImpute.setDisabled(1)
            self.boxLoess.setDisabled(1)
            self.commitBtn.setDisabled(1)

    def data(self, data):
        if data != None:
            self._data = data
##            self._dataMA = chipstat.orng2ma(data)
            self._dataMA = data.toNumpyMA("a")[0]
            # info text
            self.infoa.setText("Examples: %i profiles on %i points" % (self._dataMA.shape[0], self._dataMA.shape[1]))
            numTotalMissing = int(Numeric.multiply.reduce(self._dataMA.shape) - MA.count(self._dataMA))
            if numTotalMissing > 0:
                numValsByCol = MA.count(self._dataMA, 0)
                numEmptyCol = Numeric.add.reduce(Numeric.where(numValsByCol==0, 1, 0))
                colNonEmpty = Numeric.compress(numValsByCol!=0, Numeric.arange(self._dataMA.shape[1]))
                dataRemEmptyCol = self._dataMA.take(colNonEmpty, 1)
                self.numRowsMissing = Numeric.add.reduce(Numeric.where(MA.count(dataRemEmptyCol, 1) < dataRemEmptyCol.shape[1], 1, 0))
                s1 = ""
                s2 = ""
                if numEmptyCol > 0: s1 = "s"
                if self.numRowsMissing > 0: s2 = "s"
                self.infob.setText("missing %i values, %i column%s completely, %i row%s partially" % (numTotalMissing, numEmptyCol, s1, self.numRowsMissing, s2))
            else:
                self.infob.setText("")
        else:
            self._data = None
            self._dataMA = None
            self.infoa.setText("No examples on input")
            self.infob.setText("")
            self.numRowsMissing = 0
        self.setGuiCommonExpChip()
        if self.commitOnChange:
            self.senddata(1)


    def chipdata(self, data):
        """Input data: [(dirname0, [et0, et1, ...]), ...]
        """
        self.numRowsMissingChipData = 0
        self._chipdataMA = []
        if data != None:
            self._chipdata = data
            numValsAll = 0
            numValsNonMasked = 0
            numFiles = 0
            numExamplesList = []
            attribDict = {}
            numColMissing = 0
            for (name, etList) in data:
                numFiles += len(etList)
                self._chipdataMA.append((name,[]))
                for et in etList:
                    attribDict.update(dict(zip(map(lambda x: x.name, et.domain.attributes), et.domain.attributes)))
                    numExamplesList.append(len(et))
                    etm = et.toNumpyMA("a")[0]
                    colNonMissingInd = Numeric.compress(Numeric.not_equal(MA.count(etm, 0), 0), Numeric.arange(etm.shape[1])) # indices of columns that are not completely missing
                    numColMissing += etm.shape[1] - colNonMissingInd.shape[0]
                    self.numRowsMissingChipData += int(Numeric.add.reduce(Numeric.less(MA.count(etm.take(colNonMissingInd, 1), 1), etm.shape[1])))
                    numValsAll += int(Numeric.multiply.reduce(etm.shape))
                    numValsNonMasked += int(MA.count(etm))
                    self._chipdataMA[-1][1].append(etm)
            # info text
            self.infoc.setText("Structured Data: %i data files with %i profiles on %i points" % (numFiles, numExamplesList[0], len(attribDict)))
            numTotalMissing = numValsAll-numValsNonMasked
            if numTotalMissing > 0:
                print numTotalMissing, numColMissing, self.numRowsMissingChipData
                print type(numTotalMissing), type(numColMissing), type(self.numRowsMissingChipData)
                self.infod.setText("missing %i values, %i column%s completely, %i row%s partially" % (numTotalMissing, numColMissing, ["","s"][numColMissing!=1], self.numRowsMissingChipData, ["","s"][self.numRowsMissingChipData!=1]))
            else:
                self.infod.setText("")                
        else:
            self._chipdata = None
            self.infoc.setText("No structured data on input")
            self.infod.setText("")
        self.setGuiCommonExpChip()
        if self.commitOnChange:
            self.senddata(2)


    def senddata(self, outputSelector=0):
        """outputSelector = [0: examples + chip data | 1: examples | 2: chip data]
        """
        assert outputSelector in [0,1,2]
        # progress bar settings: 1: examples, 2: chip data
        steps = 0
        if outputSelector in [0,1]:
            if self.impute:
                steps += self.numRowsMissing
            if self.smooth and self._dataMA != None:
                steps += self._dataMA.shape[0]
        if outputSelector in [0,2]:
            if self.impute:
                steps += self.numRowsMissingChipData
            if self.smooth and self._chipdataMA:
                steps += len(self._chipdataMA[0][1][0]) * reduce(lambda x,y: x+len(y[1]), self._chipdataMA, 0) #self._chipdataMA.shape[2]
        if steps == 0: steps = 1
        pbStep = 100./steps
        self.progressBarInit()
        if outputSelector in [0,1]:
            self._sendexampledata(pbStep)
        if outputSelector in [0,2]:
            self._sendchipdata(pbStep)
        self.progressBarFinished()

    def _sendexampledata(self, pbStep):
        if self._dataMA != None:
            ma2d = self._dataMA
            if self.impute:
                ma2d = chipimpute.kNNimputeMA(ma2d, int(self.imputeK), callback=lambda: self.progressBarAdvance(pbStep))
            if self.smooth:
                ws = max(1.1, int(self.windowSize)) # fixes bug in statc.loess
                ma2d = chipimpute.loessMA(ma2d, ws, 1, callback=lambda: self.progressBarAdvance(pbStep))
            et = chipstat.ma2orng_keepClassMetas(ma2d, self._data)
            et.name = self._data.name
            self.send("Examples", et)
        else:
            self.send("Examples", None)


    def _sendchipdata(self, pbStep):
        if self._chipdataMA:
            chipdataNew = []    # [(dirname0, [et0, et1, ...]), ...]
            ws = max(1.1, int(self.windowSize)) # fixes bug in statc.loess
            for sIdx, (name, etmList) in enumerate(self._chipdataMA):
                chipdataNew.append((name,[]))
                for eIdx, etm in enumerate(etmList):
                    if self.impute:
                        etm = chipimpute.kNNimputeMA(etm, int(self.imputeK), callback=lambda: self.progressBarAdvance(pbStep))
                    if self.smooth:
                        etm = chipimpute.loessMA(etm, ws, 1, callback=lambda: self.progressBarAdvance(pbStep))
##                   chipdataNew[-1][1].append(orange.ExampleTable(et.domain, etm))
                    etNew = chipstat.ma2orng_keepClassMetas(etm, self._chipdata[sIdx][1][eIdx])
                    etNew.name = self._chipdata[sIdx][1][eIdx].name
                    chipdataNew[-1][1].append(etNew)
            self.send("Structured Data", chipdataNew)
        else:
            self.send("Structured Data", None)



if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWImputeProfiles()
    a.setMainWidget(ow)
    ow.show()
##    ow.data(orange.ExampleTable("meanExpr_ann_pkaC.tab"))
##    ow.data(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS_100_yakApufA\pufA\pufA1.1.xls.tab"))
##    ow.data(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\potato\I\I_30m_1_1042_teh1.tab"))
    ow.data(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Data Hs\2007-08-29 Banjo\- data\untr maxAbsDiff.tab"))
    a.exec_loop()
    ow.saveSettings()

##    from . import OWDataFiles
##    from Orange.orng import orngSignalManager
##    signalManager = orngSignalManager.SignalManager(0)
##    a=QApplication(sys.argv)
##    ow=OWImputeProfiles(signalManager = signalManager)
##    a.setMainWidget(ow)
##    ow.show()
##    ds = OWDataFiles.OWDataFiles(signalManager = signalManager)
####    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\2005-08 NIB Potato\potato.attrib")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\2005-08 NIB Potato\potato.sub100.avg")
##    signalManager.addWidget(ow)
##    signalManager.addWidget(ds)
##    signalManager.setFreeze(1)
##    signalManager.addLink(ds, ow, 'Structured Data', 'Structured Data', 1)
##    signalManager.setFreeze(0)
##    a.exec_loop()
##    ow.saveSettings()
