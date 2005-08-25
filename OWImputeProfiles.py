"""
<name>Impute Profiles</name>
<description>Imputation and non-parametric smoothing of expression profiles.</description>
<icon>icons/ImputeProfiles.png</icon>
<priority>1065</priority>
"""

import math
import Numeric, MA
from OWWidget import *
import OWGUI
from OWDataFiles import DataFiles
import chipstat
import chipimpute


class OWImputeProfiles(OWWidget):
    settingsList  = ['impute', 'imputeK', 'smooth', 'windowSize', 'commitOnChange']

    def __init__(self, parent=None, signalManager = None):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, signalManager, 'Impute & Loess Profiles')

        self._data = None       # exampleTable
        self._dataMA = None     # input 2d masked array
        self._chipdata = None       # [(dirname0, [et0, et1, ...]), ...]
        self._chipdataMA = None     # input 3d masked array
        self.impute = 1
        self.imputeK = 20
        self.smooth = 1
        self.windowSize = 3
        self.commitOnChange = 1
               
        # Settings
        self.loadSettings()

        # GUI
        # info
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel("No examples on input", box)
        self.infob = QLabel("", box)
        OWGUI.separator(box, 300)
        self.infoc = QLabel("No structured data on input", box)
        self.infod = QLabel("", box)
        OWGUI.separator(self.controlArea)
        # KNN impute
        self.boxImpute = QVGroupBox("Impute missing values", self.controlArea)
        OWGUI.checkBox(self.boxImpute, self, "impute", "KNN impute", tooltip="Impute missing values from values of K nearest neighbours.", callback=self.change)
##        self.sliderK = OWGUI.hSlider(self.boxImpute, self, "imputeK", box=None, minValue=1, maxValue=7744, step=1, callback=self.imputeChange, labelFormat=" K = %i", ticks=0)
##        self.sliderK = OWGUI.qwtHSlider(self.boxImpute, self, "imputeK", box=None, label="K", labelWidth=12, minValue=1, maxValue=7744, step=0.02, precision=0, callback=self.imputeChange, logarithmic=1, ticks=0, maxWidth=200)
        self.sliderK = OWGUI.qwtHSlider(self.boxImpute, self, "imputeK", box=None, label="K", labelWidth=12, minValue=1, maxValue=100, step=1, precision=0, callback=self.imputeKChange, logarithmic=0, ticks=0, maxWidth=200)
        self.boxImpute.setDisabled(1)
        OWGUI.separator(self.controlArea)
        # loess
        self.boxLoess = QVGroupBox("Smoothing", self.controlArea)
        OWGUI.checkBox(self.boxLoess, self, "smooth", "Loess smoothing", tooltip="Loess profiles, impute missing columns.", callback=self.change)
        lbl = QLabel("Window size (number of points)", self.boxLoess)
        lbl.setAlignment(Qt.AlignHCenter)
        self.sliderW = OWGUI.qwtHSlider(self.boxLoess, self, "windowSize", box=None, label="W", labelWidth=12, minValue=1, maxValue=10, step=1, precision=0, callback=self.smoothWChange, logarithmic=0, ticks=0, maxWidth=200)
        self.boxLoess.setDisabled(1)
        OWGUI.separator(self.controlArea)
        # output
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

        self.inputs = [("Examples", ExampleTable, self.data, 1), ("Structured Data", DataFiles, self.chipdata, 1)]
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
        if self._chipdataMA != None:
            maxK = min(max(self._chipdataMA.shape[0]/50, 50), self._chipdataMA.shape[0], maxK)
            maxW = min(self._chipdataMA.shape[1], maxW)
        if maxK == 1e9: maxK = 1
        if maxW == 1e9: maxW = 1
        # enable sliders and commit button
        if self._dataMA != None or self._chipdataMA != None:        
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
            self._dataMA = chipstat.orng2ma(data)
            # info text
            self.infoa.setText("Examples: %i profiles on %i points" % (self._dataMA.shape[0], self._dataMA.shape[1]))
            numTotalMissing = Numeric.multiply.reduce(self._dataMA.shape) - MA.count(self._dataMA)
            if numTotalMissing > 0:
                numValsByCol = MA.count(self._dataMA, 0)
                numEmptyCol = Numeric.add.reduce(Numeric.where(numValsByCol==0, 1, 0))
                colNonEmpty = Numeric.compress(numValsByCol!=0, Numeric.arange(self._dataMA.shape[1]))
                dataRemEmptyCol = MA.take(self._dataMA, colNonEmpty, 1)
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
        if data != None:
            self._chipdata = data
            shp = [0,0,0]
            shp[0] = len(data[0][1][0])
            shp[1] = len(data[0][1][0].domain.attributes)
            shp[2] = reduce(lambda x,y: x+len(y[1]), data, 0)
            self._chipdataMA = MA.zeros(shp, MA.Float)
            idx = 0
            for (name, etList) in data:
                for et in etList:
                    self._chipdataMA[:,:,idx] = chipstat.orng2ma(et)
                    idx += 1
            # info text
            self.infoc.setText("Structured Data: %i data files with %i profiles on %i points" % (shp[2], shp[0], shp[1]))
            numTotalMissing = Numeric.multiply.reduce(self._chipdataMA.shape) - MA.count(self._chipdataMA)
            if numTotalMissing > 0:
                numValsByCol = MA.count(self._chipdataMA, 0)
                emptyColInd = Numeric.compress(Numeric.ravel(numValsByCol)==0, Numeric.arange(Numeric.multiply.reduce(numValsByCol.shape))) # 1d indices of empty columns
                numEmptyCol = emptyColInd.shape[0]
                # fill empty columns in order to count number of rows with missing values
                dataFilledEmptyCol = MA.reshape(MA.array(self._chipdataMA), (self._chipdataMA.shape[0], self._chipdataMA.shape[1]*self._chipdataMA.shape[2]))
                for idx in emptyColInd:
                    dataFilledEmptyCol[:,idx] = 1e20
                dataFilledEmptyCol = MA.reshape(dataFilledEmptyCol, self._chipdataMA.shape)
                self.numRowsMissingChipData = 0
                for idx2 in range(dataFilledEmptyCol.shape[2]):
                    self.numRowsMissingChipData += Numeric.add.reduce(Numeric.where(MA.count(dataFilledEmptyCol[:,:,idx2], 1) < dataFilledEmptyCol.shape[1], 1, 0))
                s1 = ""
                s2 = ""
                if numEmptyCol > 0: s1 = "s"
                if self.numRowsMissingChipData > 0: s2 = "s"
                self.infod.setText("missing %i values, %i column%s completely, %i row%s partially" % (numTotalMissing, numEmptyCol, s1, self.numRowsMissingChipData, s2))
            else:
                self.infod.setText("")                
        else:
            self._chipdata = None
            self._chipdataMA = None
            self.infoc.setText("No structured data on input")
            self.infod.setText("")
            self.numRowsMissingChipData = 0
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
            if self.smooth and self._chipdataMA != None:
                steps += self._chipdataMA.shape[0] * self._chipdataMA.shape[2]
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
            self.send("Examples", chipstat.ma2orng_keepClassMetas(ma2d, self._data))
        else:
            self.send("Examples", None)

    def _sendchipdata(self, pbStep):
        if self._chipdataMA != None:
            ma3d = MA.array(self._chipdataMA)
            if self.impute:
                for idx in range(ma3d.shape[2]):
                    ma3d[:,:,idx] = chipimpute.kNNimputeMA(ma3d[:,:,idx], int(self.imputeK), callback=lambda: self.progressBarAdvance(pbStep))
            if self.smooth:
                ws = max(1.1, int(self.windowSize)) # fixes bug in statc.loess
                for idx in range(ma3d.shape[2]):
                    ma3d[:,:,idx] = chipimpute.loessMA(ma3d[:,:,idx], ws, 1, callback=lambda: self.progressBarAdvance(pbStep))
            chipdataNew = []    # [(dirname0, [et0, et1, ...]), ...]
            idxTotal = 0
            for (dirname, etList) in self._chipdata:
                etListNew = []
                for et in etList:
                    etListNew.append(chipstat.ma2orng_keepClassMetas(ma3d[:,:,idxTotal], et))
                    etListNew[-1].name = et.name
                    idxTotal += 1
                chipdataNew.append((dirname, etListNew))
            self.send("Structured Data", chipdataNew)
        else:
            self.send("Structured Data", None)


if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWImputeProfiles()
    a.setMainWidget(ow)
    ow.show()
    ow.data(orange.ExampleTable("meanExpr_ann_pkaC.tab"))
    a.exec_loop()
    ow.saveSettings()
