"""
<name>Impute & Loess Profiles</name>
<description>Imputation and non-parametric smoothing of expression profiles.</description>
<category>Genomics</category>
<icon>icons/Unknown.png</icon>
<priority>320</priority>
"""

import math
import Numeric, MA
from OWWidget import *
import OWGUI
import chipstat
import chipimpute


class OWImputeSmoothProfiles(OWWidget):
    settingsList  = ['impute', 'imputeK', 'smooth', 'windowSize', 'commitOnChange']

    def __init__(self, parent=None, name='Impute & Loess Profiles'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, name, "Imputation and non-parametric smoothing of expression profiles.")

        self._data = None       # exampleTable
        self._dataMA = None     # input 2d masked array
        self._dataOut = None    # output 2d masked array 
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
        self.infoa = QLabel("No data on input", box)
        self.infob = QLabel("", box)
        self.infoc = QLabel("", box)
        OWGUI.separator(self.controlArea)
        # KNN impute
        self.boxImpute = QVGroupBox("Impute missing values", self.controlArea)
        OWGUI.checkBox(self.boxImpute, self, "impute", "KNN impute", tooltip="Impute missing values from values of K nearest neighbours.", callback=self.change)
##        self.leImputeK = OWGUI.hSlider(self.boxImpute, self, "imputeK", box=None, minValue=1, maxValue=7744, step=1, callback=self.imputeChange, labelFormat=" K = %i", ticks=0)
##        self.leImputeK = OWGUI.qwtHSlider(self.boxImpute, self, "imputeK", box=None, label="K", labelWidth=12, minValue=1, maxValue=7744, step=0.02, precision=0, callback=self.imputeChange, logarithmic=1, ticks=0, maxWidth=200)
        self.leImputeK = OWGUI.qwtHSlider(self.boxImpute, self, "imputeK", box=None, label="K", labelWidth=12, minValue=1, maxValue=100, step=1, precision=0, callback=self.imputeKChange, logarithmic=0, ticks=0, maxWidth=130)
        self.boxImpute.setDisabled(1)
        OWGUI.separator(self.controlArea)
        # loess
        self.boxLoess = QVGroupBox("Smoothing", self.controlArea)
        OWGUI.checkBox(self.boxLoess, self, "smooth", "Loess smoothing", tooltip="Loess profiles, impute missing columns.", callback=self.change)
        lbl = QLabel("Window size (number of points)", self.boxLoess)
        lbl.setAlignment(Qt.AlignHCenter)
        self.leWindowSize = OWGUI.qwtHSlider(self.boxLoess, self, "windowSize", box=None, label="W", labelWidth=12, minValue=1, maxValue=10, step=1, precision=0, callback=self.smoothWChange, logarithmic=0, ticks=0, maxWidth=130)
        self.boxLoess.setDisabled(1)
        OWGUI.separator(self.controlArea)
        # output
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change', callback=self.change)
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

        self.inputs = [("Examples", ExampleTable, self.data, 1)]
        self.outputs = [("Examples", ExampleTable)]
        self.resize(50,50)

        # data dependent variables
        self.numRowsMissing = 0

    def change(self):
        if self.commitOnChange:
            self.senddata()

    def imputeKChange(self):
        if self.commitOnChange and self.impute:
            self.senddata()

    def smoothWChange(self):
        if self.commitOnChange and self.smooth:
            self.senddata()

    def data(self, data):
        if data != None:
            self._data = data
            self._dataMA = chipstat.orng2ma(data)
            # info text
            self.infoa.setText("Data: %i profiles on %i points" % (self._dataMA.shape[0], self._dataMA.shape[1]))
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
                self.infob.setText("Missing %i values," % (numTotalMissing))
                self.infoc.setText("%i column%s completely, %i row%s partially" % (numEmptyCol, s1, self.numRowsMissing, s2))
            # impute, loess
            self.boxImpute.setEnabled(1)
            maxK = min(max(self._dataMA.shape[0]/50, 50), self._dataMA.shape[0])
            self.leImputeK.setRange(1, maxK, 1)
            if self.imputeK > maxK:
                self.leImputeK.setValue(maxK)
                self.imputeK = maxK
            self.boxLoess.setEnabled(1)
            self.leWindowSize.setRange(1, self._dataMA.shape[1], 1)
            if self.windowSize > self._dataMA.shape[1]:
                self.leWindowSize.setValue(self._dataMA.shape[1])
                self.windowSize = self._dataMA.shape[1]
            # commit button
            self.commitBtn.setEnabled(1)
        else:
            self._data = None
            self._dataMA = None
            self.infoa.setText("No data on input")
            self.infob.setText("")
            self.infoc.setText("")
            self.boxImpute.setDisabled(1)
            self.boxLoess.setDisabled(1)
            self.commitBtn.setDisabled(1)
            self.numRowsMissing = 0
        if self.commitOnChange:
            self.senddata()

    def senddata(self):
        if self._dataMA != None:
            self._dataOut = self._dataMA
            self.progressBarInit()
            if self.numRowsMissing == 0:
                nrm = 1
            else:
                nrm = self.numRowsMissing
            if self.impute and self.smooth:
                pbStep1 = 45./nrm
                pbStep2 = 45./self._dataOut.shape[0]
                pbSet1 = 45
                pbSet2 = 90
            elif self.impute:
                pbStep1 = 90./nrm
                pbSet1 = 90
            elif self.smooth:
                pbStep2 = 90./self._dataOut.shape[0]
                pbSet2 = 90
            # impute
            if self.impute:
                self._dataOut = chipimpute.kNNimputeMA(self._dataOut, int(self.imputeK), callback=lambda: self.progressBarAdvance(pbStep1))
                self.progressBarSet(pbSet1)
            if self.smooth:
                ws = int(self.windowSize)
                if ws == 1: ws = 1.1    # fixes bug in statc.loess
                self._dataOut = chipimpute.loessMA(self._dataOut, ws, 1, callback=lambda: self.progressBarAdvance(pbStep2))
                self.progressBarSet(pbSet2)
            self.send("Examples", chipstat.ma2orng_keepClassMetas(self._dataOut, self._data))
            self.progressBarFinished()
        else:
            self.send("Examples", None)


if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWImputeSmoothProfiles()
    a.setMainWidget(ow)
    ow.show()
    ow.data(orange.ExampleTable("meanExpr_ann_pkaC.tab"))
    a.exec_loop()
    ow.saveSettings()
