## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

"""
<name>Approximate Profiles</name>
<description>Approximation of expression profiles by various kernel functions.</description>
<icon>icons/ApproximateProfiles.png</icon>
<priority>310</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

from __future__ import absolute_import

import numpy.oldnumeric as Numeric

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import chipappx, chipstat
from .OWDataFiles import DataFiles

class OWApproxProfiles(OWWidget):
    settingsList  = ['kernel', 'kernelSize', 'useSignificance', 'alpha', 'commitOnChange']

    def __init__(self, parent=None, signalManager = None, name='Approximate Profiles'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, signalManager, name)

        self._data = None       # exampleTable
        self._dataN = None      # 2d numeric array
        self._chipdata = None       # [(dirname0, [et0, et1, ...]), ...]
        self._chipdataN = None      # 3d numeric array
        self.kernel = 0
        self.kernels = ["Orthogonal polynomials", "Trigonometric functions"]
        self.kernelSize = None
        self.kernelSizes = [    ["Const"] + map(lambda x: "degree <= %i"%x, range(1,99)),
                                ["Const", "cos x", "sin x"] + reduce(lambda x,y: x+[y[0],y[1]], map(lambda i: ("cos %ix"%i, "sin %ix"%i), range(2,99)), [])
                            ]
        self.useSignificance = 0
        self.alpha = 3
        self.alphas = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
        self.commitOnChange = 1
               
        # Settings
        self.loadSettings()

        # GUI: info, comboKernelSize, cbUseSignificance, vboxSignificance, commitBtn
        # info
        box = QVGroupBox("Info", self.controlArea)
        self.infoa = QLabel("No examples on input", box)
        OWGUI.separator(box,250)
        self.infob = QLabel("No structured data on input", box)
        OWGUI.separator(self.controlArea)

        # kernel selection
        box = QVGroupBox("Kernel functions", self.controlArea)
        OWGUI.comboBox(box, self, "kernel", items=self.kernels, callback=self.kernelChange)
        OWGUI.separator(self.controlArea)

        # kernel settings
        box = QVGroupBox("Kernel settings", self.controlArea)
        self.comboKernelSize = OWGUI.comboBox(box, self, "kernelSize", callback=self.kernelSizeChange, label="Number of kernel functions", labelWidth=135, orientation="horizontal", valueType=int)
        self.comboKernelSize.setDisabled(1)
        self.cbUseSignificance = OWGUI.checkBox(box, self, "useSignificance", "Significance of coefficients (F-statistics)", callback=self.useSignificanceChange, tooltip="Use kernels with coefficients significantly different from 0.")
        self.vboxSignificance = QVBox(box)
        OWGUI.comboBox(self.vboxSignificance, self, "alpha", items = self.alphas, callback=self.alphaChange, label="p <", labelWidth=20, orientation="horizontal")
        OWGUI.separator(self.controlArea)

        # output
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

        self.inputs = [("Examples", ExampleTable, self.data), ("Structured Data", DataFiles, self.chipdata)]
        self.outputs = [("Approximated Examples", ExampleTable, Default), ("Approximation Coefficients", ExampleTable), ("Approximated Structured Data", DataFiles, Default), ("Structured Approximation Coefficients", DataFiles)]
        self.resize(200,100)

    def kernelChange(self):
        numPoints = self._getMinDataShape1()
        if self.kernel == 0 and self.kernelSize < numPoints-1:
            self.cbUseSignificance.setEnabled(1)
            self.vboxSignificance.setEnabled(1)
        else:
            self.cbUseSignificance.setDisabled(1)
            self.vboxSignificance.setDisabled(1)
        self.setGuiCommonExpChip()
        if self.commitOnChange:
            self.senddata()

    def kernelSizeChange(self):
        numPoints = self._getMinDataShape1()
        if self.kernel == 0:
            if self.kernelSize in [0, numPoints-1]:
                self.cbUseSignificance.setDisabled(1)
                self.vboxSignificance.setDisabled(1)
            else:
                self.cbUseSignificance.setEnabled(1)
                self.vboxSignificance.setEnabled(1)
        if self.commitOnChange:
            self.senddata()

    def alphaChange(self):
        if self.commitOnChange:
            self.senddata()

    def useSignificanceChange(self):
        if self.useSignificance:
            self.vboxSignificance.setEnabled(1)
        else:
            self.vboxSignificance.setDisabled(1)
        if self.commitOnChange:
            self.senddata()

    def _getMinDataShape1(self):            
        numPoints = 1e9 # number of attributes (if both inputs present: take min.)
        if self._dataN != None: numPoints = min(numPoints, self._dataN.shape[1])
        if self._chipdataN != None: numPoints = min(numPoints, self._chipdataN.shape[1])
        if numPoints == 1e9: numPoints = 0
        return numPoints

    def setGuiCommonExpChip(self):
        """Sets GUI features that depend on both inputs, e.g. max. values of sliders, commit button.
        """
        numPoints = self._getMinDataShape1()
        # set up comboKernelSize
        self.comboKernelSize.clear()
        # print "numPoints:", numPoints
        if numPoints > 0:
            for i in range(0,numPoints):
                self.comboKernelSize.insertItem(self.kernelSizes[self.kernel][i])
            if self.kernelSize==None or self.kernelSize >= numPoints:
                self.kernelSize = max([int((numPoints)/2)-2,0])
            self.comboKernelSize.setCurrentItem(self.kernelSize)
            self.comboKernelSize.setEnabled(1)
            self.commitBtn.setEnabled(1)
        else:
            self.comboKernelSize.setDisabled(1)
            #self.kernelSize = None
            self.commitBtn.setDisabled(1)

    def data(self, data):
        if data != None:
            self._data = data
            self._dataN = chipstat.orng2ma(data)
            self.infoa.setText("Examples: %i profiles on %i points" % (self._dataN.shape[0], self._dataN.shape[1]))
        else:
            self._data = None
            self._dataN = None
            self.infoa.setText("No examples on input")
        self.setGuiCommonExpChip()
        if self.commitOnChange:
            self.senddata(1)

    def chipdata(self, data):
        if data != None:
            self._chipdata = data
            shp = [0,0,0]
            shp[0] = len(data[0][1][0])
            shp[1] = len(data[0][1][0].domain.attributes)
            shp[2] = reduce(lambda x,y: x+len(y[1]), data, 0)
            self._chipdataN = Numeric.zeros(shp, Numeric.Float)
            idx = 0
            for (name, etList) in data:
                for et in etList:
                    self._chipdataN[:,:,idx] = Numeric.asarray(chipstat.orng2ma(et))
                    idx += 1
            self.infob.setText("Structured Data: %i data files with %i profiles on %i points" % (shp[2], shp[0], shp[1]))
        else:
            self._chipdata = None
            self._chipdataN = None
            self.infob.setText("No structured data on input")
        self.setGuiCommonExpChip()
        if self.commitOnChange:
            self.senddata(2)

    def senddata(self, outputSelector=0):
        """outputSelector = [0: examples + chip data | 1: examples | 2: chip data]
        """
        assert outputSelector in [0,1,2]
        # progress bar settings: 1: examples, 2: chip data
        steps = 0
        if outputSelector in [0,1] and self._dataN != None:
            steps += 3                                  # 1 step for approximator, 1 step for coefficients, 1 step for curves
        if outputSelector in [0,2] and self._chipdataN != None:
            steps += (1 + 2*self._chipdataN.shape[2])   # 1 step for approximator, for each dataset 2 steps for coeficients and curves
        if steps == 0: steps = 1
        pbStep = 100./steps
        self.progressBarInit()
        if outputSelector in [0,1]:
            self._sendexampledata(pbStep)
        if outputSelector in [0,2]:
            self._sendchipdata(pbStep)
        self.progressBarFinished()

    def _sendexampledata(self, pbStep):
        if self._dataN != None:
            # approximation
            if self.kernel == 0:
                appx = chipappx.ApproxOrthPolyBasis(range(self._dataN.shape[1]), self.kernelSize+1)
                self.progressBarAdvance(pbStep)
                if self.useSignificance and self.cbUseSignificance.isEnabled():
                    coef = appx.getAppxCoef2d_significant(self._dataN, self.kernelSize+1, self.alphas[self.alpha])
                else:
                    coef = appx.getAppxCoef(self._dataN)
                self.progressBarAdvance(pbStep)
                curve = appx.getAppxCurve(coef)
                self.progressBarAdvance(pbStep)
            elif self.kernel == 1:
                appx = chipappx.TrigonomerticBasis(self._dataN.shape[1], self.kernelSize+1)
                self.progressBarAdvance(pbStep)
                # 2007-10-11: trigonometric functions do not use getAppxCoef2d_significant
                #if self.useSignificance and self.cbUseSignificance.isEnabled():
                #    coef = appx.getAppxCoef2d_significant(self._dataN, self.kernelSize+1, self.alphas[self.alpha])
                #else:
                coef = appx.getAppxCoef(self._dataN)
##                print self._dataN
##                print "coef", coef
                self.progressBarAdvance(pbStep)
                curve = appx.getAppxCurve(coef)
##                print "curve", curve
                self.progressBarAdvance(pbStep)
            # domain with class and metas
            if self._data.domain.classVar != None:
                domainClassMetas = orange.Domain([self._data.domain.classVar])
            else:
                domainClassMetas = orange.Domain([])
            domainClassMetas.addmetas(self._data.domain.getmetas())
            # exampleTable with class and metas
            etClassMetas = orange.ExampleTable(domainClassMetas, self._data)
            # appx. curve
            etCurve = orange.ExampleTable(orange.Domain(self._data.domain.attributes, None), curve.tolist())
            etCurve = orange.ExampleTable([etCurve, etClassMetas])
            self.send("Approximated Examples", etCurve)
            # appx. coefficients
            domainCoef = orange.Domain(map(lambda x: orange.FloatVariable("C%i" % x), range(coef.shape[1])), None)
            etCoef = orange.ExampleTable(domainCoef, coef.tolist())
            etCoef = orange.ExampleTable([etCoef, etClassMetas])
            self.send("Approximation Coefficients", etCoef)
        else:
            self.send("Approximated Examples", None)
            self.send("Approximation Coefficients", None)

    def _sendchipdata(self, pbStep):
        if self._chipdataN != None:
            # approximation
            coefs = Numeric.zeros((self._chipdataN.shape[0], self.kernelSize+1, self._chipdataN.shape[2]), Numeric.Float)
            curves = Numeric.zeros(self._chipdataN.shape, Numeric.Float)
            if self.kernel == 0:
                appx = chipappx.ApproxOrthPolyBasis(range(self._chipdataN.shape[1]), self.kernelSize+1)
                self.progressBarAdvance(pbStep)
                for idx2 in range(self._chipdataN.shape[2]):
                    if self.useSignificance and self.cbUseSignificance.isEnabled():
                        coefs[:,:,idx2] = appx.getAppxCoef2d_significant(self._chipdataN[:,:,idx2], self.kernelSize+1, self.alphas[self.alpha])
                    else:
                        coefs[:,:,idx2] = appx.getAppxCoef(self._chipdataN[:,:,idx2])
                    self.progressBarAdvance(pbStep)
                    curves[:,:,idx2] = appx.getAppxCurve(coefs[:,:,idx2])
                    self.progressBarAdvance(pbStep)
            elif self.kernel == 1:
                appx = chipappx.TrigonomerticBasis(self._chipdataN.shape[1], self.kernelSize+1)
                self.progressBarAdvance(pbStep)
                for idx2 in range(self._chipdataN.shape[2]):
                    # 2007-10-11: trigonometric functions do not use getAppxCoef2d_significant
                    #if self.useSignificance and self.cbUseSignificance.isEnabled():
                    #    coefs[:,:,idx2] = appx.getAppxCoef2d_significant(self._chipdataN[:,:,idx2], self.kernelSize+1, self.alphas[self.alpha])
                    #else:
                    coefs[:,:,idx2] = appx.getAppxCoef(self._chipdataN[:,:,idx2])
                    self.progressBarAdvance(pbStep)
                    curves[:,:,idx2] = appx.getAppxCurve(coefs[:,:,idx2])
                    self.progressBarAdvance(pbStep)
            chipcoefNew = []    # [(dirname0, [etCoef0, etCoef1, ...]), ...]
            chipcurvesNew = []  # [(dirname0, [etCurves0, etCurves1, ...]), ...]
            idxTotal = 0
            for (dirname, etList) in self._chipdata:
                etCoefListNew = []
                etCurvesListNew = []
                for et in etList:
                    # domain with class and metas
                    if et.domain.classVar != None:
                        domainClassMetas = orange.Domain([et.domain.classVar])
                    else:
                        domainClassMetas = orange.Domain([])
                    domainClassMetas.addmetas(et.domain.getmetas())
                    # exampleTable with class and metas
                    etClassMetas = orange.ExampleTable(domainClassMetas, et)
                    # appx. coefficients
                    domainCoef = orange.Domain(map(lambda x: orange.FloatVariable("C%i" % x), range(coefs.shape[1])), None)
                    etCoef = orange.ExampleTable(domainCoef, coefs[:,:,idxTotal].tolist())
                    etCoefListNew.append(orange.ExampleTable([etCoef, etClassMetas]))
                    etCoefListNew[-1].name = et.name
                    # appx. curve
                    etCurve = orange.ExampleTable(orange.Domain(et.domain.attributes, None), curves[:,:,idxTotal].tolist())
                    etCurvesListNew.append(orange.ExampleTable([etCurve, etClassMetas]))
                    etCurvesListNew[-1].name = et.name
                    idxTotal += 1
                chipcoefNew.append((dirname, etCoefListNew))
                chipcurvesNew.append((dirname, etCurvesListNew))
            self.send("Approximated Structured Data", chipcurvesNew)
            self.send("Structured Approximation Coefficients", chipcoefNew)
        else:
            self.send("Approximated Structured Data", None)
            self.send("Structured Approximation Coefficients", None)


if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWApproxProfiles()
    a.setMainWidget(ow)
    ow.show()
    ow.data(orange.ExampleTable("meanExpr_ann_pkaC.tab"))
    a.exec_loop()
    ow.saveSettings()
