"""
<name>Approximate Profiles</name>
<description>Approximation of expression profiles by various kernel functions.</description>
<category>Genomics</category>
<icon>icons/Unknown.png</icon>
<priority>310</priority>
"""

from OWWidget import *
import OWGUI
import chipstat
import chipappx


class OWApproxProfiles(OWWidget):
    #settingsList  = ['kernel', 'kernelSize', 'useSignificance', 'alpha', 'discardAbove', 'commitOnChange']
    settingsList  = ['kernel', 'kernelSize', 'useSignificance', 'alpha', 'commitOnChange']

    def __init__(self, parent=None, name='Approximate Profiles'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, name, "Approximation of expression profiles by various kernel functions.")

        self._data = None       # exampleTable
        self._dataMA = None     # 2d masked array
        self.kernel = 0
        self.kernels = ["Orthogonal polynomials", "Trigonometric functions"]
        self.kernelSize = None
        self.useSignificance = 0
        self.alpha = 3
        self.alphas = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5]
##        self.discardAbove = 0
        self.commitOnChange = 1
               
        # Settings
        self.loadSettings()

        # GUI: info, comboKernelSize, cbUseSignificance, vboxSignificance, commitBtn
        # info
        box = QVGroupBox("Info", self.controlArea)
        self.info = QLabel("No data on input", box)
        OWGUI.separator(self.controlArea)

        # kernel selection
        box = QVGroupBox("Kernel functions", self.controlArea)
        OWGUI.comboBox(box, self, "kernel", items=self.kernels, callback=self.kernelChange)
##        OWGUI.comboBox(box, self, "kernel", items=self.kernels, callback=self.selectionChange)
        OWGUI.separator(self.controlArea)

        # kernel settings
        box = QVGroupBox("Kernel settings", self.controlArea)
        self.comboKernelSize = OWGUI.comboBox(box, self, "kernelSize", sendSelectedValue=1, callback=self.selectionChange, label="Number of kernel functions", labelWidth=135, orientation="horizontal")
        self.comboKernelSize.setDisabled(1)
        self.cbUseSignificance = OWGUI.checkBox(box, self, "useSignificance", "Significance of coefficients", callback=self.useSignificanceChange, tooltip="Use kernels with coefficients significantly different from 0.")
        self.vboxSignificance = QVBox(box)
        OWGUI.comboBox(self.vboxSignificance, self, "alpha", items = self.alphas, callback=self.selectionChange, label="p <", labelWidth=20, orientation="horizontal")
##        OWGUI.checkBox(self.vboxSignificance, self, "discardAbove", "Discard kernel functions above rejected", callback=self.selectionChange)
        OWGUI.separator(self.controlArea)

        # output
        box = QVGroupBox("Output", self.controlArea)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.commitBtn = OWGUI.button(box, self, "Commit", callback=self.senddata, disabled=1)

        self.inputs = [("Examples", ExampleTable, self.data, 1)]
        self.outputs = [("Approximated Examples", ExampleTable), ("Approximation Coefficients", ExampleTable)]
        self.resize(200,100)

    def kernelChange(self):
        if self.kernel == 0:
            self.cbUseSignificance.setEnabled(1)
            self.vboxSignificance.setEnabled(1)
        elif self.kernel == 1:
            self.cbUseSignificance.setDisabled(1)
            self.vboxSignificance.setDisabled(1)
        if self.commitOnChange:
            self.senddata()

    def selectionChange(self):
        if self.commitOnChange:
            self.senddata()

    def useSignificanceChange(self):
        if self.useSignificance:
            self.vboxSignificance.setEnabled(1)
        else:
            self.vboxSignificance.setDisabled(1)
        if self.commitOnChange:
            self.senddata()

    def data(self, data):
        if data != None:
            self._data = data
            self._dataMA = chipstat.orng2ma(data)
            # info text
            self.info.setText("Data: %i profiles on %i points" % (self._dataMA.shape[0], self._dataMA.shape[1]))
            # combo kernel size
            self.comboKernelSize.clear()
            for i in range(0,self._dataMA.shape[1]):
                self.comboKernelSize.insertItem(str(i+1))
            self.comboKernelSize.setCurrentItem(max([int((self._dataMA.shape[1])/2)-2,0]))
            self.comboKernelSize.setEnabled(1)
            self.kernelSize = int(str(self.comboKernelSize.currentText()))
            # commit button
            self.commitBtn.setEnabled(1)
        else:
            self._data = None
            self._dataMA = None
            self.info.setText("No data on input")
            self.comboKernelSize.setDisabled(1)
            self.commitBtn.setDisabled(1)
        if self.commitOnChange:
            self.senddata()

    def senddata(self):
        if self._dataMA != None:
            # approximation
            self.progressBarInit()
            if self.kernel == 0:
                appx = chipappx.ApproxOrthPolyBasis(range(self._dataMA.shape[1]), int(self.kernelSize))
                self.progressBarSet(30)
                if self.useSignificance:
                    coef = appx.getAppxCoef2d_significant(self._dataMA, int(self.kernelSize), self.alphas[self.alpha])
                else:
                    coef = appx.getAppxCoef(self._dataMA)
                self.progressBarSet(60)
                curve = appx.getAppxCurve(coef)
                self.progressBarSet(90)
            elif self.kernel == 1:
                appx = chipappx.TrigonomerticBasis(self._dataMA.shape[1], int(self.kernelSize))
                self.progressBarSet(30)
##                if self.useSignificance:
##                    coef = appx.getAppxCoef2d_significant(self._dataMA, int(self.kernelSize), self.alphas[self.alpha])
##                else:
                coef = appx.getAppxCoef(self._dataMA)
                self.progressBarSet(60)
                curve = appx.getAppxCurve(coef)
                self.progressBarSet(90)
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
            self.progressBarFinished()
        else:
            self.send("Approximated Examples", None)
            self.send("Approximation Coefficients", None)
            


if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWApproxProfiles()
    a.setMainWidget(ow)
    ow.show()
    ow.data(orange.ExampleTable("meanExpr_ann_pkaC.tab"))
    a.exec_loop()
    ow.saveSettings()
