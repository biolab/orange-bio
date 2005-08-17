"""
<name>ANOVA</name>
<description>One/Two Way Analysis of Variance.</description>
<icon>icons/ChipANOVA.png</icon>
<priority>1070</priority>
"""

import Numeric, MA
from OWWidget import *
import OWGUI
import qwt
from OWDataFiles import DataFiles, ExampleSelection
import Anova


class OWANOVA(OWWidget):
    settingsList  = ["anovaType", "interaction", "selectorA", "selectorB", "selectorI", "alphaA", "alphaB", "alphaI", "autoUpdateSelName", "commitOnChange"]

    def __init__(self, parent=None, signalManager = None):
        print "__init__ S"
        OWWidget.__init__(self, parent, signalManager, 'ANOVA')
        # input / output data: [("name1", [orange.ExampleTable1a,...]), ("name2", [orange.ExampleTable2a,...])]
        self.inputs = [("Structured Data", DataFiles, self.onDataInput, 1)]
        self.outputs = [("Example Selection", ExampleSelection), ("Selected Structured Data", DataFiles), ("Other Structured Data", DataFiles)]

        # data, p-values, selected examples
        self.dataStructure = None   # input data
        self.numExamples = 0
        self.attrNameList = []      # names of attributes
        self.ps = None              # p-values: 2D Numeric.array of shape (3, numExamples)
        self.selectorName = ""      # for Example Selection output: (self.selectorName, [0,1,0,...])

        # Settings
        self.anovaType = 0      # 0: one-way (A), 1: one-way (B), 2: two-way
        self.interaction = 0    # 0: no interaction, 1: test for interaction effect
        self.selectorA = True
        self.selectorB = False
        self.selectorI = False
        self.alphaA = "0.05"
        self.alphaB = "0.05"
        self.alphaI = "0.05"
        self.autoUpdateSelName = 1
        self.commitOnChange = 0
        self.loadSettings()

        # GUI
        self.mainArea.setFixedWidth(0)
        ca=QFrame(self.controlArea)
        gl=QGridLayout(ca,4,1,5)
        
        # info
        box = QVGroupBox("Info", ca)
        gl.addWidget(box,0,0)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('', box)
        self.infoc = QLabel('', box)

        # ANOVA type
        self.boxAnovaType = QVButtonGroup("ANOVA Type", ca)
        gl.addWidget(self.boxAnovaType,1,0)
        self.boxAnovaType.setDisabled(1)
        self.rbgAnovaType = OWGUI.radioButtonsInBox(self.boxAnovaType, self, value="anovaType", btnLabels=["Single-factor (A, attributes)", "Single-factor (B, data sets)", "Two-factor"], callback=self.onAnovaType)
        self.cbInteraction = OWGUI.checkBox(self.boxAnovaType, self, value="interaction", label="Test interaction", callback=self.onInteraction)
        
        # selection of examples
        self.boxSelection = QVGroupBox("Example Selection", ca)
        gl.addWidget(self.boxSelection,2,0)
        self.lblNumGenes = []   # list of labels
        # selector A
        self.boxSelectorA = QVBox(self.boxSelection)
        self.cbSelectorA = OWGUI.checkBox(self.boxSelectorA, self, "selectorA", "Factor A (attributes)", callback=self.onSelectionChange,
                                          tooltip='H0: The mean does not depend on factor A (represented by attributes).')
        frmA = QFrame(self.boxSelectorA)
        glA = QGridLayout(frmA,1,3,5)
        leA = OWGUI.lineEdit(frmA, self, "alphaA", orientation="horizontal", callback=lambda x=0: self.onAlphaChange(x))
        glA.addWidget(leA,0,1) # Qt.AlignRight
##        le.setFixedSize(50, le.sizeHint().height())
##        self.connect(le, SIGNAL("clearFocus()"), self.onAlpha)
        glA.addWidget(QLabel("     p < ", frmA), 0,0)
        self.lblNumGenes.append(QLabel('', frmA))
        glA.addWidget(self.lblNumGenes[-1],0,2) # Qt.AlignRight | 0x22

        # selector B
        self.boxSelectorB = QVBox(self.boxSelection)
        self.cbSelectorB = OWGUI.checkBox(self.boxSelectorB, self, "selectorB", "Factor B (datasets)", callback=self.onSelectionChange,
                                          tooltip='H0: The mean does not depend on factor B (represented by datasets).')
        frmB = QFrame(self.boxSelectorB)
        glB = QGridLayout(frmB,1,3,5)
        leB = OWGUI.lineEdit(frmB, self, "alphaB", orientation="horizontal", callback=lambda x=1: self.onAlphaChange(x))
        glB.addWidget(leB,0,1)
        glB.addWidget(QLabel("     p < ", frmB), 0,0)
        self.lblNumGenes.append(QLabel('', frmB))
        glB.addWidget(self.lblNumGenes[-1],0,2)
        
        # selector I
        self.boxSelectorI = QVBox(self.boxSelection)
        self.cbSelectorI = OWGUI.checkBox(self.boxSelectorI, self, "selectorI", "Interaction (attributes * datasets)", callback=self.onSelectionChange,
                                          tooltip='H0: There is no interaction between factor A and factor B.')
        frmI = QFrame(self.boxSelectorI)
        glI = QGridLayout(frmI,1,3,5)
        leI = OWGUI.lineEdit(frmI, self, "alphaI", orientation="horizontal", callback=lambda x=2: self.onAlphaChange(x))
        ## slider could be used to replace lineEdit (but not sensitive enough)
        ##        self.alphaIf = 0.05
        ##        leI = OWGUI.qwtHSlider(self.boxSelectorI, self, "alphaIf", box="", label="      p < ", labelWidth=None, minValue=0.0001, maxValue=1.0, step=0.1, precision=3, callback=lambda x=2: self.onAlphaChange(x), logarithmic=1, ticks=0, maxWidth=None)
        glI.addWidget(leI,0,1)
        glI.addWidget(QLabel("     p < ", frmI), 0,0)
        self.lblNumGenes.append(QLabel('', frmI))
        glI.addWidget(self.lblNumGenes[-1],0,2)
        
        # output
        box = QVGroupBox("Output", ca)
        gl.addWidget(box,3,0)
        self.leSelectorName = OWGUI.lineEdit(box, self, 'selectorName', label='Selector Name: ')
        self.leSelectorName.setReadOnly(self.autoUpdateSelName)
        OWGUI.checkBox(box, self, 'autoUpdateSelName', 'Automatically update selector name', callback=self.onAutoUpdateSelNameChange)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change')
        self.btnCommit = OWGUI.button(box, self, "Commit", callback=self.onCommit)

        # enable/disable anova type box, example selection box, commit button, update the number of examples for individual selectors
        self.updateSettings()
        self.updateSelectorInfos()
        self.updateSelectorName()
        self.resize(350, self.sizeHint().height())
        print "__init__ F"


    def onDataInput(self, structuredData):
        """handles input data; sets self.dataStructure, self.numExamples, self.attrNameList and self.ps;
        updates info, calls updateSettings(), runs ANOVA and sends out new data.
        """
        print "onDataInput S"
        self.dataStructure = structuredData
        self.numExamples = 0
        self.attrNameList = []
        self.ps = None
        if structuredData:
            numFiles = reduce(lambda a,b: a+len(b[1]), structuredData, 0)
            lenSD = len(structuredData)
            self.infoa.setText("%d set%s, total of %d data file%s." % (lenSD, ["","s"][lenSD!=1], numFiles, ["","s"][numFiles!=1]))
            numExamplesList = []
            # construct a list of ExampleTable lengths and a list of attribute names
            for (name, etList) in structuredData:
                for et in etList:
                    numExamplesList.append(len(et))
                    for attr in et.domain.attributes:
                        if attr.name not in self.attrNameList:
                            self.attrNameList.append(attr.name)
            # test that all ExampleTables consist of equal number of examples
            if len(numExamplesList) == 0 or Numeric.add.reduce(Numeric.equal(numExamplesList, numExamplesList[0])) != len(numExamplesList):
                self.dataStructure = None
                self.numExamples = -1
                self.attrNameList = []
                self.infob.setText("Error: data files contain unequal number of examples, aborting ANOVA computation.")
                self.infoc.setText('')
            else:
                self.numExamples = numExamplesList[0]
                numAttributes = len(self.attrNameList)
                self.infob.setText("%d attribute%s, %d example%s in each file." % (numAttributes, ["","s"][numAttributes!=1], self.numExamples, ["","s"][self.numExamples!=1]))
                if self.numExamples > 0:
                    self.infoc.setText('Press Commit button to start ANOVA computation.')
                else:
                    self.infoc.setText('')
                self.boxAnovaType.setEnabled(1)
                self.boxSelection.setEnabled(1)
                self.btnCommit.setEnabled(True)
        else:
            self.infoa.setText('No data on input.')
            self.infob.setText('')
            self.infoc.setText('')
        # enable/disable anova type selection depending on the type of input data
        self.updateSettings()
        if self.autoUpdateSelName:
            self.updateSelectorName()
        # run ANOVA
        if self.commitOnChange:
            self.runANOVA()
            self.senddata()
        self.updateSelectorInfos()
        print "onDataInput F"

        
    def runANOVA(self):
        """converts structured data [(name, [orngET1, orngET2, ...]),...] to a 3D masked array
        with the following axes: 0: examples, 1: attributes, 2: ExampleTables;
        runs ANOVA computations and sets self.ps;
        """
        print "runAnova S"
        if self.dataStructure and self.numExamples > 0:
            ma3d = MA.zeros((self.numExamples, len(self.attrNameList), reduce(lambda a,b: a+len(b[1]), self.dataStructure, 0)), MA.Float) * MA.masked
            groupLens = []
            attrNameDict = dict(zip(self.attrNameList, range(len(self.attrNameList))))  # key: attrName, val: idx
            etIdx = 0
            for dsName, etList in self.dataStructure:
                for et in etList:
                    etm = et.toMA("a")[0]
                    if [attr.name for attr in et.domain.attributes] == self.attrNameList:
                        ma3d[:,:,etIdx] = etm
                    else:
                        print dsName + ": data copied for individual attributes"
                        for attrIdx, attr in enumerate(et.domain.attributes):
                            ma3d[:,attrNameDict[attr.name],etIdx] = etm[:,attrIdx]
                    etIdx += 1
                groupLens.append(len(etList))

            # run ANOVA
            self.infoc.setText('ANOVA computation started...')
            self.progressBarInit()
            pbStep = 100./self.numExamples
            self.ps = Numeric.ones((3, self.numExamples), Numeric.Float)
            if self.anovaType == 2:
                ps = self.anova2(ma3d, groupLens, self.interaction, repMeasuresOnA=False, callback=lambda: self.progressBarAdvance(pbStep))
                for rIdx in range(ps.shape[0]):
                    self.ps[rIdx] = ps[rIdx]
            elif self.anovaType == 1:
                self.ps[1] = self.anova1B(ma3d, groupLens, repMeasures=False, callback=lambda: self.progressBarAdvance(pbStep))
            elif self.anovaType == 0:
                self.ps[0] = self.anova1A(ma3d, repMeasures=False, callback=lambda: self.progressBarAdvance(pbStep))
            self.progressBarFinished()
        print "runAnova F"


    def anova1A(self, ma3d, repMeasures, callback):
        """conducts one-way ANOVA on individual examples wrt factor A (attributes, ma3d axis 1);
        returns Numeric array of p-values in shape (1, numExamples).
        """
        ps = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
        if repMeasures:
            fAnova = Anova.AnovaRM12LR
        else:
            fAnova = Anova.Anova1wayLR_2D
        for eIdx in range(ma3d.shape[0]):
            an = fAnova(MA.transpose(ma3d[eIdx]))
            ps[eIdx] = an.Fprob
            callback()
        return ps

    def anova1B(self, ma3d, groupLens, repMeasures, callback):
        """conducts one-way ANOVA on individual examples wrt factor B (datasets);
        ma3d axis 2 also contains replicas according to groupLens;
        returns Numeric array of p-values in shape (1, numExamples).
        WARNING: works slower than anova1A because it requires to copy 1D array to 2D array
                 although we could use Anova1wayLR instead of Anova1wayLR_2D, but not for repeated measures
                 additionaly, Anova1wayLR_2D handles missing factor levels correctly, which is not the case for Anova1wayLR
        """
        ps = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
        # groupLens [2,3,4] -> groupInd [[0,1],[2,3,4],[5,6,7,8]]
        if repMeasures:
            fAnova = Anova.AnovaRM12LR
        else:
            fAnova = Anova.Anova1wayLR_2D
        grpLensAcc = Numeric.concatenate([[0],Numeric.add.accumulate(groupLens)])
        grpInd = map(lambda i,j: range(i, j), grpLensAcc[:-1], grpLensAcc[1:])
        for eIdx in range(ma3d.shape[0]):
            m2 = MA.zeros((max(groupLens)*ma3d.shape[1], len(groupLens)), MA.Float) * MA.masked # axis0: replicas, axis1: factor B levels
            for groupIdx,takeInd in enumerate(grpInd):
                m2[:groupLens[groupIdx]*ma3d.shape[1], groupIdx] = MA.ravel(MA.take(ma3d[eIdx], takeInd, 1))
            an = fAnova(m2)
            ps[eIdx] = an.Fprob
            callback()
        return ps

    def anova2(self, ma3d, groupLens, addInteraction, repMeasuresOnA, callback):
        """Conducts two-way ANOVA on individual examples;
        returns a Numeric array of p-values in shape (2, numExamples) or (3, numExamples), depending whether we test for interaction;
        Note: levels of factors A and B that cause empty cells are removed prior to conducting ANOVA.
        """
        groupLens = Numeric.asarray(groupLens)
        # arrays to store p-vals
        if addInteraction:
            ps = Numeric.ones((3, ma3d.shape[0]), Numeric.Float)
        else:
            ps = Numeric.ones((2, ma3d.shape[0]), Numeric.Float)
        # decide between non-repeated / repeated measures ANOVA for factor time
        if repMeasuresOnA:
            fAnova = Anova.AnovaRM12LR
        else:
            fAnova = Anova.Anova2wayLR
        # check for empty cells for all genes at once and remove them
        tInd2rem = []
        ax2Ind = Numeric.concatenate(([0], Numeric.add.accumulate(groupLens)))
        for aIdx in range(ma3d.shape[1]):
            for rIdx in range(groupLens.shape[0]):
                if Numeric.add.reduce(MA.count(ma3d[:,aIdx,ax2Ind[rIdx]:ax2Ind[rIdx+1]],1)) == 0:
                    tInd2rem.append(aIdx)
                    break
        if len(tInd2rem) > 0:
            print "Warning: removing time indices %s for all genes" % (str(tInd2rem))
            tInd2keep = range(ma3d.shape[1])
            for aIdx in tInd2rem:
                tInd2keep.remove(aIdx)
            ma3d = MA.take(ma3d, tInd2keep, 1)
        # for each gene...
        for eIdx in range(ma3d.shape[0]):
            # faster check for empty cells for that gene -> remove time indices with empty cells
            ma2d = ma3d[eIdx]
            cellCount = MA.zeros((ma2d.shape[0], groupLens.shape[0]), MA.Int)
            for g,(i0,i1) in enumerate(zip(ax2Ind[:-1], ax2Ind[1:])):
                cellCount[:,g] = MA.count(ma2d[:,i0:i1], 1)
            ma2dTakeInd = Numeric.logical_not(Numeric.add.reduce(Numeric.equal(cellCount,0),1)) # 1 where to take, 0 where not to take
            if Numeric.add.reduce(ma2dTakeInd) != ma2dTakeInd.shape[0]:
                print "Warning: removing time indices %s for gene %i" % (str(Numeric.compress(ma2dTakeInd == 0, Numeric.arange(ma2dTakeInd.shape[0]))), eIdx)
                ma2d = MA.compress(ma2dTakeInd, ma2d, 0)
            # conduct ANOVA
            # 2    def __init__(self, arr2d, groupLens, addInteraction=0, allowReductA=True, allowReductB=False):
            # 2RM  def __init__(self, arr2d, groupLens, addInteraction=0, allowReductA=True, allowReductB=False):
            an = fAnova(ma2d, groupLens, addInteraction, allowReductA=True, allowReductB=True)
            ps[:,eIdx] = an.ps
            callback()
        return ps


    def updateSettings(self):
        """enables/disables: - anova type selection box;
                             - two-way anova and interaction checkbox;
                             - example selection box;
                             - factor A/B and interaction selectors
                             - Commit button
        """
        print "updateSettings S"
        if self.dataStructure and self.numExamples > 0:
            # enable anova type box, example selection box, commit button
            self.boxAnovaType.setEnabled(1)
            self.boxSelection.setEnabled(1)
            self.btnCommit.setEnabled(1)
            # enable/disable: two-way anova radio button, interaction checkbox
            if len(self.dataStructure) == 1 and len(self.attrNameList) > 1:
                # switch to single-factor (A) ANOVA, disable other ANOVAs
                self.anovaType = 0
                self.rbgAnovaType.buttons[1].setDisabled(1)
                self.rbgAnovaType.buttons[2].setDisabled(1)
                self.cbInteraction.setDisabled(1)
            elif len(self.dataStructure) > 1 and len(self.attrNameList) == 1:
                # switch to single-factor (B) ANOVA, disable other ANOVAs
                self.anovaType = 1
                self.rbgAnovaType.buttons[0].setDisabled(1)
                self.rbgAnovaType.buttons[2].setDisabled(1)
                self.cbInteraction.setDisabled(1)
            elif len(self.dataStructure) > 1 and len(self.attrNameList) > 1:
                # enable single and two-factor ANOVAs
                self.rbgAnovaType.buttons[0].setEnabled(1)
                self.rbgAnovaType.buttons[1].setEnabled(1)
                self.rbgAnovaType.buttons[2].setEnabled(1)
                self.cbInteraction.setEnabled(1)
            # enable/disable selectors
            if self.anovaType == 0:
                self.boxSelectorA.setEnabled(1)
                self.boxSelectorB.setDisabled(1)
                self.boxSelectorI.setDisabled(1)
            elif self.anovaType == 1:
                self.boxSelectorA.setDisabled(1)
                self.boxSelectorB.setEnabled(1)
                self.boxSelectorI.setDisabled(1)
            elif self.anovaType == 2:
                self.boxSelectorA.setEnabled(1)
                self.boxSelectorB.setEnabled(1)
                self.boxSelectorI.setEnabled(self.interaction)
        else:
            # disable anova type box, example selection box, commit button
            self.boxAnovaType.setDisabled(1)
            self.boxSelection.setDisabled(1)
            self.btnCommit.setDisabled(1)
        print "updateSettings F"


    def updateSelectorInfos(self, selectorIdx=None):
        """updates the number of examples that match individual selectors;
        if selectorIdx is given, updates only the corresponding info.
        """
        print "updateSelectorInfos S:" + str(selectorIdx)
        if not selectorIdx:
            selectorInd = range(3)
        else:
            selectorInd = [selectorIdx]
        alphas = [self.alphaA, self.alphaB, self.alphaI]
        for si in selectorInd:
            try:
                alpha = float(alphas[si])
                ps = self.ps[si]
            except:
                alpha = None
                ps = None
            if ps !=None and alpha != None and self.selectorA and self.anovaType in [[0,2],[1,2],[2]][si] and (self.interaction or [1,1,0][si]):
                numSelected = Numeric.add.reduce(Numeric.less(self.ps[si], alpha))
                self.lblNumGenes[si].setText('  (%d example%s)' % (numSelected, ['', 's'][numSelected!=1]))
            else:
                self.lblNumGenes[si].setText('  (no examples)')
        print "updateSelectorInfos F:" + str(selectorIdx)


    def senddata(self):
        """computes selectionList, partitions the examples and updates infoc;
        sends out selectionList and selected/other dataStructure or None;
        """
        print "senddata S"
        if self.dataStructure and self.ps:
            # set selectionList
            selectionList = Numeric.ones((self.numExamples,))
            print selectionList
            if self.selectorA and self.anovaType in [0,2]:
                try:
                    alpha = float(self.alphaA)
                except:
                    alpha = 0.0
                selectionList = Numeric.logical_and(selectionList, Numeric.less(self.ps[0], alpha))
                print selectionList
            if self.selectorB and self.anovaType in [1,2]:
                try:
                    alpha = float(self.alphaA)
                except:
                    alpha = 0.0
                selectionList = Numeric.logical_and(selectionList, Numeric.less(self.ps[1], alpha))
                print selectionList
            if self.selectorI and self.anovaType == 2 and self.interaction:
                try:
                    alpha = float(self.alphaA)
                except:
                    alpha = 0.0
                selectionList = Numeric.logical_and(selectionList, Numeric.less(self.ps[2], alpha))
                print selectionList
            numExamples = Numeric.add.reduce(Numeric.greater(selectionList, 0))
            self.infoc.setText('Total of %d example%s match criteria.' % (numExamples, ['', 's'][numExamples!=1]))
            selectionList = selectionList.tolist()
            self.send("Example Selection", (self.selectorName, selectionList))
            # partition dataStructure
            dataStructS = []
            dataStructN = []
            for (dsName, etList) in self.dataStructure:
                etListS = [et.select(selectionList) for et in etList]
                etListN = [et.select(selectionList, negate=1) for et in etList]
                for i in range(len(etList)):
                    etListS[i].name = etListN[i].name = etList[i].name
                dataStructS.append((dsName, etListS))
                dataStructN.append((dsName, etListN))
            self.send("Selected Structured Data", dataStructS)
            self.send("Other Structured Data", dataStructN)
        else:
            self.send("Example Selection", None)
            self.send("Selected Structured Data", None)
            self.send("Other Structured Data", None)
        print "senddata F"
            

    def updateSelectorName(self):
        """update selector name shown in selector edit box
        """
        s = 'ANOVA'
        if self.dataStructure:
            s += " (%s)" % reduce(lambda a,b: a + ", " + b[0], self.dataStructure, "")[2:]
            if self.selectorA and self.anovaType in [0,2]:
                s += ", pA<%s" % self.alphaA
            if self.selectorB and self.anovaType in [1,2]:
                s += ", pB<%s" % self.alphaB
            if self.selectorI and self.anovaType == 2 and self.interaction:
                s += ", pI<%s" % self.alphaI
        self.selectorName = s.strip()


    #==========================================================================
    # Event handlers
    #==========================================================================

    def onAnovaType(self):
        """handles changes of ANOVA type:
            - resets self.ps;
            - calls updateSettings()
            - updates infoc
            - calls updateSelectorInfos()
        runs ANOVA and sends out new data;
        """
        print "onAnovaType S"
        self.ps = None
        self.updateSettings()
        if self.autoUpdateSelName:
            self.updateSelectorName()
        if self.commitOnChange:
            self.runANOVA()
            self.senddata()
        elif self.dataStructure and self.numExamples > 0:
            self.infoc.setText('Press Commit button to start ANOVA computation.')
        self.updateSelectorInfos()
        print "onAnovaType F"


    def onInteraction(self):
        """handles clicks on interaction checkbox:
            - resets self.ps;
            - enables/disables selector I,
            - updates infoc
        runs ANOVA and sends out new data
        """
        print "onInteraction S"
        self.ps = None
        self.boxSelectorI.setEnabled(self.interaction)
        if self.autoUpdateSelName:
            self.updateSelectorName()
        if self.commitOnChange:
            self.runANOVA()
            self.senddata()
        elif self.dataStructure and self.numExamples > 0:
            self.infoc.setText('Press Commit button to start ANOVA computation.')
        self.updateSelectorInfos()
        print "onInteraction F"


    def onSelectionChange(self):
        """handles changes in example selector checkboxes;
        sends out new data;
        """
        print "onSelectionChange S"
        if self.autoUpdateSelName:
            self.updateSelectorName()
        if self.commitOnChange:
            self.senddata()
        print "onSelectionChange F"


    def onAlphaChange(self, selectorIdx):
        """handles changes in example selector alphas;
        prints number of selected examples for individual selectors and sends out new data;
        """
        print "onAlphaChange S"
        if self.autoUpdateSelName:
            self.updateSelectorName()
        if self.commitOnChange:
            self.senddata()
        self.updateSelectorInfos(selectorIdx)
        print "onAlphaChange F"


    def onAutoUpdateSelNameChange(self):
        """handles clicks on auto update selector name checkbox
        """
        self.leSelectorName.setReadOnly(self.autoUpdateSelName)

            
    def onCommit(self):
        """handles Commit clicks; runs ANOVA (if not already computed) and sends out data;
        """
        print "onCommit S"
        if self.dataStructure:
            if not self.ps:
                self.runANOVA()
            self.senddata()
        self.updateSelectorInfos()
        print "onCommit F"
        print self.height(), self.width()


if __name__=="__main__":
    import OWDataFiles, orngSignalManager
    signalManager = orngSignalManager.SignalManager(0)
    a=QApplication(sys.argv)
    ow=OWANOVA(signalManager = signalManager)
    a.setMainWidget(ow)
    ow.show()
    ds = OWDataFiles.OWDataFiles(signalManager = signalManager)
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS_10_yakApufA")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS_10_yakApufA_time0")
    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS_10_yakApufA_time0_swappedAB")
    signalManager.addWidget(ow)
    signalManager.addWidget(ds)
    signalManager.setFreeze(1)
    signalManager.addLink(ds, ow, 'Structured Data', 'Structured Data', 1)
    signalManager.setFreeze(0)
    a.exec_loop()
    ow.saveSettings()

### test: compare 1-way and 2-way ANOVA
##    d1d = Numeric.array([5,4,6, 4,7,6, 3,4,5])
##    rgi = [[0,1,2],[3,4,5],[6,7,8]]
##    a1 = Anova.Anova1wayLR(d1d, rgi)
##
##    d2d = Numeric.array([[5,4,6],[4,7,6],[3,4,5]])
##    gl = [3]
##    a2 = Anova.Anova2wayLR(d2d, gl, 1, False, False)
