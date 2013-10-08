## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

"""
<name>ANOVA</name>
<description>Single Sample T-test, One/Two Way Analysis of Variance.</description>
<icon>icons/ChipANOVA.png</icon>
<priority>1070</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

from __future__ import absolute_import

import numpy.oldnumeric as Numeric, numpy.oldnumeric.ma as MA
import scipy.stats

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import Anova
from .OWDataFiles import DataFiles, ExampleSelection

class OWANOVA(OWWidget):
    settingsList  = ["anovaType", "compareToValue", "_interaction", "selectorA", "selectorB", "selectorI", "alphaA", "alphaB", "alphaI", "autoUpdateSelName", "sendNotSelectedData", "sendProbabilities", "commitOnChange"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'ANOVA')
        # input / output data: [("name1", [orange.ExampleTable1a,...]), ("name2", [orange.ExampleTable2a,...])]
        self.inputs = [("Structured Data", DataFiles, self.onDataInput)]
        self.outputs = [("Example Selection", ExampleSelection, Default), ("Selected Structured Data", DataFiles, Default), ("Other Structured Data", DataFiles)]

        # data, p-values, selected examples
        self.dataStructure = None   # input data
        self.numExamples = 0
        self.numVariables = 0
        self.ps = None              # p-values: 2D Numeric.array of shape (3, numExamples)
        self.selectorName = ""      # for Example Selection output: (self.selectorName, [0,1,0,...])

        # Settings
        self.anovaType = 0      # 0: single-sample t-test, 1: one-way (A), 2: one-way (B), 3: two-way (A,B), 4: ful factorial (A, B, A*B)
        self.compareToValue = 0 # single sample t-test, value to compare to
        self._interaction = 0    # 0: no interaction, 1: test for interaction effect (set this value manually !!!)
        self.selectorA = True
        self.selectorB = False
        self.selectorI = False
        self.alphaA = "0.05"
        self.alphaB = "0.05"
        self.alphaI = "0.05"
        self.autoUpdateSelName = 1
        self.sendNotSelectedData = 1
        self.sendProbabilities = 0
        self.commitOnChange = 0
        self.loadSettings()

        # GUI
        self.mainArea.setFixedWidth(0)
        ca = self.controlArea
       
        # info
        box = OWGUI.widgetBox(ca, "Info")
        #gl.addWidget(box,0,0)
        self.infoa = OWGUI.label(box, self, 'No data on input.')
        self.infob = OWGUI.label(box, self, "")
        self.infoc = OWGUI.label(box, self, "")

        # ANOVA type
        # group selection
        anovaTypes = ["Single sample t-test", 
            "Single-factor (A, variables)", 
            "Single-factor (B, data sets)", 
            "Two-factor", 
            "Two-factor with interaction effect"]

        self.boxAnovaType = OWGUI.widgetBox(ca, "Anova Type")
        self.anovaTypeS = OWGUI.radioButtonsInBox(self.boxAnovaType, self, "anovaType", btnLabels=anovaTypes)
        
        self.boxAnovaType.setDisabled(1)

        self.boxCompareTo = OWGUI.widgetBox(self.boxAnovaType)
        OWGUI.lineEdit(self.boxCompareTo, self, "compareToValue", callback=self.onCompareToChange, label="compare to")

        # selection of examples
        self.boxSelection = OWGUI.widgetBox(ca, "Example Selection")
        self.lblNumGenes = []   # list of labels

        # selector A
        self.boxSelectorA = OWGUI.widgetBox(self.boxSelection)
        self.cbSelectorA = OWGUI.checkBox(self.boxSelectorA, self, "selectorA", "Factor A (variables)", callback=self.onSelectionChange,
                                          tooltip='H0: The mean does not depend on factor A (represented by variables).')

        frmA = OWGUI.widgetBox(self.boxSelectorA)
        leA = OWGUI.lineEdit(frmA, self, "alphaA", orientation="horizontal", callback=lambda x=0: self.onAlphaChange(x), label= "p <= ")
        self.lblNumGenes.append(OWGUI.label(frmA, self, ""))

        # selector B
        self.boxSelectorB = OWGUI.widgetBox(self.boxSelection)
        self.cbSelectorB = OWGUI.checkBox(self.boxSelectorB, self, "selectorB", "Factor B (data sets)", callback=self.onSelectionChange,
                                          tooltip='H0: The mean does not depend on factor B (represented by data sets).')
 
        frmB = OWGUI.widgetBox(self.boxSelectorB)
        leB = OWGUI.lineEdit(frmB, self, "alphaB", orientation="horizontal", callback=lambda x=1: self.onAlphaChange(x), label= "p <= ")
        self.lblNumGenes.append(OWGUI.label(frmB, self, ""))

        # selector I
        self.boxSelectorI = OWGUI.widgetBox(self.boxSelection)
        self.cbSelectorI = OWGUI.checkBox(self.boxSelectorI, self, "selectorI", "Interaction (variables * data sets)", callback=self.onSelectionChange,
                                          tooltip='H0: There is no interaction between factor A and factor B.')
 
        frmI = OWGUI.widgetBox(self.boxSelectorI)
        leI = OWGUI.lineEdit(frmI, self, "alphaI", orientation="horizontal", callback=lambda x=2: self.onAlphaChange(x), label= "p <= ")
        self.lblNumGenes.append(OWGUI.label(frmI, self, ""))

        # output
        box = OWGUI.widgetBox(ca, "Output")
        self.leSelectorName = OWGUI.lineEdit(box, self, 'selectorName', label='Selector Name: ')
        self.leSelectorName.setReadOnly(self.autoUpdateSelName)
        OWGUI.checkBox(box, self, 'autoUpdateSelName', 'Automatically update selector name', callback=self.onAutoUpdateSelNameChange)
        OWGUI.checkBox(box, self, 'sendNotSelectedData', 'Send not selected data', callback=self.onSendNotSelectedChange)
        OWGUI.checkBox(box, self, 'sendProbabilities', 'Show p-values', callback=self.onSendProbabilitiesChange)
        OWGUI.checkBox(box, self, 'commitOnChange', 'Commit data on selection change', callback=lambda: self.onCommit(self.commitOnChange))
        self.btnCommit = OWGUI.button(box, self, "Commit", callback=self.onCommit)

        # enable/disable anova type box, example selection box, commit button, update the number of examples for individual selectors
        self.updateAnovaTypeBox()
        self.updateSelectorBox()
        self.updateSelectorInfos()
        self.updateSelectorName()

        self.resize(283, self.sizeHint().height())


    def onDataInput(self, structuredData):
        """handles input data; sets self.dataStructure, self.numExamples, self.numVariables and self.ps;
        updates info, calls updateAnovaTypeBox(), runs ANOVA and sends out new data.
        """
        self.dataStructure = structuredData
        self.numExamples = 0
        self.numVariables = 0
        self.ps = None
        if structuredData:
            numFiles = reduce(lambda a,b: a+len(b[1]), structuredData, 0)
            lenSD = len(structuredData)
            self.infoa.setText("%d set%s, total of %d data file%s." % (lenSD, ["","s"][lenSD!=1], numFiles, ["","s"][numFiles!=1]))
            numExamplesList = []
            numVariablesList = []
            # construct a list of ExampleTable lengths and a list of number of variables
            for (name, etList) in structuredData:
                for et in etList:
                    numExamplesList.append(len(et))
                    numVariablesList.append(len(et.domain.variables))
            # test that all ExampleTables consist of equal number of examples and variables
            if len(numExamplesList) == 0 or Numeric.add.reduce(Numeric.equal(numExamplesList, numExamplesList[0])) != len(numExamplesList):
                self.dataStructure = None
                self.numExamples = -1
                self.infob.setText("Error: data files contain unequal number of examples, aborting ANOVA computation.")
                self.infoc.setText('')
            elif len(numVariablesList) == 0 or Numeric.add.reduce(Numeric.equal(numVariablesList, numVariablesList[0])) != len(numVariablesList):
                self.dataStructure = None
                self.numVariables = -1
                self.infob.setText("Error: data files contain unequal number of variables, aborting ANOVA computation.")
                self.infoc.setText('')
            else:
                self.numExamples = numExamplesList[0]
                self.numVariables = numVariablesList[0]
                self.infob.setText("%d variable%s, %d example%s in each file." % (self.numVariables, ["","s"][self.numVariables!=1], self.numExamples, ["","s"][self.numExamples!=1]))
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
        self.updateAnovaTypeBox()
        self.updateSelectorBox()
        if self.autoUpdateSelName:
            self.updateSelectorName()
        # run ANOVA
        if self.commitOnChange:
            self.runANOVA()
            self.senddata()
        self.updateSelectorInfos()

        
    def runANOVA(self):
        """converts structured data [(name, [orngET1, orngET2, ...]),...] to a 3D masked array
        with the following axes: 0: examples, 1: variables, 2: ExampleTables;
        runs ANOVA computations and sets self.ps;
        """
        if self.dataStructure and self.numExamples > 0:
            ma3d = MA.zeros((self.numExamples, self.numVariables, reduce(lambda a,b: a+len(b[1]), self.dataStructure, 0)), Numeric.Float) * MA.masked
            groupLens = []
            etIdx = 0
            for dsName, etList in self.dataStructure:
                for et in etList:
                    ma3d[:,:,etIdx] = et.toNumpyMA("ac")[0]
                    etIdx += 1
                groupLens.append(len(etList))

            #print "ma3d SHAPE", ma3d.shape
            #print "ma3d from top", ma3d[0,:,:]
            # run ANOVA
            self.infoc.setText('ANOVA computation started...')
            self.progressBarInit()
            pbStep = 100./self.numExamples
            self.ps = Numeric.ones((3, self.numExamples), Numeric.Float)
            if self.anovaType >= 3:
                ps = self.anova2(ma3d, groupLens, self.anovaType==4, repMeasuresOnA=False, callback=lambda: self.progressBarAdvance(pbStep))
                for rIdx in range(ps.shape[0]):
                    self.ps[rIdx] = ps[rIdx]
            elif self.anovaType == 2:
                self.ps[1] = self.anova1B(ma3d, groupLens, repMeasures=False, callback=lambda: self.progressBarAdvance(pbStep))
            elif self.anovaType == 1:
                self.ps[0] = self.anova1A(ma3d, repMeasures=False, callback=lambda: self.progressBarAdvance(pbStep))
            elif self.anovaType == 0:
                try:
                    compToVal = float(self.compareToValue)
                except:
                    print "Warning: cannot convert %s to float, using 0" % str(self.compareToValue)
                    self.compareToValue = 0
                    compToVal = 0
                self.ps[0] = self.ttest_ssmpl(ma3d, compToVal, callback=lambda: self.progressBarAdvance(pbStep))
            self.progressBarFinished()


    def ttest_ssmpl(self, ma3d, compToVal, callback):
        """conducts single-sample t-test on individual examples wrt factor A (variables, ma3d axis 1);
        returns Numeric array of p-values in shape (1, numExamples).
        """
        ps = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
        for eIdx in range(ma3d.shape[0]):
            data = Numeric.asarray(MA.transpose(ma3d[eIdx]).compressed())
            if len(data) >= 2:
                try:
                    ps[eIdx] = scipy.stats.ttest_1samp(data, compToVal)[1]
                except:
                    print "Warning: zero variance, check the example %i:" % eIdx, data
                    ps[eIdx] = 1.0
            else:
                ps[eIdx] = 1.0
            callback()
        return ps

    def anova1A(self, ma3d, repMeasures, callback):
        """conducts one-way ANOVA on individual examples wrt factor A (variables, ma3d axis 1);
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
        """conducts one-way ANOVA on individual examples wrt factor B (data sets);
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
            m2 = MA.zeros((max(groupLens)*ma3d.shape[1], len(groupLens)), Numeric.Float) * MA.masked # axis0: replicas, axis1: factor B levels
            for groupIdx,takeInd in enumerate(grpInd):
                m2[:groupLens[groupIdx]*ma3d.shape[1], groupIdx] = MA.ravel(ma3d[eIdx].take(takeInd, 1))
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
            ma3d = ma3d.take(tInd2keep, 1)
        # for each gene...
        for eIdx in range(ma3d.shape[0]):
            # faster check for empty cells for that gene -> remove time indices with empty cells
            ma2d = ma3d[eIdx]
            cellCount = MA.zeros((ma2d.shape[0], groupLens.shape[0]), Numeric.Int)
            for g,(i0,i1) in enumerate(zip(ax2Ind[:-1], ax2Ind[1:])):
                cellCount[:,g] = MA.count(ma2d[:,i0:i1], 1)
            ma2dTakeInd = Numeric.logical_not(Numeric.add.reduce(Numeric.equal(cellCount,0),1)) # 1 where to take, 0 where not to take
            if Numeric.add.reduce(ma2dTakeInd) != ma2dTakeInd.shape[0]:
                print "Warning: removing time indices %s for gene %i" % (str(Numeric.compress(ma2dTakeInd == 0, Numeric.arange(ma2dTakeInd.shape[0]))), eIdx)
                ma2d = MA.compress(ma2dTakeInd, ma2d, 0)
            an = fAnova(ma2d, groupLens, addInteraction, allowReductA=True, allowReductB=True)
            ps[:,eIdx] = an.ps
            callback()
        return ps


    def updateAnovaTypeBox(self):
        """enables/disables: - anova type selection box;
                             - example selection box;
                             - selectors A, B and I
                             - Commit button
        """
        if self.dataStructure and self.numExamples > 0:
            # enable anova type box and commit button
            self.boxAnovaType.setEnabled(1)
            self.btnCommit.setEnabled(1)
            # select appropriate anova type
            if len(self.dataStructure) == 1 and self.numVariables == 1:
                # single-sample t-test (factor A)
                self.anovaType = 0
            elif len(self.dataStructure) == 1 and self.numVariables > 1:
                # single-factor (A) ANOVA
                self.anovaType = 1
            elif len(self.dataStructure) > 1 and self.numVariables == 1:
                # single-factor (B) ANOVA
                self.anovaType = 2
            elif len(self.dataStructure) > 1 and self.numVariables > 1:
                # two-factor ANOVA
                self.anovaType = int(self._interaction) + 3
            # enable/disable appropriate anova type radio buttons
            if self.anovaType <= 2:
                for i in range(5):
                    self.boxAnovaType.buttons[i].setEnabled(self.anovaType == i)
            else:
                self.boxAnovaType.buttons[0].setEnabled(0)
                for i in range(1,5):
                    self.boxAnovaType.buttons[i].setEnabled(1)
            # enable/disable compareTo lineEdit
            self.boxCompareTo.setEnabled(self.anovaType == 0)
        else:
            # disable anova type box and commit button
            self.boxAnovaType.setDisabled(1)
            self.btnCommit.setDisabled(1)


    def updateSelectorBox(self):
        """enables / disables individual selectors
        """
        if self.dataStructure and self.numExamples > 0:
            # enable example selection box
            self.boxSelection.setEnabled(1)
            # enable/disable selectors A, B and I
            self.boxSelectorA.setEnabled(self.anovaType != 2)
            self.boxSelectorB.setEnabled(self.anovaType >= 2)
            self.boxSelectorI.setEnabled(self.anovaType == 4)
        else:
            # disable example selection box
            self.boxSelection.setDisabled(1)


    def updateSelectorInfos(self, selectorIdx=None):
        """updates the number of examples that match individual selectors;
        if selectorIdx is given, updates only the corresponding info.
        """
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
            if ps != None and alpha != None and self.anovaType in [[0,1,3,4],[2,3,4],[4]][si]:
                numSelected = Numeric.add.reduce(Numeric.less_equal(self.ps[si], alpha))
                self.lblNumGenes[si].setText('  (%d example%s)' % (numSelected, ['', 's'][numSelected!=1]))
            else:
                self.lblNumGenes[si].setText('  (no examples)')


    def senddata(self):
        """computes selectionList, partitions the examples and updates infoc;
        sends out selectionList and selected/other dataStructure or None;
        """
##        if self.dataStructure and self.ps:
        if self.dataStructure and self.ps.shape[1]:
            # set selectionList
            alphas = [self.alphaA, self.alphaB, self.alphaI]
            selectors = [self.selectorA, self.selectorB, self.selectorI]
            selectionList = Numeric.ones((self.numExamples,))
            for si in range(3):
                try:
                    if selectors[si] and self.anovaType in [[0,1,3,4],[2,3,4],[4]][si]:
                        selectionList = Numeric.logical_and(selectionList, Numeric.less_equal(self.ps[si], float(alphas[si])))
                except:
                    pass
            self.infoc.setText('Sending out data...')
            
            if self.sendProbabilities:
                # create example table with probabilities
                print self.ps
                print Numeric.transpose(self.ps).shape
                etProb = orange.ExampleTable(orange.Domain([orange.FloatVariable("Factor A p-val"),orange.FloatVariable("Factor B p-val"),orange.FloatVariable("Interaction p-val")]), Numeric.transpose(self.ps))
                # in etProb, convert p-val to meta attribute
                domProb = orange.Domain([])
                domProb.addmetas(dict(zip([orange.newmetaid(),orange.newmetaid(),orange.newmetaid()], etProb.domain.variables)))
                etProb = orange.ExampleTable(domProb, etProb)
            else:
                # create new etProb without attributes/metas and of length equal to etProb
                etProb = orange.ExampleTable(orange.Domain([]), Numeric.zeros((selectionList.shape[0],0)))

            # partition dataStructure and send out data
            selectionList = selectionList.tolist()
            self.send("Example Selection", (self.selectorName, selectionList))
            dataStructS = []
            dataStructN = []
            self.progressBarInit()

            if self.sendNotSelectedData:
                pbStep = 50./len(self.dataStructure)
            else:
                pbStep = 100./len(self.dataStructure)

            for (dsName, etList) in self.dataStructure:
                etListS = [et.select(selectionList) for et in etList]
                for i in range(len(etList)):
                    # append probabilities (if etProb not empty)
                    etListS[i] = orange.ExampleTable([etListS[i], etProb.select(selectionList)])
                    # add name
                    etListS[i].name = etList[i].name
                dataStructS.append((dsName, etListS))
                self.progressBarAdvance(pbStep)
            self.send("Selected Structured Data", dataStructS)

            if self.sendNotSelectedData:
                for (dsName, etList) in self.dataStructure:
                    etListN = [et.select(selectionList, negate=1) for et in etList]
                    for i in range(len(etList)):
                        # append probabilities (if etProb not empty)
                        etListN[i] = orange.ExampleTable([etListN[i], etProb.select(selectionList, negate=1)])
                        # add name
                        etListN[i].name = etList[i].name
                    dataStructN.append((dsName, etListN))
                    self.progressBarAdvance(pbStep)
                self.send("Other Structured Data", dataStructN)
            else:
                self.send("Other Structured Data", None)

            self.progressBarFinished()
            # report the number of selected examples
            numExamples = Numeric.add.reduce(Numeric.greater(selectionList, 0))
            self.infoc.setText('Total of %d example%s match criteria.' % (numExamples, ['', 's'][numExamples!=1]))
        else:
            self.send("Example Selection", None)
            self.send("Selected Structured Data", None)
            self.send("Other Structured Data", None)
            

    def updateSelectorName(self):
        """update selector name shown in selector edit box
        """
        if self.dataStructure:
            if self.anovaType == 0:
                s = '1 smpl. t-test, compared to %s' % self.compareToValue
            else:
                s = 'ANOVA'
            s += " (%s)" % reduce(lambda a,b: a + ", " + b[0], self.dataStructure, "")[2:]
            if self.selectorA and self.anovaType in [0,1,3,4]:
                s += ", pA<%s" % self.alphaA
            if self.selectorB and self.anovaType in [2,3,4]:
                s += ", pB<%s" % self.alphaB
            if self.selectorI and self.anovaType == 4:
                s += ", pI<%s" % self.alphaI
            self.selectorName = s.strip()
        else:
            self.selectorName = ""


    #==========================================================================
    # Event handlers
    #==========================================================================

    def onCompareToChange(self):
        """handles changes of ANOVA type:
            - resets self.ps;
            - updates infoc
            - calls updateSelectorInfos()
        runs ANOVA and sends out new data;
        """
        if self.anovaType == 0:
            self.ps = None
            if self.autoUpdateSelName:
                self.updateSelectorName()
            if self.commitOnChange:
                self.runANOVA()
                self.senddata()
            elif self.dataStructure and self.numExamples > 0:
                self.infoc.setText('Press Commit button to start ANOVA computation.')
            self.updateSelectorInfos()


    def onAnovaType(self):
        """handles changes of ANOVA type:
            - resets self.ps;
            - calls updateSelectorBox()
            - updates infoc
            - calls updateSelectorInfos()
        runs ANOVA and sends out new data;
        """
        #print "self.anovaType", self.anovaType
        self._interaction = self.anovaType == 4
        self.ps = None
        self.updateSelectorBox()
        if self.autoUpdateSelName:
            self.updateSelectorName()
        if self.commitOnChange:
            self.runANOVA()
            self.senddata()
        elif self.dataStructure and self.numExamples > 0:
            self.infoc.setText('Press Commit button to start ANOVA computation.')
        self.updateSelectorInfos()


    def onInteraction(self):
        """handles clicks on interaction checkbox:
            - resets self.ps;
            - enables/disables selector I,
            - updates infoc
        runs ANOVA and sends out new data
        """
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


    def onSelectionChange(self):
        """handles changes in example selector checkboxes;
        sends out new data;
        """
        if self.autoUpdateSelName:
            self.updateSelectorName()
        if self.commitOnChange:
            self.senddata()


    def onAlphaChange(self, selectorIdx):
        """handles changes in example selector alphas;
        prints number of selected examples for individual selectors and sends out new data;
        """
        if self.autoUpdateSelName:
            self.updateSelectorName()
        if self.commitOnChange:
            self.senddata()
        self.updateSelectorInfos(selectorIdx)


    def onAutoUpdateSelNameChange(self):
        """handles clicks on auto update selector name checkbox
        """
        self.leSelectorName.setReadOnly(self.autoUpdateSelName)

            
    def onSendNotSelectedChange(self):
        """handles clicks on sendNotSelectedData checkbox
        """
        if self.commitOnChange:
            self.senddata()

    def onSendProbabilitiesChange(self):            
        """handles clicks on show p-values checkbox
        """
        if self.commitOnChange:
            self.senddata()
            
    def onCommit(self, commit=True):
        """handles Commit clicks; runs ANOVA (if not already computed) and sends out data;
        """
        if commit:
            if self.dataStructure:
                if self.ps == None:
                    self.runANOVA()
                self.senddata()
            self.updateSelectorInfos()


if __name__=="__main__":
    from . import OWDataFiles, OWDataFilesSelector
    from Orange.orng import orngSignalManager
    from Orange.OrangeWidgets.Data import OWDataTable

    signalManager = orngSignalManager.SignalManager(0)
    a=QApplication(sys.argv)
    an=OWANOVA(signalManager = signalManager)

    an.show()
    
    df = OWDataFiles.OWDataFiles(signalManager = signalManager)
    df.loadData("/home/marko/anova/smallchipdata")
    df.show()

    signalManager.addWidget(an)
    signalManager.addWidget(df)
    signalManager.setFreeze(1)
    signalManager.addLink(df, an, 'Structured Data', 'Structured Data', 1)

    # data files selector, data table
    dfs = OWDataFilesSelector.OWDataFilesSelector(signalManager = signalManager)
    signalManager.addWidget(dfs)
    dfs.show()
    signalManager.addLink(an, dfs, 'Selected Structured Data', 'Structured Data', 1)
    dt = OWDataTable.OWDataTable(signalManager = signalManager)
    signalManager.addWidget(dt)
    signalManager.addLink(dfs, dt, 'Examples', 'Examples', 1)
    signalManager.setFreeze(0)
    dt.show()

    a.exec_()
    an.saveSettings()
