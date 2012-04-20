## Automatically adapted for numpy.oldnumeric Oct 04, 2007 by 

"""
<name>Hypotheses test</name>
<description>Single Sample / Related Samples t-test, One/Two Way Analysis of Variance, Local-Pooled-Error Two Sample Z-test.</description>
<icon>icons/ChipANOVA.png</icon>
<priority>1071</priority>
<contact>Peter Juvan (peter.juvan@fri.uni-lj.si)</contact>
<prototype>1</prototype>
"""

from __future__ import absolute_import

import qt, qwt

import numpy.oldnumeric as Numeric, numpy.oldnumeric.ma as MA

import scipy.stats

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .. import Anova
from .OWDataFiles import DataFiles, ExampleSelection

class OWHypTest(OWWidget):
    settingsList  = ["anovaType", "popMean", "_interaction", "selectorA", "selectorB", "selectorI", "alphaA", "alphaB", "alphaI", "autoUpdateSelName", "sendNotSelectedData", "sendProbabilities", "commitOnChange"]

    StNames = ["Single sample t-test", "Related samples t-test", "[not impl.: LPE (local-pooled-error two sample z-test)]", "Single-factor ANOVA (A, variables)", "Single-factor ANOVA (B, data sets)", "Two-factor ANOVA", "Two-factor ANOVA with interaction effect"]
    StSST  = 0   # single sample t-test
    StRST  = 1   # related seamples t-test (either by A or B)
    StLPE  = 2   # Local-pool-error two sample z-test
    St1A   = 3   # One-way ANOVA by variables (A)
    St1B   = 4   # One-way ANOVA by datasets (B)
    St2AB  = 5   # Two-way ANOVA
    St2ABI = 6   # Two-way ANOVA with interaction

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'ANOVA')
        # input / output data: [("name1", [orange.ExampleTable1a,...]), ("name2", [orange.ExampleTable2a,...])]
        self.inputs = [("Structured Data", DataFiles, self.onDataInput)]
        self.outputs = [("Example Selection", ExampleSelection, Default), ("Selected Structured Data", DataFiles, Default), ("Other Structured Data", DataFiles)]

        # data, p-values, selected examples
        self.dataStructure = None                       # input data
        self.numExamples = 0
        self.numVariables = 0
        self.ps = Numeric.ones((3,0), Numeric.Float)    # p-values: 2D Numeric.array of shape (3, numExamples)
        self.selectorName = ""                          # for Example Selection output: (self.selectorName, [0,1,0,...])

        # Settings
        self.anovaType = OWHypTest.StSST
        self.popMean = 0         # single sample t-test, value to compare to
        self.useFactors = [0,0,0]       # [use factor A, use factor B, use interaction]
        self._interaction = 0           # to store last setting: 0: no interaction, 1: test for interaction effect (set this value manually !!!)
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
        ca=QFrame(self.controlArea)
        gl=QGridLayout(ca,4,1,5)
        
        # info
        box = QVGroupBox("Info", ca)
        gl.addWidget(box,0,0)
        self.infoa = QLabel('No data on input.', box)
        self.infob = QLabel('', box)
        self.infoc = QLabel('', box)

        # ANOVA type
        self.boxAnovaType = QVButtonGroup("Statistics", ca)
        gl.addWidget(self.boxAnovaType,1,0)
        self.boxAnovaType.setDisabled(1)

        self.boxAnovaType.setRadioButtonExclusive(1)
        self.boxAnovaType.buttons = []
##        for i,lbl in enumerate(OWHypTest.StNames):
##            w = QRadioButton(lbl, self.boxAnovaType)
##            w.setOn(self.anovaType == i)
##            self.boxAnovaType.buttons.append(w)
##            if i == OWHypTest.StSST:
##                self.boxPopMean = QHBox(self.boxAnovaType)
##                QLabel("             population mean  ", self.boxPopMean)
##                OWGUI.lineEdit(self.boxPopMean, self, "popMean", callback=self.onPopMeanChange)
        for i,lbl in enumerate(OWHypTest.StNames):
            w = QRadioButton(lbl, self.boxAnovaType)
            w.setOn(self.anovaType == i)
            self.boxAnovaType.buttons.append(w)
            if i == OWHypTest.StSST:
                self.boxPopMean = QHBox(self.boxAnovaType)
                QLabel("       population mean  ", self.boxPopMean)
                OWGUI.lineEdit(self.boxPopMean, self, "popMean", callback=self.onPopMeanChange)
        OWGUI.connectControl(self.boxAnovaType, self, "anovaType", self.onAnovaType, "clicked(int)", OWGUI.CallFront_radioButtons(self.boxAnovaType))
        
        # selection of examples
        self.boxSelection = QVGroupBox("Example Selection", ca)
        gl.addWidget(self.boxSelection,2,0)
        self.lblNumGenes = []   # list of labels
        # selector A
        self.boxSelectorA = QVBox(self.boxSelection)
        self.cbSelectorA = OWGUI.checkBox(self.boxSelectorA, self, "selectorA", "Factor A (variables)", callback=self.onSelectionChange,
                                          tooltip='H0: The mean does not depend on factor A (represented by variables).')
        frmA = QFrame(self.boxSelectorA)
        glA = QGridLayout(frmA,1,3,5)
        leA = OWGUI.lineEdit(frmA, self, "alphaA", orientation="horizontal", controlWidth=None, callback=lambda x=0: self.onAlphaChange(x))
        glA.addWidget(leA,0,1) # Qt.AlignRight
        glA.addWidget(QLabel("     p <= ", frmA), 0,0)
        self.lblNumGenes.append(QLabel('', frmA))
        glA.addWidget(self.lblNumGenes[-1],0,2) # Qt.AlignRight | 0x22

        # selector B
        self.boxSelectorB = QVBox(self.boxSelection)
        self.cbSelectorB = OWGUI.checkBox(self.boxSelectorB, self, "selectorB", "Factor B (data sets)", callback=self.onSelectionChange,
                                          tooltip='H0: The mean does not depend on factor B (represented by data sets).')
        frmB = QFrame(self.boxSelectorB)
        glB = QGridLayout(frmB,1,3,5)
        leB = OWGUI.lineEdit(frmB, self, "alphaB", orientation="horizontal", callback=lambda x=1: self.onAlphaChange(x))
        glB.addWidget(leB,0,1)
        glB.addWidget(QLabel("     p <= ", frmB), 0,0)
        self.lblNumGenes.append(QLabel('', frmB))
        glB.addWidget(self.lblNumGenes[-1],0,2)
        
        # selector I
        self.boxSelectorI = QVBox(self.boxSelection)
        self.cbSelectorI = OWGUI.checkBox(self.boxSelectorI, self, "selectorI", "Interaction (variables * data sets)", callback=self.onSelectionChange,
                                          tooltip='H0: There is no interaction between factor A and factor B.')
        frmI = QFrame(self.boxSelectorI)
        glI = QGridLayout(frmI,1,3,5)
        leI = OWGUI.lineEdit(frmI, self, "alphaI", orientation="horizontal", callback=lambda x=2: self.onAlphaChange(x))
        ## slider could be used to replace lineEdit (but not sensitive enough)
        ##        self.alphaIf = 0.05
        ##        leI = OWGUI.qwtHSlider(self.boxSelectorI, self, "alphaIf", box="", label="      p < ", labelWidth=None, minValue=0.0001, maxValue=1.0, step=0.1, precision=3, callback=lambda x=2: self.onAlphaChange(x), logarithmic=1, ticks=0, maxWidth=None)
        glI.addWidget(leI,0,1)
        glI.addWidget(QLabel("     p <= ", frmI), 0,0)
        self.lblNumGenes.append(QLabel('', frmI))
        glI.addWidget(self.lblNumGenes[-1],0,2)
        
        # output
        box = QVGroupBox("Output", ca)
        gl.addWidget(box,3,0)
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
        self.ps = Numeric.ones((3,0), Numeric.Float)
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
            # run ANOVA
            self.infoc.setText('ANOVA computation started...')
            self.progressBarInit()
            pbStep = 100./self.numExamples
            self.ps = Numeric.ones((3, self.numExamples), Numeric.Float)

            if self.anovaType == OWHypTest.St2AB or self.anovaType == OWHypTest.St2ABI:
                ps = self.anova2(ma3d, groupLens, addInteraction=self.anovaType==OWHypTest.St2ABI, repMeasuresOnA=False, callback=lambda: self.progressBarAdvance(pbStep))
                for rIdx in range(ps.shape[0]):
                    self.ps[rIdx] = ps[rIdx]

            elif self.anovaType == OWHypTest.St1B:
                self.ps[1] = self.anova1B(ma3d, groupLens, repMeasures=False, callback=lambda: self.progressBarAdvance(pbStep))

            elif self.anovaType == OWHypTest.St1A:
                self.ps[0] = self.anova1A(ma3d, repMeasures=False, callback=lambda: self.progressBarAdvance(pbStep))

            elif self.anovaType == OWHypTest.StSST:
                try:
                    popMeanVal = float(self.popMean)
                except ValueError:
                    print "Warning: cannot convert %s to float, using 0" % str(self.popMean)
                    self.popMean = 0
                    popMeanVal = 0
                self.ps[0] = self.ttest_ssmpl(ma3d, popMeanVal, callback=lambda: self.progressBarAdvance(pbStep))

            elif self.anovaType == OWHypTest.StLPE:
               raise Exception, "NOT IMPLEMENTED"
               if self.numVariables == 2:
                  self.ps[0] = self.lpeA(ma3d, callback=lambda: self.progressBarAdvance(pbStep))
               elif self.numVariables == 1:
                  self.ps[1] = self.lpeB(ma3d, groupLens, callback=lambda: self.progressBarAdvance(pbStep))
               else:
                  raise RuntimeError, "%s: expected 2 variables and 1 group, or 1 variable and 2 groups, got %s variables and %s groups" % (OWHypTest.StNames[self.anovaType], self.numVariables, len(groupLens))

            elif self.anovaType == OWHypTest.StRST:
               if self.numVariables == 2 and len(groupLens) == 1:
                  self.ps[0] = self.ttest_rsmplA(ma3d, callback=lambda: self.progressBarAdvance(pbStep))
               elif self.numVariables == 1 and len(groupLens) == 2 and groupLens[0] == groupLens[1]:
                  self.ps[1] = self.ttest_rsmplB(ma3d, groupLens, callback=lambda: self.progressBarAdvance(pbStep))
               else:
                  raise RuntimeError, "%s: expected 2 variables and 1 group, or 1 variable and 2 groups of equal length, got %s variables and %s groups of length %s" % (OWHypTest.StNames[self.anovaType], self.numVariables, len(groupLens), str(groupLens))
                  
            self.progressBarFinished()


    def lpeA(self, ma3d, callback):
        """conducts local-pooled-error test;
        reference: Jain et.al. (2003) Bioinformatics.
        """
        raise Exception, "NOT IMPLEMENTED"

    def lpeB(self, ma3d, callback):
        """conducts local-pooled-error test;
        reference: Jain et.al. (2003) Bioinformatics.
        """
        raise Exception, "NOT IMPLEMENTED"


    def ttest_ssmpl(self, ma3d, popMeanVal, callback):
        """conducts single-sample t-test on individual examples wrt factor A (variables, ma3d axis 1);
        returns Numeric array of p-values in shape (1, numExamples).
        """
        ps = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
        for eIdx in range(ma3d.shape[0]):
            data = Numeric.asarray(MA.transpose(ma3d[eIdx]).compressed())
            if len(data) >= 2:
                try:
                    ps[eIdx] = scipy.stats.ttest_1samp(data, popMeanVal)[1]
                except:
                    print "Warning: zero variance, check the example %i:" % eIdx, data
                    ps[eIdx] = 1.0
            else:
##                print "Warning: removing example %i:\n%s\n%s\n" % (eIdx, str(data))
                print "Warning: removing example %i:" % eIdx, str(data)
                ps[eIdx] = 1.0
            callback()
        return ps


    def ttest_rsmplA(self, ma3d, callback):
        """conducts related samples t-test on individual examples wrt factor A (variables, ma3d axis 1);
        returns Numeric array of p-values in shape (1, numExamples).
        """
        ps = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
        for eIdx in range(ma3d.shape[0]):
            a = ma3d[eIdx][0]
            b = ma3d[eIdx][1]
            cond = Numeric.logical_not(Numeric.logical_or(MA.getmaskarray(a), MA.getmaskarray(b)))
            a = Numeric.asarray(MA.compress(cond, a))
            b = Numeric.asarray(MA.compress(cond, b))
            if len(a) >= 2:
                try:
                    ps[eIdx] = scipy.stats.ttest_rel(a,b)[1]
                except Exception, inst:
                    print "Warning: %s" % str(inst)
                    print "Example %i:\n%s\n%s\n" % (eIdx, str(a), str(b))
                    ps[eIdx] = 1.0
            else:
                print "Warning: removing example %i:\n%s\n%s\n" % (eIdx, str(a), str(b))
                ps[eIdx] = 1.0
            callback()
        return ps


    def ttest_rsmplB(self, ma3d, groupLens, callback):
        """conducts related samples t-test on individual examples wrt factor B (datasets, ma3d axis 2);
        ma3d axis 2 also contains replicas according to groupLens;
        returns Numeric array of p-values in shape (1, numExamples).
        """
        ps = -1*Numeric.ones((ma3d.shape[0],), Numeric.Float)
        for eIdx in range(ma3d.shape[0]):
            a = ma3d[eIdx][0][:groupLens[0]]
            b = ma3d[eIdx][0][groupLens[0]:]
            cond = Numeric.logical_not(Numeric.logical_or(MA.getmaskarray(a), MA.getmaskarray(b)))
            a = Numeric.asarray(MA.compress(cond, a))
            b = Numeric.asarray(MA.compress(cond, b))
            if len(a) >= 2:
                try:
                    ps[eIdx] = scipy.stats.ttest_rel(a,b)[1]
                except Exception, inst:
                    print "Warning: %s" % str(inst)
                    print "Example %i:\n%s\n%s\n" % (eIdx, str(a), str(b))
                    ps[eIdx] = 1.0
            else:
                print "Warning: removing example %i:\n%s\n%s\n" % (eIdx, str(a), str(b))
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
        """Sets self.anovaType according to the size of the dataset;
           enables/disables: - anova type selection box;
                             - individual anova type radio buttons;
                             - compareTo box;
                             - commit button
        """
        if self.dataStructure and self.numExamples > 0:
            # enable anova type box and commit button
            self.boxAnovaType.setEnabled(1)
            self.btnCommit.setEnabled(1)
            # disable all radio buttions
            for i in range(len(OWHypTest.StNames)):
                self.boxAnovaType.buttons[i].setEnabled(0)
            # select appropriate anova type and enable corresponding radio buttons
            if len(self.dataStructure) == 1 and self.numVariables == 1:
                # single-sample t-test (factor A)
                self.anovaType = OWHypTest.StSST
                self.boxAnovaType.buttons[self.anovaType].setEnabled(1)
                self.boxAnovaType.buttons[self.anovaType].setChecked(1)
            elif len(self.dataStructure) == 1 and self.numVariables == 2:
                # RST or LPE or ANOVA1 by variables (A)
                self.boxAnovaType.buttons[OWHypTest.StRST].setEnabled(1)
                self.boxAnovaType.buttons[OWHypTest.StLPE].setEnabled(1)
                self.boxAnovaType.buttons[OWHypTest.St1A].setEnabled(1)
                if not self.boxAnovaType.buttons[OWHypTest.StLPE].isChecked() and not self.boxAnovaType.buttons[OWHypTest.St1A].isChecked():
                    self.anovaType = OWHypTest.StRST
                    self.boxAnovaType.buttons[self.anovaType].setChecked(1)
            elif len(self.dataStructure) == 2 and self.numVariables == 1:
                # RST or LPE or ANOVA1 by datasets (B)
                self.boxAnovaType.buttons[OWHypTest.StRST].setEnabled(1)
                self.boxAnovaType.buttons[OWHypTest.StLPE].setEnabled(1)
                self.boxAnovaType.buttons[OWHypTest.St1B].setEnabled(1)
                if not self.boxAnovaType.buttons[OWHypTest.StLPE].isChecked() and not self.boxAnovaType.buttons[OWHypTest.St1B].isChecked():
                    self.anovaType = OWHypTest.StRST
                    self.boxAnovaType.buttons[self.anovaType].setChecked(1)
            elif len(self.dataStructure) == 1 and self.numVariables > 1:
                # single-factor (A) ANOVA
                self.anovaType = OWHypTest.St1A
                self.boxAnovaType.buttons[self.anovaType].setEnabled(1)
                self.boxAnovaType.buttons[self.anovaType].setChecked(1)
            elif len(self.dataStructure) > 1 and self.numVariables == 1:
                # single-factor (B) ANOVA
                self.anovaType = OWHypTest.St1B
                self.boxAnovaType.buttons[self.anovaType].setEnabled(1)
                self.boxAnovaType.buttons[self.anovaType].setChecked(1)
            elif len(self.dataStructure) > 1 and self.numVariables > 1:
                # two-factor ANOVA
                self.boxAnovaType.buttons[OWHypTest.St1A].setEnabled(1)
                self.boxAnovaType.buttons[OWHypTest.St1B].setEnabled(1)
                self.boxAnovaType.buttons[OWHypTest.St2AB].setEnabled(1)
                self.boxAnovaType.buttons[OWHypTest.St2ABI].setEnabled(1)
                if not self.boxAnovaType.buttons[OWHypTest.St1A].isChecked() and not self.boxAnovaType.buttons[OWHypTest.St1B].isChecked() and not self.boxAnovaType.buttons[OWHypTest.St2AB].isChecked() and not self.boxAnovaType.buttons[OWHypTest.St2ABI].isChecked():
                    if self._interaction:
                        self.anovaType = OWHypTest.St2ABI
                    else:
                        self.anovaType = OWHypTest.St2AB
                    self.boxAnovaType.buttons[self.anovaType].setChecked(1)
            # enable/disable compareTo lineEdit
            self.boxPopMean.setEnabled(self.anovaType == OWHypTest.StSST)
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
            self.boxSelectorA.setEnabled(self.anovaType == OWHypTest.StSST or
                                         self.anovaType == OWHypTest.St1A or
                                         (self.anovaType == OWHypTest.StRST and self.numVariables == 2) or
                                         (self.anovaType == OWHypTest.StLPE and self.numVariables == 2) or
                                         self.anovaType == OWHypTest.St2AB or
                                         self.anovaType == OWHypTest.St2ABI)
            self.boxSelectorB.setEnabled(self.anovaType == OWHypTest.St1B or
                                         (self.anovaType == OWHypTest.StRST and self.numVariables == 1) or
                                         (self.anovaType == OWHypTest.StLPE and self.numVariables == 1) or
                                         self.anovaType == OWHypTest.St2AB or
                                         self.anovaType == OWHypTest.St2ABI)
            self.boxSelectorI.setEnabled(self.anovaType == OWHypTest.St2ABI)
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
        boxSelectors = [self.boxSelectorA, self.boxSelectorB, self.boxSelectorI]
        for si in selectorInd:
            try:
                alpha = float(alphas[si])
                ps = self.ps[si]
            except ValueError:
                alpha = None
                ps = None
##            if ps != None and alpha != None and self.anovaType in [[0,1,3,4],[2,3,4],[4]][si]:
            if ps != None and alpha != None and boxSelectors[si].isEnabled():                
                numSelected = Numeric.add.reduce(Numeric.less(self.ps[si], alpha))
                self.lblNumGenes[si].setText('  (%d example%s)' % (numSelected, ['', 's'][int(numSelected!=1)]))
            else:
                self.lblNumGenes[si].setText('  (no examples)')


    def senddata(self):
        """computes selectionList, partitions the examples and updates infoc;
        sends out selectionList and selected/other dataStructure or None;
        """
        if self.dataStructure and self.ps.shape[1]:
            # set selectionList
            alphas = [self.alphaA, self.alphaB, self.alphaI]
            selectors = [self.selectorA, self.selectorB, self.selectorI]
            selectionList = Numeric.ones((self.numExamples,))
            boxSelectors = [self.boxSelectorA, self.boxSelectorB, self.boxSelectorI]
            for si in range(3):
                try:
##                    if selectors[si] and self.anovaType in [[0,1,3,4],[2,3,4],[4]][si]:
                    if selectors[si] and boxSelectors[si].isEnabled():
                        selectionList = Numeric.logical_and(selectionList, Numeric.less(self.ps[si], float(alphas[si])))
                except ValueError:
                    print "Warning: cannot convert %s to float" % str(alphas[si])
                    pass
            self.infoc.setText('Sending out data...')
            
            if self.sendProbabilities:
                # create example table with probabilities
##                print self.ps
##                print Numeric.transpose(self.ps).shape
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
            self.infoc.setText('Total of %d example%s match criteria.' % (numExamples, ['', 's'][int(numExamples!=1)]))
        else:
            self.send("Example Selection", None)
            self.send("Selected Structured Data", None)
            self.send("Other Structured Data", None)
            

##    def updateSelectorName(self):
##        """update selector name shown in selector edit box
##        """
##        if self.dataStructure:
##            if self.anovaType == 0:
##                s = '1 smpl. t-test, compared to %s' % self.popMean
##            else:
##                s = 'ANOVA'
##            s += " (%s)" % reduce(lambda a,b: a + ", " + b[0], self.dataStructure, "")[2:]
##            if self.selectorA and self.anovaType in [0,1,3,4]:
##                s += ", pA<%s" % self.alphaA
##            if self.selectorB and self.anovaType in [2,3,4]:
##                s += ", pB<%s" % self.alphaB
##            if self.selectorI and self.anovaType == 4:
##                s += ", pI<%s" % self.alphaI
##            self.selectorName = s.strip()
##        else:
##            self.selectorName = ""
    def updateSelectorName(self):
        """update selector name shown in selector edit box
        """
        if self.dataStructure:
            s = OWHypTest.StNames[self.anovaType]
            if self.anovaType == OWHypTest.StSST:
               s += " (pop. mean: %s)" % self.popMean
            s += " (%s)" % reduce(lambda a,b: a + ", " + b[0], self.dataStructure, "")[2:]
            if self.selectorA and self.boxSelectorA.isEnabled():
                s += ", pA<%s" % self.alphaA
            if self.selectorB and self.boxSelectorB.isEnabled():
                s += ", pB<%s" % self.alphaB
            if self.selectorI and self.boxSelectorI.isEnabled():
                s += ", pI<%s" % self.alphaI
            self.selectorName = s.strip()
        else:
            self.selectorName = ""


    #==========================================================================
    # Event handlers
    #==========================================================================

    def onPopMeanChange(self):
        """handles changes of ANOVA type:
            - resets self.ps;
            - updates infoc
            - calls updateSelectorInfos()
        runs ANOVA and sends out new data;
        """
        if self.anovaType == OWHypTest.StSST:
            self.ps = Numeric.ones((3,0), Numeric.Float)
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
        # store info whether we have tested for interaction effect
        self._interaction = self.anovaType == OWHypTest.St2ABI
        self.ps = Numeric.ones((3,0), Numeric.Float)
        self.updateSelectorBox()
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
                if self.ps.shape[1] == 0:
                    self.runANOVA()
                self.senddata()
            self.updateSelectorInfos()


if __name__=="__main__":
    from . import OWDataFiles
    from Orange.orng import orngSignalManager
    signalManager = orngSignalManager.SignalManager(0)
    a=QApplication(sys.argv)
    ow=OWHypTest(signalManager = signalManager)
    a.setMainWidget(ow)
    ow.show()
    ds = OWDataFiles.OWDataFiles(signalManager = signalManager)
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\potato.sub100")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\potato.sub1000")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS_10_yakApufA")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS_10_yakApufA_time0")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS_10_yakApufA_time0_swappedAB")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\DictyChipData_BR_ACS")

##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\_one-sample t-test")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\_factor A")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\_factor B")
##    ds.loadData(r"C:\Documents and Settings\peterjuv\My Documents\Orange\ANOVA\_factors A and B")

    signalManager.addWidget(ow)
    signalManager.addWidget(ds)
    signalManager.setFreeze(1)
    signalManager.addLink(ds, ow, 'Structured Data', 'Structured Data', 1)
    signalManager.setFreeze(0)
    a.exec_loop()
    ow.saveSettings()
