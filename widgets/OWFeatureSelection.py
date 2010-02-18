"""
<name>Gene Selection</name>
<description>Gene differential expression scoring and selection.</description>
<priority>1010</priority>
<icon>icons/FeatureSelection.png</icon>
"""

from __future__ import with_statement
import orange

from obiExpression import *

from OWGraph import *
from OWGraphTools import PolygonCurve
from OWWidget import *
from OWHist import OWInteractiveHist
from OWToolbars import ZoomSelectToolbar
from obiGEO import transpose
from collections import defaultdict
import numpy as np
import numpy.ma as ma

import OWGUI
        
class ExpressionSignificance_TTest_PValue(ExpressionSignificance_TTest):
    def __call__(self, *args, **kwargs):
        return [(key, pval) for key, (t, pval) in ExpressionSignificance_TTest.__call__(self, *args, **kwargs)]
    
class ExpressionSignificance_TTest_T(ExpressionSignificance_TTest):
    def __call__(self, *args, **kwargs):
        return [(key, t) for key, (t, pval) in ExpressionSignificance_TTest.__call__(self, *args, **kwargs)]
    
class ExpressionSignificance_ANOVA_PValue(ExpressionSignificance_ANOVA):
    def __call__(self, *args, **kwargs):
        return [(key, pval) for key, (t, pval) in ExpressionSignificance_ANOVA.__call__(self, *args, **kwargs)]
    
class ExpressionSignificance_ANOVA_F(ExpressionSignificance_ANOVA):
    def __call__(self, *args, **kwargs):
        return [(key, f) for key, (f, pval) in ExpressionSignificance_ANOVA.__call__(self, *args, **kwargs)]
    
class ExpressionSignificance_Log2FoldChange(ExpressionSignificance_FoldChange):
    def __call__(self, *args, **kwargs):
        return [(key, math.log(fold, 2.0) if fold > 1e-300 and fold < 1e300 else 0.0) for key, fold in ExpressionSignificance_FoldChange.__call__(self, *args, **kwargs)]
    
class ExpressionSignigicance_MannWhitneyu_U(ExpressionSignificance_MannWhitneyu):
    def __call__(self, *args, **kwargs):
        return [(key, u) for key, (u, p_val) in ExpressionSignificance_MannWhitneyu.__call__(self, *args, **kwargs)]    
    
class ScoreHist(OWInteractiveHist):
    def __init__(self, master, parent=None, type="hiTail"):
        OWInteractiveHist.__init__(self, parent, type=type)
        self.master = master
        self.setAxisTitle(QwtPlot.xBottom, "Score")
        self.setAxisTitle(QwtPlot.yLeft, "Frequency")
        
    def setBoundary(self, low, hi):
        OWInteractiveHist.setBoundary(self, low, hi)
        self.master.UpdateSelectedInfoLabel(low, hi)
        self.master.CommitIf()
            
        
class OWFeatureSelection(OWWidget):
    settingsList = ["methodIndex", "dataLabelIndex", "computeNullDistribution", "permutationsCount", "selectPValue", "autoCommit"]
##    contextHandlers={"":DomainContextHandler("",[])}
    def __init__(self, parent=None, signalManager=None, name="Gene selection"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, showSaveGraph=True)
        self.inputs = [("Examples", ExampleTable, self.SetData)]
        self.outputs = [("Example table with selected genes", ExampleTable), ("Example table with remaining genes", ExampleTable), ("Selected genes", ExampleTable)]

        self.methodIndex = 0
        self.dataLabelIndex = 0
        self.computeNullDistribution = False
        self.permutationsCount = 10
        self.autoCommit = False
        self.selectNBest = 20
        self.selectPValue = 0.01
        self.dataChangedFlag = False
        self.addScoresToOutput = True

        self.oneTailTestHi = oneTailTestHi = lambda array, low, hi: array >= hi
        self.oneTailTestLow = oneTailTestLow = lambda array, low, hi: array <= low
        self.twoTailTest = twoTailTest = lambda array, low, hi: (array >= hi) | (array <= low)
        self.middleTest = middleTest = lambda array, low, hi: (array <= hi) | (array >= low)
        
        self.histType = {oneTailTestHi:"hiTail", oneTailTestLow:"lowTail", twoTailTest:"twoTail", middleTest:"middle"}

        self.scoreMethods = [("fold change", ExpressionSignificance_FoldChange, twoTailTest),
                             ("log2 fold change", ExpressionSignificance_Log2FoldChange, twoTailTest),
                             ("t-test", ExpressionSignificance_TTest_T, twoTailTest),
                             ("t-test p-value", ExpressionSignificance_TTest_PValue, oneTailTestLow),
                             ("anova", ExpressionSignificance_ANOVA_F, oneTailTestHi),
                             ("anova p-value", ExpressionSignificance_ANOVA_PValue, oneTailTestLow),
                             ("signal to noise ratio", ExpressionSignificance_SignalToNoise, twoTailTest),
                             ("info gain", ExpressionSignificance_Info, oneTailTestHi),
                             ("chi-square", ExpressionSignificance_ChiSquare, oneTailTestHi),
                             ("mann-whitney", ExpressionSignigicance_MannWhitneyu_U, oneTailTestLow)]

        boxHistogram = OWGUI.widgetBox(self.mainArea)
        self.histogram = ScoreHist(self, boxHistogram)
        boxHistogram.layout().addWidget(self.histogram)
        self.histogram.show()
        
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.dataInfoLabel = OWGUI.widgetLabel(box, "\n\n")
        self.selectedInfoLabel = OWGUI.widgetLabel(box, "")

        #self.testRadioBox = OWGUI.radioButtonsInBox(self.controlArea, self, "methodIndex", [sm[0] for sm in self.scoreMethods], box="Scoring Method", callback=self.Update, addSpace=True)
        # OWGUI.comboBoxWithCaption(box, self, "ptype", items=nth(self.permutationTypes, 0), \
        #             tooltip="Permutation type.", label="Permutate")
        box1 = OWGUI.widgetBox(self.controlArea, "Scoring Method")
        self.testRadioBox = OWGUI.comboBox(box1, self, "methodIndex", items=[sm[0] for sm in self.scoreMethods], callback=self.Update)
        self.dataLabelComboBox = OWGUI.comboBox(box1, self, "dataLabelIndex", "Attribute labels", callback=self.Update, tooltip="Use attribute labels for score computation")
##        self.useClassCheck = OWGUI.checkBox(box1, self, "useClass", "Use class information", callback=self.Update, tooltip="Use class information for score computation", disables=[self.da)
    
        ZoomSelectToolbar(self, self.controlArea, self.histogram, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])
        OWGUI.separator(self.controlArea)
        
        box = OWGUI.widgetBox(self.controlArea, "Selection", addSpace=True)

        callback = self.SetBoundary
        self.upperBoundarySpin = OWGUI.doubleSpin(box, self, "histogram.upperBoundary", min=-1e6, max=1e6, step= 1e-6, label="Upper threshold:", callback=callback, callbackOnReturn=True)
        self.lowerBoundarySpin = OWGUI.doubleSpin(box, self, "histogram.lowerBoundary", min=-1e6, max=1e6, step= 1e-6, label="Lower threshold:", callback=callback, callbackOnReturn=True)
        check = OWGUI.checkBox(box, self, "computeNullDistribution", "Compute null distribution", callback=self.Update)

        check.disables.append(OWGUI.spin(box, self, "permutationsCount", min=1, max=10, label="Permutations:", callback=self.Update, callbackOnReturn=True))

        box1 = OWGUI.widgetBox(box, orientation='horizontal')
        check.disables.append(OWGUI.doubleSpin(box1, self, "selectPValue" , min=2e-7, max=1.0, step=1e-7, label="P-value:"))
        check.disables.append(OWGUI.button(box1, self, "Select", callback=self.SelectPBest))
        check.makeConsistent()

        box1 = OWGUI.widgetBox(box, orientation='horizontal')
        OWGUI.spin(box1, self, "selectNBest", 0, 10000, step=1, label="Best Ranked:")
        OWGUI.button(box1, self, "Select", callback=self.SelectNBest)

        box = OWGUI.widgetBox(self.controlArea, "Output")
        b = OWGUI.button(box, self, "&Commit", callback=self.Commit)
        cb = OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.setStopper(self, b, cb, "dataChangedFlag", self.Commit)
        OWGUI.checkBox(box, self, "addScoresToOutput", "Add gene scores to output", callback=self.CommitIf) 
        
        OWGUI.rubber(self.controlArea)

        self.connect(self.graphButton, SIGNAL("clicked()"), self.histogram.saveToFile)
        
        self.loadSettings()

        self.data = None
        self.discData = None
        self.scoreCache = {}
        self.nullDistCache = {}
        self.cuts = {}
        self.discretizer = orange.EquiNDiscretization(numberOfIntervals=5)
        self.transposedData = False
        self.nullDistribution = []
        self.targets = []
        self.scores = {}

        self.resize(800, 600)
        
    def SetData(self, data):
        self.error(0)
        self.warning(0)
        self.scoreCache = {}
        self.nullDistCache = {}
        self.discData = None
        self.data = data
        self.transposedData = None
        disabled = []
        if self.data:
            self.dataLabels = reduce(lambda dict, tags: [dict[key].add(value) for key, value in tags.items()] and False or dict,
                                       [attr.attributes for attr in self.data.domain.attributes],
                                       defaultdict(set))
            self.dataLabels = [key for key, value in self.dataLabels.items() if len(value) > 1]
        else:
            self.dataLabels = []
        self.dataLabelComboBox.clear()
        if self.data and data.domain.classVar:
            self.dataLabels = ["(None)"] + self.dataLabels
        self.dataLabelComboBox.addItems(self.dataLabels)
        
        self.dataLabelComboBox.setDisabled(len(self.dataLabels) <= 1)
        self.dataLabelIndex = max(min(self.dataLabelIndex, len(self.dataLabels) - 1), 0)
        
        self.Update()
        if not self.data:
            self.UpdateDataInfoLabel()
            self.send("Example table with selected genes", None)
            self.send("Example table with remaining genes", None)
            self.send("Selected genes", None)
    
    def ComputeScores(self, data, scoreFunc, useAttributeLabels, target=None, advance=lambda :None):
        scoreFunc = scoreFunc(data, useAttributeLabels)
        advance()
        score = scoreFunc(target=target)
        score = [(key, val) for key, val in score if val is not ma.masked]
        return score
    
    def ComputeNullDistribution(self, data, scoreFunc, useAttributes, target=None, permCount=10, advance=lambda: None):
        scoreFunc = scoreFunc(data, useAttributes)
        dist = scoreFunc.null_distribution(permCount, target, advance=advance)
        return [score for run in dist for k, score in run]
        
    def Update(self):
        if not self.data:
            self.histogram.removeDrawingCurves()
            self.histogram.clear()
            return
        targetLabel = self.dataLabels[self.dataLabelIndex]
        if targetLabel == "(None)":
            target = [self.data.domain.classVar(0), self.data.domain.classVar(1)]
            self.targets = targets = list(self.data.domain.classVar.values)
            self.genesInColumns = False
        else:
            self.targets = targets = list(set([attr.attributes.get(targetLabel) for attr in self.data.domain.attributes]) - set([None]))
            target = [(targetLabel, targets[0]), (targetLabel, targets[1])]
            self.genesInColumns = True
        if self.methodIndex in [4, 5]: ## ANOVA
            target = targets if targetLabel == "(None)" else [(targetLabel, t) for t in targets]
        scoreFunc = self.scoreMethods[self.methodIndex][1] 
        pb = OWGUI.ProgressBar(self, 4 + self.permutationsCount if self.computeNullDistribution else 3)
        self.scores = dict(self.ComputeScores(self.data, scoreFunc, self.genesInColumns, target, advance=pb.advance))
        pb.advance()
        if self.computeNullDistribution:
            self.nullDistribution = self.ComputeNullDistribution(self.data, scoreFunc, self.genesInColumns, target, self.permutationsCount, advance=pb.advance)
        pb.advance()
        self.histogram.type = self.histType[self.scoreMethods[self.methodIndex][2]]
        if self.scores:
            self.histogram.setValues(self.scores.values())
#        self.histogram.setBoundary(*self.cuts.get(None, (self.histogram.minx if self.histogram.type in ["lowTail", "twoTail"] else self.histogram.maxx,
#                                                                     self.histogram.maxx if self.histogram.type in ["hiTail", "twoTail"] else self.histogram.minx)))
            self.histogram.setBoundary(self.histogram.minx if self.histogram.type in ["lowTail", "twoTail"] else self.histogram.maxx,
                                       self.histogram.maxx if self.histogram.type in ["hiTail", "twoTail"] else self.histogram.minx)
            if self.computeNullDistribution:
                nullY, nullX = numpy.histogram(self.nullDistribution, bins=self.histogram.xData)
                self.histogram.nullCurve = self.histogram.addCurve("nullCurve", Qt.black, Qt.black, 6, symbol = QwtSymbol.NoSymbol, style = QwtPlotCurve.Steps, xData = nullX, yData = nullY/self.permutationsCount)
                
                minx = min(min(nullX), self.histogram.minx)
                maxx = max(max(nullX), self.histogram.maxx)
                miny = min(min(nullY/self.permutationsCount), self.histogram.miny)
                maxy = max(max(nullY/self.permutationsCount), self.histogram.maxy)

                self.histogram.setAxisScale(QwtPlot.xBottom, minx - (0.05 * (maxx - minx)), maxx + (0.05 * (maxx - minx)))
                self.histogram.setAxisScale(QwtPlot.yLeft, miny - (0.05 * (maxy - miny)), maxy + (0.05 * (maxy - miny)))
            state = dict(hiTail=(False, True), lowTail=(True, False), twoTail=(True, True))
            for spin, visible in zip((self.upperBoundarySpin, self.lowerBoundarySpin), state[self.histogram.type]):
                spin.setVisible(visible)
            
##            if self.methodIndex in [2, 3, 5, 6]:
            if self.methodIndex in [0, 2, 6]:
                if self.methodIndex == 0: ## fold change is centered on 1.0
                    x1, y1 = (self.histogram.minx + 1) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx + 1) / 2 , self.histogram.maxy
                else:
                    x1, y1 = (self.histogram.minx) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx) / 2 , self.histogram.maxy
                self.histogram.addMarker(self.targets[0], x1, y1)
                self.histogram.addMarker(self.targets[1], x2, y2)
            self.warning(0)
        else:
            self.warning(0, "No scores obtained.")
        self.histogram.replot()
        pb.advance()
        pb.finish()
        self.UpdateDataInfoLabel()
            
            
    def UpdateDataInfoLabel(self):
        if self.data:
            samples, genes = len(self.data), len(self.data.domain.attributes)
            if self.genesInColumns:
                samples, genes = genes, samples
            text = "%i samples, %i genes\n" % (samples, genes)
            text += "Sample label: %s" % self.targets[0]
        else:
            text = "No data on input\n"
        self.dataInfoLabel.setText(text)

    def UpdateSelectedInfoLabel(self, cutOffLower=0, cutOffUpper=0):
        self.cuts[self.methodIndex] = (cutOffLower, cutOffUpper)
        if self.data:
            scores = np.array(self.scores.values())
            test = self.scoreMethods[self.methodIndex][2]
            self.selectedInfoLabel.setText("%i selected genes" % len(np.nonzero(test(scores, cutOffLower, cutOffUpper))[0]))
        else:
            self.selectedInfoLabel.setText("0 selected genes")

    def SelectNBest(self):
        scores = self.scores.items()
        scores.sort(lambda a,b:cmp(a[1], b[1]))
        if not scores:
            return
        if self.scoreMethods[self.methodIndex][2]==self.oneTailTestHi:
            scores = scores[-max(self.selectNBest, 1):]
            self.histogram.setBoundary(scores[0][1], scores[0][1])
        elif self.scoreMethods[self.methodIndex][2]==self.oneTailTestLow:
            scores = scores[:max(self.selectNBest,1)]
            self.histogram.setBoundary(scores[-1][1], scores[-1][1])
        else:
            scoresHi = scores[-max(min(self.selectNBest, len(scores)/2), 1):]
            scoresLo = scores[:max(min(self.selectNBest, len(scores)/2), 1)]
            scores = [(abs(score), 1) for attr, score in scoresHi] + [(abs(score), -1) for attr, score in scoresLo]
            if self.scoreMethods[self.methodIndex][0]=="fold change": ## comparing fold change on a logaritmic scale
                scores =  [(abs(math.log(max(min(score, 1e300), 1e-300), 2.0)), sign) for score, sign in scores]
            scores.sort()
            scores = scores[-max(self.selectNBest, 1):]
            countHi = len([score for score, sign in scores if sign==1])
            countLo = len([score for score, sign in scores if sign==-1])
            cutHi = scoresHi[-countHi][1] if countHi else scoresHi[-1][1] + 1e-7
            cutLo = scoresLo[countLo-1][1] if countLo else scoresLo[0][1] - 1e-7
            self.histogram.setBoundary(cutLo, cutHi)

    def SetBoundary(self):
        if self.data != None and self.scores.items():
            if self.scoreMethods[self.methodIndex][2]==self.oneTailTestHi:
                self.histogram.setBoundary(self.histogram.lowerBoundary, self.histogram.lowerBoundary)
            elif self.scoreMethods[self.methodIndex][2]==self.oneTailTestLow:
                self.histogram.setBoundary(self.histogram.upperBoundary, self.histogram.upperBoundary)
            else:
                self.histogram.setBoundary(self.histogram.lowerBoundary, self.histogram.upperBoundary)

    def SelectPBest(self):
        if not self.nullDistribution:
            return
        nullDist = sorted(self.nullDistribution)
        test = self.scoreMethods[self.methodIndex][2]
        count = int(len(nullDist)*self.selectPValue)
        if test == self.oneTailTestHi:
            cut = nullDist[-count] if count else nullDist[-1] # + 1e-7
            self.histogram.setBoundary(cut, cut)
        elif test == self.oneTailTestLow:
            cut = nullDist[count - 1] if count else nullDist[0] # - 1e-7
            self.histogram.setBoundary(cut, cut)
        elif count:
            scoresHi = nullDist[-count:]
            scoresLo = nullDist[:count]
            scores = [(abs(score), 1) for score in scoresHi] + [(abs(score), -1) for score in scoresLo]
            if self.scoreMethods[self.methodIndex][0] == "fold change": ## fold change is on a logaritmic scale
                scores =  [(abs(math.log(max(min(score, 1e300), 1e-300), 2.0)), sign) for score, sign in scores]
            scores = sorted(scores)[-count:]
            countHi = len([score for score, sign in scores if sign==1])
            countLo = len([score for score, sign in scores if sign==-1])
            cutHi = scoresHi[-countHi] if countHi else scoresHi[-1] + 1e-7
            cutLo = scoresLo[countLo-1] if countLo else scoresLo[0] - 1e-7
            self.histogram.setBoundary(cutLo, cutHi)
        else:
            self.histogram.setBoundary(nullDist[0] - 1e-7, nullDist[-1] + 1e-7)
            

    def CommitIf(self):
        if self.autoCommit:
            self.Commit()
        else:
            self.dataChangedFlag = True
        
    def Commit(self):
        if not self.data:
            return
        test = self.scoreMethods[self.methodIndex][2]
        
        cutOffUpper = self.histogram.upperBoundary
        cutOffLower = self.histogram.lowerBoundary
        
        scores = np.array(self.scores.items())
        scores[:, 1] = test(np.array(scores[:, 1], dtype=float), cutOffLower, cutOffUpper)
        selected = set([key for key, test in scores if test])
        remaining = set([key for key, test in scores if not test])
        if self.data and self.genesInColumns:
            selected = sorted(selected)
            newdata = orange.ExampleTable(orange.Domain(self.data.domain), [self.data[int(i)] for i in selected])
            if self.addScoresToOutput:
                scoreAttr = orange.FloatVariable(self.scoreMethods[self.methodIndex][0])
                mid = orange.newmetaid()
                newdata.domain.addmeta(mid, scoreAttr)
                for ex, key in zip(newdata, selected):
                    ex[mid] = self.scores[key]
            self.send("Example table with selected genes", newdata)
            remaining = sorted(remaining)
            newdata = orange.ExampleTable(orange.Domain(self.data.domain), [self.data[int(i)] for i in remaining]) #ex in enumerate(self.data) if i in remaining])
            if self.addScoresToOutput:
                newdata.domain.addmeta(mid, scoreAttr)
                for ex, key in zip(newdata, remaining):
                    ex[mid] = self.scores[key]
            self.send("Example table with remaining genes", newdata)
            
        elif self.data and not self.genesInColumns:
            
            selectedAttrs = [attr for attr in self.data.domain.attributes if attr in selected or attr.varType == orange.VarTypes.String]
            newdomain = orange.Domain(selectedAttrs, self.data.domain.classVar)
            if self.addScoresToOutput:
                for attr in newdomain.attributes:
                    attr.attributes[self.scoreMethods[self.methodIndex][0]] = str(self.scores[attr])
            newdomain.addmetas(self.data.domain.getmetas())
            newdata = orange.ExampleTable(newdomain, self.data)
            self.send("Example table with selected genes", newdata if selectedAttrs else None)
            
            remainingAttrs = [attr for attr in self.data.domain.attributes if attr in remaining]
            newdomain = orange.Domain(remainingAttrs, self.data.domain.classVar)
            if self.addScoresToOutput:
                for attr in newdomain.attributes:
                    attr.attributes[self.scoreMethods[self.methodIndex][0]] = str(self.scores[attr])
            newdomain.addmetas(self.data.domain.getmetas())
            newdata = orange.ExampleTable(newdomain, self.data)
            self.send("Example table with remaining genes", newdata if remainingAttrs else None)
            
            domain = orange.Domain([orange.StringVariable("label"), orange.FloatVariable(self.scoreMethods[self.methodIndex][0])], False)
            self.send("Selected genes", orange.ExampleTable([orange.Example(domain, [attr.name, self.scores.get(attr, 0)]) for attr in selectedAttrs]) if selectedAttrs else None)
            
        else:
            self.send("Example table with selected genes", None)
            self.send("Example table with remaining genes", None)
            self.send("Selected genes", None)
        self.dataChangedFlag = False

if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
#    data = orange.ExampleTable("E:\\out1.tab")
#    data = orange.ExampleTable("/home/marko/t2.tab")
    data = orange.ExampleTable("../../../doc/datasets/brown-selected")
    w = OWFeatureSelection()
    w.show()
    w.SetData(data)
    app.exec_()
    w.saveSettings()
