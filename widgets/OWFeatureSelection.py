"""
<name>Feature Selection</name>
<description>Feature selection</description>
<priority>230</priority>
<icon>icons/FeatureSelection.png</icon>
"""
import orange

from obiExpression import *

from OWGraph import *
from OWGraphTools import PolygonCurve
from OWWidget import *
from OWHist import OWInteractiveHist
from OWToolbars import ZoomSelectToolbar
from obiGEO import transpose

import OWGUI

class ScoreHist(OWInteractiveHist):
    def __init__(self, master, parent=None, type="hiTail"):
        OWInteractiveHist.__init__(self, parent, type=type)
        self.master = master
        self.setAxisTitle(QwtPlot.xBottom, "Score")
        self.setAxisTitle(QwtPlot.yLeft, "Frequency")
        
    def setBoundary(self, low, hi):
        OWInteractiveHist.setBoundary(self, low, hi)
        self.master.UpdateSelectedInfoLabel(low, hi)
        if self.master.autoCommit:
            self.master.Commit()
        
class OWFeatureSelection(OWWidget):
    settingsList=["methodIndex", "computeNullDistribution", "permutationsCount", "selectPValue", "autoCommit"]
##    contextHandlers={"":DomainContextHandler("",[])}
    def __init__(self, parent=None, signalManager=None, name="Feature selection"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, showSaveGraph=True)
        self.inputs = [("Examples", ExampleTable, self.SetData)]
        self.outputs = [("Examples with selected attributes", ExampleTable), ("Examples with remaining attributes", ExampleTable), ("Selected attributes", ExampleTable)]

        self.methodIndex = 0
        self.computeNullDistribution = False
        self.permutationsCount = 10
        self.autoCommit = False
        self.selectNBest = 20
        self.selectPValue = 0.01

        self.oneTailTestHi = oneTailTestHi = lambda attr, low, hi:self.scores.get(attr,0)>=hi
        self.oneTailTestLow = oneTailTestLow = lambda attr, low, hi:self.scores.get(attr,0)<=low
        self.twoTailTest = twoTailTest = lambda attr, low, hi:self.scores.get(attr,0)>=hi or self.scores.get(attr,0)<=low
        self.middleTest = middleTest = lambda attr, low, hi:self.scores.get(attr,0)<=hi and self.scores.get(attr,0)>=low
        self.histType = {oneTailTestHi:"hiTail", oneTailTestLow:"lowTail", twoTailTest:"twoTail", middleTest:"middle"}
##        self.scoreMethods = [("chi-square", orange.MeasureAttribute_chiSquare, oneTailTestHi),
##                             ("info gain", orange.MeasureAttribute_info, oneTailTestHi),
##                             ("signal to noise ratio", lambda attr, data: MA_signalToNoise()(attr, data), twoTailTest),
##                             ("t-test",lambda attr, data: MA_t_test()(attr, data), twoTailTest),
##                             ("t-test p-value",lambda attr, data: MA_t_test(prob=True)(attr, data), oneTailTestLow),
##                             ("fold change", lambda attr, data: MA_fold_change()(attr, data), twoTailTest),
##                             ("log2 fold change", lambda attr, data: math.log(max(min(MA_fold_change()(attr, data), 1e300), 1e-300), 2.0), twoTailTest),
##                             ("anova", lambda attr, data: MA_anova()(attr, data), oneTailTestHi),
##                             ("anova p-value", lambda attr, data: MA_anova(prob=True)(attr, data), oneTailTestLow)]
        self.scoreMethods = [("fold change", lambda attr, data: MA_fold_change()(attr, data), twoTailTest),
                             ("log2 fold change", lambda attr, data: math.log(max(min(MA_fold_change()(attr, data), 1e300), 1e-300), 2.0), twoTailTest),
                             ("t-test",lambda attr, data: MA_t_test()(attr, data), twoTailTest),
                             ("t-test p-value",lambda attr, data: MA_t_test(prob=True)(attr, data), oneTailTestLow),
                             ("anova", lambda attr, data: MA_anova()(attr, data), oneTailTestHi),
                             ("anova p-value", lambda attr, data: MA_anova(prob=True)(attr, data), oneTailTestLow),
                             ("signal to noise ratio", lambda attr, data: MA_signalToNoise()(attr, data), twoTailTest),
                             ("info gain", orange.MeasureAttribute_info, oneTailTestHi),
                             ("chi-square", orange.MeasureAttribute_chiSquare, oneTailTestHi)]

        boxHistogram = OWGUI.widgetBox(self.mainArea)
        self.histogram = ScoreHist(self, boxHistogram)
        boxHistogram.layout().addWidget(self.histogram)
##        self.histogram = ScoreHist(self)
##        self.mainArea.layout().addWidget(self.histogram)
##        self.histogram.setSizePolicy(QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum))
        self.histogram.show()
        
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.dataInfoLabel = OWGUI.widgetLabel(box, "\n\n")
        self.selectedInfoLabel = OWGUI.widgetLabel(box, "")
        self.testRadioBox = OWGUI.radioButtonsInBox(self.controlArea, self, "methodIndex", [sm[0] for sm in self.scoreMethods], box="Score Method", callback=self.Update, addSpace=True)
        ZoomSelectToolbar(self, self.controlArea, self.histogram, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])
        OWGUI.separator(self.controlArea)
        
        box = OWGUI.widgetBox(self.controlArea, "Selection", addSpace=True)
        box2 = OWGUI.widgetBox(box, "Threshold")
        callback = lambda: self.histogram.setBoundary(self.histogram.lowerBoundary, self.histogram.upperBoundary)
        self.upperBoundarySpin = OWGUI.doubleSpin(box2, self, "histogram.upperBoundary", min=-1e6, max=1e6, step= 1e-6, label="Upper:", callback=callback, callbackOnReturn=True)
        self.lowerBoundarySpin = OWGUI.doubleSpin(box2, self, "histogram.lowerBoundary", min=-1e6, max=1e6, step= 1e-6, label="Lower:", callback=callback, callbackOnReturn=True)
        check = OWGUI.checkBox(box, self, "computeNullDistribution", "Compute null distribution", callback=self.Update)
        check.disables.append(OWGUI.spin(box, self, "permutationsCount", min=1, max=10, label="Repetitions:", callback=self.Update, callbackOnReturn=True))
        check.disables.append(OWGUI.button(box, self, "Select w.r.t null distribution", callback=self.SelectPBest))
        check.disables.append(OWGUI.doubleSpin(box, self, "selectPValue" , min=2e-7, max=1.0, step=1e-7, label="p-value:"))
        check.makeConsistent()
        OWGUI.button(box, self, "Select n best features", callback=self.SelectNBest)
        OWGUI.spin(box, self, "selectNBest", 0, 10000, step=1, label="n:")
        box = OWGUI.widgetBox(self.controlArea, "Commit")
        OWGUI.button(box, self, "&Commit", callback=self.Commit)
        OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
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

        self.resize(700, 600)        
        
    def SetData(self, data):
        self.error(0)
        self.warning(0)
        self.scoreCache = {}
        self.nullDistCache = {}
        self.discData = None
        self.data = data
        self.transposedData = False
        disabled = []
        if self.data and not data.domain.classVar:
            try:
                self.data = data = transpose(data)
                self.transposedData = True
            except Exception, ex:
                self.error(0, "Class var missing! " + str(ex))
                self.data = None
        if self.data and len(self.data.domain.classVar.values) == 2:
            disabled = [4, 5]
        elif self.data and len(self.data.domain.classVar.values) > 2:
           disabled = [0, 1, 2, 3, 6]
        for i, button in enumerate(self.testRadioBox.buttons):
            button.setDisabled(i in disabled)
        self.UpdateDataInfoLabel()
        self.Update()
        self.Commit()

    def ComputeAttributeScore(self, data, scoreFunc, useCache=True, progressCallback=None):
        if scoreFunc in self.scoreCache and useCache:
            return self.scoreCache[scoreFunc]
        attributes = data.domain.attributes
##        if self.methodIndex in [2,3,4,5,6] and len(data.domain.classVar.values)>2:
        if self.methodIndex in [0, 1, 2, 3, 6] and len(data.domain.classVar.values) > 2:
            self.warning(0, self.scoreMethods[self.methodIndex][0]+" works only on two class data (using only first two class values for computation)")
        else:
            self.warning(0)
##        self.progressBarInit()
        if scoreFunc==orange.MeasureAttribute_info or scoreFunc==orange.MeasureAttribute_chiSquare:
            if self.discData and useCache:
                data = self.discData
                newAttrs = data.domain.attributes
            else:
                newAttrs = [self.discretizer(attr, data) for attr in attributes]
##                data = data.select(newAttrs + [data.domain.classVar])
                newDomain = orange.Domain(newAttrs + [data.domain.classVar])
                table = orange.ExampleTable(newDomain)
                for i, ex in enumerate(data):
                    table.append(orange.Example(newDomain, ex))
                    self.progressBarSet(100.0*i/len(data))
                if useCache:
                    self.discData = data = table
        else:
            newAttrs = attributes
        scores = {}
        milestones = set(range(0, len(attributes), max(len(attributes)/100, 1)))
        for i, (attr, newAttr) in enumerate(zip(attributes, newAttrs)):
            scores[attr] = scoreFunc(newAttr, data)
            if progressCallback and i in milestones:
                progressCallback(100.0*i/len(attributes))
##        self.progressBarFinished()
        if useCache:
            self.scoreCache[scoreFunc] = scores
        return scores

    def ComputeNullDistribution(self, data, scoreFunc, progressCallback=None):
        if (scoreFunc, self.permutationsCount) in self.nullDistCache:
            return self.nullDistCache[scoreFunc, self.permutationsCount]

        originalClasses = [ex.getclass() for ex in data]
        scores = []
        import random
        for i in range(self.permutationsCount):
            permClasses = list(originalClasses)
            random.shuffle(permClasses)
            for ex, class_ in zip(data, permClasses):
                ex.setclass(class_)
            _progressCallback = lambda val: progressCallback(100.0*i/self.permutationsCount + float(val)/self.permutationsCount) if \
                               progressCallback else None
            scores.extend(self.ComputeAttributeScore(data, scoreFunc, useCache=False, progressCallback=_progressCallback).values())

        for ex, class_ in zip(data, originalClasses):
            ex.setclass(class_)

        self.nullDistCache[scoreFunc, self.permutationsCount] = scores
        return scores
        
    def Update(self):
        if self.data and self.data.domain.classVar:
            self.progressBarInit()
            self.scores = self.ComputeAttributeScore(self.data, self.scoreMethods[self.methodIndex][1], progressCallback=self.progressBarSet)
            if self.computeNullDistribution:
                self.nullDistribution = self.ComputeNullDistribution(self.data, self.scoreMethods[self.methodIndex][1], progressCallback=self.progressBarSet)
            self.progressBarFinished()
            self.histogram.type = self.histType[self.scoreMethods[self.methodIndex][2]]
            self.histogram.setValues(self.scores.values())
            self.histogram.setBoundary(*self.cuts.get(self.methodIndex, (0, 0)))
            if self.computeNullDistribution:
                nullY, nullX = numpy.histogram(self.nullDistribution, bins=100)
                self.histogram.nullCurve = self.histogram.addCurve("nullCurve", Qt.black, Qt.black, 6, symbol = QwtSymbol.NoSymbol, style = QwtPlotCurve.Steps, xData = nullX, yData = nullY/self.permutationsCount)
            state = dict(hiTail=(False, True), lowTail=(True, False), twoTail=(True, True))
            for spin, visible in zip((self.upperBoundarySpin, self.lowerBoundarySpin), state[self.histogram.type]):
                spin.setVisible(visible)
            
##            if self.methodIndex in [2, 3, 5, 6]:
            if self.methodIndex in [0, 2, 4, 6]:
                classValues = self.data.domain.classVar.values
                if self.methodIndex == 0: ## fold change is centered on 1.0
                    x1, y1 = (self.histogram.minx + 1) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx + 1) / 2 , self.histogram.maxy
                else:
                    x1, y1 = (self.histogram.minx) / 2 , self.histogram.maxy
                    x2, y2 = (self.histogram.maxx) / 2 , self.histogram.maxy
                self.histogram.addMarker(classValues[1], x1, y1)
                self.histogram.addMarker(classValues[0], x2, y2)
            self.histogram.replot()
            
        else:
            self.histogram.clear()
            
    def UpdateDataInfoLabel(self):
        if self.data:
            text = "%i data instances with\n%i attributes\n" % (len(self.data), len(self.data.domain.attributes))
            if self.data.domain.classVar:
                text = text+"Class var: %s" % self.data.domain.classVar.name
            else:
                text = text+"Class var missing"
        else:
            text = "No data on input\n\n"
        self.dataInfoLabel.setText(text)

    def UpdateSelectedInfoLabel(self, cutOffLower=0, cutOffUpper=0):
        self.cuts[self.methodIndex] = (cutOffLower, cutOffUpper)
        if self.data:
            test = self.scoreMethods[self.methodIndex][2]
            self.selectedInfoLabel.setText("%i selected attributes" %len([attr for attr in self.data.domain.attributes if test(attr, cutOffLower, cutOffUpper)])) #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]))
        else:
            self.selectedInfoLabel.setText("0 selected attributes")

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
            if self.scoreMethods[self.methodIndex][0]=="fold change": ## fold change is on a logaritmic scale
                scores =  [(abs(math.log(max(min(score, 1e300), 1e-300), 2.0)), sign) for score, sign in scores]
            scores.sort()
            scores = scores[-max(self.selectNBest, 1):]
            countHi = len([score for score, sign in scores if sign==1])
            countLo = len([score for score, sign in scores if sign==-1])
            cutHi = scoresHi[-countHi][1] if countHi else scoresHi[-1][1] + 1e-7
            cutLo = scoresLo[countLo-1][1] if countLo else scoresLo[0][1] - 1e-7
            self.histogram.setBoundary(cutLo, cutHi)

    def SelectPBest(self):
        if not self.nullDistribution:
            return
        nullDist = sorted(self.nullDistribution)
        test = self.scoreMethods[self.methodIndex][2]
        count = int(len(nullDist)*self.selectPValue)
        if test == self.oneTailTestHi:
            cut = nullDist[-count] if count else nullDist[-1] + 1e-7
            self.histogram.setBoundary(cut, cut)
            print cut
        elif test == self.oneTailTestLow:
            cut = nullDist[count - 1] if count else nullDist[0] - 1e-7
            self.histogram.setBoundary(cut, cut)
            print cut
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
        
    def Commit(self):
        if self.data and self.data.domain.classVar:
            cutOffUpper = self.histogram.upperBoundary
            cutOffLower = self.histogram.lowerBoundary
            test = self.scoreMethods[self.methodIndex][2]
            
            selectedAttrs = [attr for attr in self.data.domain.attributes if  test(attr, cutOffLower, cutOffUpper)] #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]
            data = self.data.select(selectedAttrs+[self.data.domain.classVar])
            if self.transposedData and selectedAttrs:
                data = transpose(data)
            self.send("Examples with selected attributes", data)
            
            remainingAttrs = [attr for attr in self.data.domain.attributes if  not test(attr, cutOffLower, cutOffUpper)] #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]
            data = self.data.select(remainingAttrs+[self.data.domain.classVar])
            if self.transposedData and remainingAttrs:
                data = transpose(data)
            self.send("Examples with remaining attributes", data)
            
            domain = orange.Domain([orange.StringVariable("label"), orange.FloatVariable(self.scoreMethods[self.methodIndex][0])], False)
            self.send("Selected attributes", orange.ExampleTable([orange.Example(domain, [attr.name, self.scores.get(attr, 0)]) for attr in selectedAttrs]) if selectedAttrs else None)
            
        else:
            self.send("Examples with selected attributes", None)
            self.send("Examples with remaining attributes", None)
            self.send("Selected attributes", None)

if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    data = orange.ExampleTable("E://leukemia.tab")
    w = OWFeatureSelection()
    w.show()
    w.SetData(data)
    app.exec_()
    w.saveSettings()
