"""
<name>FeatureSelection</name>
<description>Feature selection</description>
<priority>230</priority>
"""
import orange

from obiExpression import *

from OWGraph import *
from OWGraphTools import PolygonCurve
from OWWidget import *
from OWHist import OWInteractiveHist
from OWToolbars import ZoomSelectToolbar

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
    settingsList=["methodIndex", "autoCommit"]
##    contextHandlers={"":DomainContextHandler("",[])}
    def __init__(self, parent=None, signalManager=None, name="Feature selection"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, showSaveGraph=True)
        self.inputs = [("Examples", ExampleTable, self.SetData)]
        self.outputs = [("Examples with selected attributes", ExampleTable), ("Examples with remaining attributes", ExampleTable)]

        self.methodIndex = 0
        self.autoCommit = False
        self.selectNBest = 20
##        self.infoLabel = "No data on input"        

        self.oneTailTestHi = oneTailTestHi = lambda attr, low, hi:self.scores.get(attr,0)>=hi
        self.oneTailTestLow = oneTailTestLow = lambda attr, low, hi:self.scores.get(attr,0)<=low
        self.twoTailTest = twoTailTest = lambda attr, low, hi:self.scores.get(attr,0)>=hi or self.scores.get(attr,0)<=low
        self.middleTest = middleTest = lambda attr, low, hi:self.scores.get(attr,0)<=hi and self.scores.get(attr,0)>=low
        self.histType = {oneTailTestHi:"hiTail", oneTailTestLow:"lowTail", twoTailTest:"twoTail", middleTest:"middle"}
        self.scoreMethods = [("chi-square", orange.MeasureAttribute_chiSquare, oneTailTestHi),
                             ("info gain", orange.MeasureAttribute_info, oneTailTestHi),
                             ("signal to noise ratio", lambda attr, data:MA_signalToNoise()(attr, data), twoTailTest),
                             ("t-test",lambda attr, data:MA_t_test()(attr, data), twoTailTest),
                             ("t-test p-value",lambda attr, data:1.0 - MA_t_test(prob=True)(attr, data), oneTailTestLow),
                             ("fold change", lambda attr, data:MA_fold_change()(attr, data), twoTailTest),
                             ("log2 fold change", lambda attr, data:math.log(max(min(MA_fold_change()(attr, data), 1e300), 1e-300)), twoTailTest),
                             ("anova", lambda attr, data:MA_anova()(attr, data), oneTailTestHi),
                             ("anova p-value", lambda attr, data:1.0-MA_anova(prob=True)(attr, data), oneTailTestLow)]

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
        OWGUI.radioButtonsInBox(self.controlArea, self, "methodIndex", [sm[0] for sm in self.scoreMethods], box="Score Method", callback=self.Update, addSpace=True)
        ZoomSelectToolbar(self, self.controlArea, self.histogram, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])
        OWGUI.separator(self.controlArea)
        box = OWGUI.widgetBox(self.controlArea, "Selection", addSpace=True)
        OWGUI.button(box, self, "Select n best features", callback=self.SelectNBest)
        OWGUI.spin(box, self, "selectNBest", 0, 10000, step=1, label="n:")
        OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.button(box, self, "&Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)

        self.connect(self.graphButton, SIGNAL("clicked()"), self.histogram.saveToFile)
        
        self.loadSettings()

        self.data = None
        self.discData = None
        self.scoreCache = {}
        self.cuts = {}
        self.discretizer = orange.EquiNDiscretization(numberOfIntervals=5)
        
    def SetData(self, data):
        self.error(0)
        self.warning(0)
        self.scoreCache = {}
        self.discData = None
        self.data = data
        if self.data and not data.domain.classVar:
            self.error(0, "Class var missing!")
        self.UpdateDataInfoLabel()
        self.Update()
        self.Commit()

    def ComputeAttributeScore(self, data, scoreFunc):
        if scoreFunc in self.scoreCache:
            return self.scoreCache[scoreFunc]
        attributes = data.domain.attributes
        if self.methodIndex in [2,3,4,5,6] and len(data.domain.classVar.values)>2:
            self.warning(0, self.scoreMethods[self.methodIndex][0]+" works only on two class data (using only first two class values for computation)")
        else:
            self.warning(0)
        self.progressBarInit()
        if scoreFunc==orange.MeasureAttribute_info or scoreFunc==orange.MeasureAttribute_chiSquare:
            if self.discData:
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
                self.discData = data = table
        else:
            newAttrs = attributes
        scores = {}
        milestones = set(range(0, len(attributes), max(len(attributes)/100, 1)))
        for i, (attr, newAttr) in enumerate(zip(attributes, newAttrs)):
            scores[attr] = scoreFunc(newAttr, data)
            if i in milestones:
                self.progressBarSet(100.0*i/len(attributes))
        self.progressBarFinished()
        self.scoreCache[scoreFunc] = scores
        return scores
        
    def Update(self):
        if self.data and self.data.domain.classVar:
            self.scores = self.ComputeAttributeScore(self.data, self.scoreMethods[self.methodIndex][1])
            self.histogram.type = self.histType[self.scoreMethods[self.methodIndex][2]]
            self.histogram.setValues(self.scores.values())
            self.histogram.setBoundary(*self.cuts.get(self.methodIndex, (0, 0)))
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
            scoresHi = scores[-max(self.selectNBest/2, 1):]
            scoresLo = scores[:max(self.selectNBest/2+self.selectNBest%2, 1)]
            self.histogram.setBoundary(scoresLo[-1][1], scoresHi[0][1])
        
    def Commit(self):
        if self.data and self.data.domain.classVar:
            cutOffUpper = self.histogram.upperBoundary
            cutOffLower = self.histogram.lowerBoundary
            test = self.scoreMethods[self.methodIndex][2]
            selectedAttrs = [attr for attr in self.data.domain.attributes if  test(attr, cutOffLower, cutOffUpper)] #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]
            self.send("Examples with selected attributes", self.data.select(selectedAttrs+[self.data.domain.classVar]))
            remainingAttrs = [attr for attr in self.data.domain.attributes if  not test(attr, cutOffLower, cutOffUpper)] #self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]
            self.send("Examples with remaining attributes", self.data.select(remainingAttrs+[self.data.domain.classVar]))
        else:
            self.send("Examples with selected attributes", None)
            self.send("Examples with remaining attributes", None)

if __name__=="__main__":
    import sys
    app = QApplication(sys.argv)
    data = orange.ExampleTable("E://leukemia.tab")
    w = OWFeatureSelection()
    w.show()
    w.SetData(data)
    app.setMainWidget(w)
    app.exec_loop()
    w.saveSettings()