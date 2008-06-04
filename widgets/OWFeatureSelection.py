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

import OWGUI

class ScoreHist(OWInteractiveHist):
    def __init__(self, master, parent=None):
        OWInteractiveHist.__init__(self, parent)
        self.master = master
        self.setAxisTitle(QwtPlot.xBottom, "Score")
        self.setAxisTitle(QwtPlot.yLeft, "Frequency")

    def updateData(self):
        OWInteractiveHist.updateData(self)
##        self.upperTailShadeKey = self.insertCurve(PolygonCurve(self, pen=QPen(Qt.blue), brush=QBrush(Qt.red))) # self.addCurve("upperTailShade", Qt.blue, Qt.blue, 6, symbol = QwtSymbol.None, style = QwtPlotCurve.Sticks)
        self.upperTailShadeKey = self.addCurve("upperTailShade", Qt.blue, Qt.blue, 6, symbol = QwtSymbol.None, style = QwtPlotCurve.Steps)
##        self.lowerTailShadeKey = self.insertCurve(PolygonCurve(self, pen=QPen(Qt.blue), brush=QBrush(Qt.red))) #self.addCurve("lowerTailShade", Qt.blue, Qt.blue, 6, symbol = QwtSymbol.None, style = QwtPlotCurve.Sticks)
        self.lowerTailShadeKey = self.addCurve("lowerTailShade", Qt.blue, Qt.blue, 6, symbol = QwtSymbol.None, style = QwtPlotCurve.Steps)
##        self.setCurveStyle(self.upperTailShadeKey, QwtPlotCurve.Steps)
##        self.setCurveStyle(self.lowerTailShadeKey, QwtPlotCurve.Steps)
        self.setCurveBrush(self.upperTailShadeKey, QBrush(Qt.red))
        self.setCurveBrush(self.lowerTailShadeKey, QBrush(Qt.red))

    def shadeTails(self):
        index = max(min(int(100*(self.upperBoundary-self.minx)/(self.maxx-self.minx)), 100)-1, 0)
        x = [self.upperBoundary] + list(self.xData[index:])
        y = [self.yData[index]] + list(self.yData[index:])
        self.setCurveData(self.upperTailShadeKey, x, y)
        if self.type == 1:
            index = max(min(int(100*(self.lowerBoundary-self.minx)/(self.maxx-self.minx)),99), 0)
            x = list(self.xData[:index]) + [self.lowerBoundary]
            y = list(self.yData[:index]) + [self.yData[index]]
            self.setCurveData(self.lowerTailShadeKey, x, y)
        else:
            self.setCurveData(self.lowerTailShadeKey, [], [])
        
    def setBoundary(self, low, hi):
        OWInteractiveHist.setBoundary(self, low, hi)
        self.master.UpdateSelectedInfoLabel(low, hi)
        self.shadeTails()
        self.replot()
        if self.master.autoCommit:
            self.master.Commit()
        
class OWFeatureSelection(OWWidget):
    settingsList=["methodIndex", "autoCommit"]
##    contextHandlers={"":DomainContextHandler("",[])}
    def __init__(self, parent=None, signalManager=None, name="Feature selection"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Examples", ExampleTable, self.SetData)]
        self.outputs = [("Examples with selected attributes", ExampleTable), ("Examples with remaining attributes", ExampleTable)]

        self.methodIndex = 0
        self.autoCommit = False
        self.selectNBest = 20
##        self.infoLabel = "No data on input"        

        self.oneTailTest = oneTailTest = lambda attr, low, hi:self.scores.get(attr,0)>=hi
        self.twoTailTest = twoTailTest = lambda attr, low, hi:self.scores.get(attr,0)>=hi or self.scores.get(attr,0)<=low
        self.scoreMethods = [("chi-squared", orange.MeasureAttribute_chiSquare, oneTailTest),
                             ("info gain", orange.MeasureAttribute_info, oneTailTest),
                             ("signal to noise ratio", lambda attr, data:MA_signalToNoise()(attr, data), twoTailTest),
                             ("t-test",lambda attr, data:MA_t_test()(attr, data), twoTailTest),
                             ("t-test p-value",lambda attr, data:MA_t_test(prob=True)(attr, data), oneTailTest),
                             ("fold test", lambda attr, data:MA_fold_test()(attr, data), twoTailTest),
                             ("log2 fold test", lambda attr, data:math.log(max(min(MA_fold_test()(attr, data), 1e300), 1e-300)), twoTailTest)]

        boxLayout = QVBoxLayout(self.mainArea)
        self.histogram = ScoreHist(self, self.mainArea)
        boxLayout.addWidget(self.histogram)
        self.histogram.setSizePolicy(QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum))
        self.histogram.show()
        
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.dataInfoLabel = OWGUI.widgetLabel(box, "\n\n")
        self.selectedInfoLabel = OWGUI.widgetLabel(box, "")
        OWGUI.radioButtonsInBox(self.controlArea, self, "methodIndex", [sm[0] for sm in self.scoreMethods], box="Score Method", callback=self.Update, addSpace=True)
        box = OWGUI.widgetBox(self.controlArea, "Selection", addSpace=True)
        OWGUI.button(box, self, "Select n best features", callback=self.SelectNBest)
        OWGUI.spin(box, self, "selectNBest", 0, 10000, step=1, label="n:")
        OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.button(box, self, "&Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)
        
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
            self.histogram.type = 0 if self.scoreMethods[self.methodIndex][2]==self.oneTailTest else 1
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
        if self.scoreMethods[self.methodIndex][2]==self.oneTailTest:
            scores = scores[-max(self.selectNBest, 1):]
            self.histogram.setBoundary(scores[0][1], scores[0][1])
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