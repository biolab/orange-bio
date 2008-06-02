"""
<name>FeatureSelection</name>
<description>Feature selection</description>
"""
import orange

from obiExpression import MA_signalToNoise, MA_t_test

from OWGraph import *
from OWWidget import *
from OWHist import OWHist

import OWGUI

class ScoreHist(OWHist):
    def __init__(self, master, parent=None):
        OWHist.__init__(self, parent)
        self.master = master
        self.setAxisTitle(QwtPlot.xBottom, "Score")
        self.setAxisTitle(QwtPlot.yLeft, "Frequency")

    def setBoundary(self, button, cut):
        if button == QMouseEvent.LeftButton:
            low, hi = cut, self.upperBoundary
        else:
            low, hi = self.lowerBoundary, cut
        if low > hi:
            low, hi = hi, low
        if self.type=1:
            OWHist.setBoundary(self, low, hi)
        else:
            OWHist.setBoundary(self, cut, cut)
        self.master.UpdateSelectedInfoLabel(low, hi)
        if self.master.autoCommit:
            self.master.Commit()
        
    def onMousePressed(self, e):
        cut = self.invTransform(QwtPlot.xBottom, e.x())
        self.mouseCurrentlyPressed = 1
        self.setBoundary(e.button(), cut)
        
    def onMouseMoved(self, e):
##        OWHist.onMouseMoved(self, e)
        if self.mouseCurrentlyPressed:
            cut = self.invTransform(QwtPlot.xBottom, e.x())
            self.setBoundary(e.state(), cut)

    def onMouseReleased(self, e):
##        OWHist.onMouseRelesed(self, e)
        cut = self.invTransform(QwtPlot.xBottom, e.x())
        self.mouseCurrentlyPressed = 0
        self.setBoundary(e.button(), cut)

class OWFeatureSelection(OWWidget):
    settingsList=["methodIndex", "autoCommit"]
##    contextHandlers={"":DomainContextHandler("",[])}
    def __init__(self, parent=None, signalManager=None, name="Feature selection"):
        OWWidget.__init__(self, parent, signalManager, name)
        self.inputs = [("Examples", ExampleTable, self.SetData)]
        self.outputs = [("Examples", ExampleTable)]

        self.methodIndex = 0
        self.autoCommit = False
##        self.infoLabel = "No data on input"

        self.scoreMethods = [("chi-squared", orange.MeasureAttribute_chiSquare),
                             ("info gain", orange.MeasureAttribute_info),
                             ("signal to noise ratio", lambda attr, data:MA_signalToNoise()(attr, data)),
                             ("t-test",lambda attr, data:MA_t_test()(attr, data))]

        boxLayout = QVBoxLayout(self.mainArea)    
        self.histogram = ScoreHist(self, self.mainArea)
        boxLayout.addWidget(self.histogram)
        self.histogram.setSizePolicy(QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum))
        self.histogram.show()
        
        box = OWGUI.widgetBox(self.controlArea, "Info", addSpace=True)
        self.dataInfoLabel = OWGUI.widgetLabel(box, "\n\n")
        self.selectedInfoLabel = OWGUI.widgetLabel(box, "")
        OWGUI.radioButtonsInBox(self.controlArea, self, "methodIndex", [sm[0] for sm in self.scoreMethods], box="Score Method", callback=self.Update, addSpace=True)
        box = OWGUI.widgetBox(self.controlArea, self, "Selection", addSpace=True)
        OWGUI.checkBox(box, self, "autoCommit", "Commit on change")
        OWGUI.button(box, self, "&Commit", callback=self.Commit)
        OWGUI.rubber(self.controlArea)
        
        self.loadSettings()

        self.data = None
        self.discData = None
        self.scoreCache = {}
        self.discretizer = orange.EquiNDiscretization(numberOfIntervals=5)
        
    def SetData(self, data):
        self.error(0)
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
            self.histogram.type = self.methodIndex in [0, 1] and 0 or 1
            self.histogram.setValues(self.scores.values())
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
        if self.data:
            self.selectedInfoLabel.setText("%i selected attributes" %len([attr for attr in self.data.domain.attributes if self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]))
        else:
            self.selectedInfoLabel.setText("0 selected attributes")
            
    def Commit(self):
        if self.data and self.data.domain.classVar:
            cutOffUpper = self.histogram.upperBoundary
            cutOffLower = self.histogram.lowerBoundary
            selectedAttrs = [attr for attr in self.data.domain.attributes if self.scores.get(attr,0)>cutOffUpper or self.scores.get(attr,0)<cutOffLower]
            self.send("Examples", self.data.select(selectedAttrs+[self.data.domain.classVar]))
        else:
            self.send("Examples", None)

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