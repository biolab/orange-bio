"""
<name>Vulcano Plot</name>
<description>Plots fold change vs. p-value.)</description>
<priority>1020</priority>
<contact>Ales Erjavec (ales.erjavec@fri.uni-lj.si)</contact>
<icon>icons/VulcanoPlot.png</icon>
"""

from OWWidget import *
from OWGraph import *
import OWGUI
import orange
from math import log
from statc import mean, ttest_ind
from obiGEO import transpose
import obiExpression

from orngDataCaching import data_hints

from OWToolbars import ZoomSelectToolbar

class VulcanoGraph(OWGraph):
    def __init__(self, master, *args, **kwargs):
        OWGraph.__init__(self, *args, **kwargs)
        self.master = master
        self.cutoffX = 2.0
        self.cutoffY = 3.0
        self.symbolSize = 5
        
        self.selectedCurve = self.addCurve("", brushColor=Qt.red)
        self.unselectedCurve = self.addCurve("", brushColor=Qt.blue)
        
        self.leftSelectionCurve = self.addCurve("", style=QwtPlotCurve.Lines, penColor=Qt.red, symbol=QwtSymbol.NoSymbol)
        
        self.rightSelectionCurve = self.addCurve("", style=QwtPlotCurve.Lines,  penColor=Qt.red, symbol=QwtSymbol.NoSymbol)
        
        self.plotValues = {}

        self.setAxisAutoScale(QwtPlot.xBottom)
        self.setAxisAutoScale(QwtPlot.yLeft)

        self.updatingAxes = False, None

#    def splitSelected(self):
#        items = self.plotValues.items()
#        return ([key for key, (x, y) in items if abs(x) >= self.cutoffX and y >= self.cutoffY],
#                [key for key, (x, y) in items if abs(x) < self.cutoffX or y < self.cutoffY])
        
    def splitSelected(self):
        test = (numpy.abs(self.plotData[:, 0]) >= self.cutoffX) & (self.plotData[:, 1] >= self.cutoffY)
        return (self.plotData[numpy.nonzero(test)], self.plotData[numpy.nonzero(~test)])
                
    
    def setPlotValues(self, values):
        self.plotValues = values
        self.plotData = numpy.array(values.values())
        self.replot_(setScale=True)

    def updateSelectionArea(self):
        x = numpy.array([self.maxX, self.cutoffX, self.cutoffX])
        y = numpy.array([self.cutoffY, self.cutoffY, self.maxY])
        self.leftSelectionCurve.setData(x, y)
        self.rightSelectionCurve.setData(-x, y)

    def replot_(self, setScale=False):
        if self.plotValues:
            data = self.plotData #numpy.array(self.plotValues.values())
            self.maxX = numpy.max(numpy.abs(data[:,0]))
            self.maxY = numpy.max(data[:, 1])
            if setScale:
                self.setAxisScale(QwtPlot.xBottom, -self.maxX, self.maxX)
                self.setAxisScale(QwtPlot.yLeft, 0.0, self.maxY)
            
            selected, unselected = self.splitSelected()
            getData = lambda keys, dim: [self.plotValues[attr][dim] for attr in keys]
            
#            self.selectedCurve.setData(getData(selected, 0), getData(selected, 1))
            self.selectedCurve.setData(selected[:,0], selected[:,1])
            self.selectedCurve.setBrush(QBrush(Qt.blue))
#            self.unselectedCurve.setData(getData(unselected, 0), getData(unselected, 1))
            self.unselectedCurve.setData(unselected[:, 0], unselected[:, 1])
            self.updateSelectionArea()
            self.master.infoLabel2.setText("%i selected genes" % len(selected))
        else:
            for curve in [self.selectedCurve, self.unselectedCurve, self.leftSelectionCurve, self.rightSelectionCurve]:
                curve.setData([],[])
            self.master.infoLabel2.setText("0 selected genes")
        self.replot()

    def getSelectionAxes(self, pos):
        cx1 = self.transform(QwtPlot.xBottom, self.cutoffX)
        cx2 = self.transform(QwtPlot.xBottom, -self.cutoffX)
        cy = self.transform(QwtPlot.yLeft, self.cutoffY)
        x = self.canvas().mapFrom(self, pos).x()
        y = self.canvas().mapFrom(self, pos).y()
##        print cx1, cx2, x, cy, y
        offset = 3
        if y < cy - offset:
            if min(abs(x - cx1), abs(x - cx2)) < offset:
                return QwtPlot.xBottom
        if x > cx1 + offset or x < cx2 - offset:
            if abs(y - cy) < offset:
                return QwtPlot.yLeft
        if abs(y - cy) <= offset and min(abs(x - cx1), abs(x - cx2)) <= offset:
            return QwtPlot.xBottom ## TODO both axes
        else:
            return -1

    def updateCutoff(self, axis, pos):
        if axis == QwtPlot.xBottom:
            self.cutoffX = abs(self.invTransform(axis, self.canvas().mapFrom(self, pos).x()))
        elif axis == QwtPlot.yLeft:
            self.cutoffY = self.invTransform(axis, self.canvas().mapFrom(self, pos).y())
        else:
            self.cutoffX = abs(self.invTransform(QwtPlot.xBottom, self.canvas().mapFrom(self, pos).x()))
            self.cutoffY = self.invTransform(QwtPlot.yLeft, self.canvas().mapFrom(self, pos).y())
        self.replot_()

    def mousePressEvent(self, event):
        if self.state == SELECT:
            axes = self.getSelectionAxes(event.pos())
            self.updateCutoff(axes, event.pos())
            self.updatingAxes = True, axes
##            print axes
        else:
            OWGraph.mousePressEvent(self, event)

    def mouseMoveEvent(self, event):
        if self.state == SELECT:
            drag, axes = self.updatingAxes
            if drag:
                self.updateCutoff(axes, event.pos())
            else:
                axes = self.getSelectionAxes(event.pos())
                if axes == QwtPlot.xBottom:
                    self.canvas().setCursor(Qt.SizeHorCursor)
                elif axes == QwtPlot.yLeft:
                    self.canvas().setCursor(Qt.SizeVerCursor)
                else:
                    self.canvas().setCursor(self._cursor)
        else:
            OWGraph.mouseMoveEvent(self, event)

    def mouseReleaseEvent(self, event):
        if self.state == SELECT:
            self.updatingAxes = False, -1
            self.master.commitIf()
        else:
            OWGraph.mouseReleaseEvent(self, event)

    def updateSymbolSize(self):
##        self.selectedCurve.setSymbol(QwtSymbol(QwtSymbol.Ellipse, 
        self.selectedCurve.symbol().setSize(self.symbolSize)
        self.unselectedCurve.symbol().setSize(self.symbolSize)
        self.replot()

class OWVulcanoPlot(OWWidget):
    settingsList =["targetClass", "graph.cutoffX", "graph.cutoffY", "graph.symbolSize", "showXTitle", "showYTitle"]
    contextHandlers = {"":DomainContextHandler("", [ContextField("targetClass"), ContextField("graph.cutoffX"),
                                                    ContextField("graph.cutoffY")])}
    def __init__(self, parent=None, signalManager=None, name="Vulcano Plot"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True)
        
        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs =[("Examples with selected attributes", ExampleTable)]

        self.genesInColumns = False
        self.targetClass = 0

        self.showXTitle = True
        self.showYTitle = True

        self.autoCommit = False
        self.selectionChangedFlag = False

        self.graph = VulcanoGraph(self)
        self.mainArea.layout().addWidget(self.graph)

        ## GUI
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.infoLabel = OWGUI.label(box, self, "")
        self.infoLabel.setText("No data on input\n")
        self.infoLabel2 = OWGUI.label(box, self, "")
        self.infoLabel2.setText("0 selected genes")
        
        box = OWGUI.widgetBox(self.controlArea, "Target Label (Class)")
        self.genesInColumnsCheck = OWGUI.checkBox(box, self, "genesInColumns", "Genes in columns", callback=[self.setTargetCombo, self.plot])
        self.targetClassCombo = OWGUI.comboBox(box, self, "targetClass", callback=self.plot)

        box = OWGUI.widgetBox(self.controlArea, "Settings")
        OWGUI.hSlider(box, self, "graph.symbolSize", label="Symbol size:   ", minValue=2, maxValue=20, step=1, callback = self.graph.updateSymbolSize)
        OWGUI.checkBox(box, self, "showXTitle", "X axis title", callback=self.setAxesTitles)
        OWGUI.checkBox(box, self, "showYTitle", "Y axis title", callback=self.setAxesTitles)
        
        ZoomSelectToolbar(self, self.controlArea, self.graph, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])

        box = OWGUI.widgetBox(self.controlArea, "Commit")
        b = OWGUI.button(box, self, "Commit", callback=self.commit)
        cb = OWGUI.checkBox(box, self, "autoCommit", "Commit automatically")
        OWGUI.setStopper(self, b, cb, "selectionChangedFlag", self.commitIf)

        self.connect(self.graphButton, SIGNAL("clicked()"), self.graph.saveToFile)
        
        OWGUI.rubber(self.controlArea)

        self.data = None
        
        self.resize(800, 600)

    def setData(self, data=None):
        self.closeContext()
        self.data = data
        self.targetClassCombo.clear()
        self.targetClass = 0
        self.error(0)
        if data:
            self.genesInColumns = not bool(data.domain.classVar)
            self.genesInColumnsCheck.setDisabled(not bool(data.domain.classVar))
            if self.genesInColumns:
                self.genesInColumns = not data_hints.get_hint(data, "genesinrows", not self.genesInColumns) 
            self.setTargetCombo()
            self.error()
            if not self.targets:
                self.error(0, "Data set with no column labels (attribute tags) or row labels (classes).")
        else:
            self.infoLabel.setText("No data on input\n")
            self.targets = []
        self.openContext("", data)
        self.plot()
        
    def setTargetCombo(self):
        if self.genesInColumns:
            self.targets = sorted(reduce(set.union, [attr.attributes.items() for attr in (self.data.domain.attributes if self.data else [])], set()))
            measurements = [attr.attributes.items() for attr in (self.data.domain.attributes if self.data else [])]
            targets = ["%s: %s" % t for t in self.targets]
            
        else:
            self.targets = list(self.data.domain.classVar.values if self.data else [])
            measurements = [set([str(ex.getclass())]) for ex in (self.data if self.data else [])]
            targets = self.targets
                                
        self.targetMeasurements = [len([m for m in measurements if target in m]) for target in self.targets]
        
        self.targetClassCombo.clear()
        self.targetClassCombo.addItems(targets)
#        self.targetClass = min(self.targetClass, len(self.targets) - 1)
                             
    def plot(self):
##        self.graph.clear()
        self.values = {}
        if self.data and self.targets:
            self.warning(0)
            targetClassIndex = min(self.targetClass, len(self.targets) - 1)
            if self.targetMeasurements[targetClassIndex] < 2 or sum(self.targetMeasurements) - self.targetMeasurements[targetClassIndex] < 2:
                self.warning(0, "Insufficient data to compute statistics. More than one measurement per class should be provided")
            targetClass = self.targets[targetClassIndex]
            self.progressBarInit()
            tt = obiExpression.ExpressionSignificance_TTest(self.data, useAttributeLabels=self.genesInColumns)(targetClass)
            self.progressBarSet(25)
            fold = obiExpression.ExpressionSignificance_FoldChange(self.data, useAttributeLabels=self.genesInColumns)(targetClass)
            self.progressBarSet(50)
            self.infoLabel.setText("%i genes on input" % len(fold))
            import numpy
            invalid = set([key for (key, (t, p)), (_, f) in zip(tt, fold) if any(v is numpy.ma.masked for v in [t, p, f]) or f==0.0])
            tt = [t for t in tt if t[0] not in invalid]
            fold = [f for f in fold if f[0] not in invalid]
            self.progressBarSet(75)
            logratio = numpy.log2(numpy.abs([v for k, v in fold]))
            logpval = -numpy.log10([p for k, (t, p) in tt])
            self.values = dict(zip([k for k, v in tt], zip(logratio, logpval)))
            self.progressBarFinished()
        self.graph.setPlotValues(self.values)
        self.setAxesTitles()
        self.updateTooltips()

    def setAxesTitles(self):
        self.graph.setAxisTitle(QwtPlot.xBottom, "log2 (ratio)" if self.showXTitle else "")
        self.graph.setAxisTitle(QwtPlot.yLeft, "-log10 (p_value)" if self.showYTitle else "")

    def updateTooltips(self):
        self.graph.tips.removeAll()
        for key, (logratio, logpval) in self.values.items():
            self.graph.tips.addToolTip(logratio, logpval, "<b>%s</b><hr>log2(ratio): %.5f<br>p-value: %.5f" \
                                       %(str(key) if self.genesInColumns else key.name, logratio, math.pow(10, -logpval)))

    def commit(self):
#        def passAttributes(src, dst, names):
#            for name in names:
#                if hasattr(src, name):
#                    setattr(dst, name, getattr(src, name))
                    
        check = lambda x,y:abs(x) >= self.graph.cutoffX and y >= self.graph.cutoffY
        if self.data and self.genesInColumns:
            selected = [self.data[i] for i in range(len(self.data)) if i in self.values and check(*self.values[i])]
            if selected:
                data = orange.ExampleTable(self.data.domain, selected)
            else:
                data = None
        elif self.data:
            check = lambda x,y:abs(x) >= self.graph.cutoffX and y >= self.graph.cutoffY
            selected = [attr for attr in self.data.domain.attributes if attr in self.values and check(*self.values[attr]) or attr.varType==orange.VarTypes.String]
            newdomain = orange.Domain(selected + [self.data.domain.classVar])
            newdomain.addmetas(self.data.domain.getmetas())
            data = orange.ExampleTable(newdomain, self.data)
        else:
            data = None
#            if self.transposedData:
#                data = transpose(data)
#        if data:
#            passAttributes(self.data, data, ["taxid", "genesinrows"])
        self.send("Examples with selected attributes", data)
        self.selectionChangedFlag = False

    def commitIf(self):
        if self.autoCommit:
            self.commit()
        else:
            self.selectionChangedFlag = True
        
if __name__ == "__main__":
    ap = QApplication(sys.argv)
    w = OWVulcanoPlot()
##    d = orange.ExampleTable("E:\\affy(HD-CC)_GS_C2cpC5.tab")
#    d = orange.ExampleTable("E:\\steroltalk-smallchip.tab")
    d = orange.ExampleTable("../../../doc/datasets/brown-selected.tab")
    w.setData(d)
    w.show()
    ap.exec_()
    w.saveSettings()

                
        
