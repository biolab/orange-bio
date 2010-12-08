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

class GraphSelections(QObject):
    def __init__(self, parent):
        QObject.__init__(self, parent)
        self.selection = []
        
    def getPos(self, event):
        graph = self.parent()
        pos = graph.canvas().mapFrom(graph, event.pos())
        x = graph.invTransform(QwtPlot.xBottom, pos.x())
        y = graph.invTransform(QwtPlot.yLeft, pos.y())
        return QPointF(x, y)
        
    def start(self, event):
        pos = self.getPos(event)
        if event.modifiers() & Qt.ControlModifier:
            self.selection.append((pos, pos))
        else:
            self.selection = [(pos, pos)]
        self.emit(SIGNAL("selectionGeometryChanged()"))
    
    def update(self, event):
        pos = self.getPos(event)
        self.selection[-1] = self.selection[-1][:-1] + (pos,)
        self.emit(SIGNAL("selectionGeometryChanged()"))
    
    def end(self, event):
        self.update(event)
        
    def testSelection(self, data):
        if len(data) == 0:
            return []
        data = numpy.asarray(data)
        region = QPainterPath()
        for p1, p2 in self.selection:
            region.addRect(QRectF(p1, p2).normalized())
        def test(point):
            return region.contains(QPointF(point[0], point[1]))
        test = numpy.apply_along_axis(test, 1, data)
        return test
        
class SymetricSelections(GraphSelections):
    def __init__(self, parent, x=3, y=3):
        GraphSelections.__init__(self, parent)
        max = 10000
        self.selection = [(QPointF(-max, max), QPointF(-x, y)), (QPointF(max, max), QPointF(x, y))]
        self.updateAxes = None
        
    def updateSelection(self, axes, pos):
        if axes == QwtPlot.xBottom or axes == -1:
            self.selection[0][1].setX(-abs(pos.x()))
            self.selection[1][1].setX(abs(pos.x()))
        if axes == QwtPlot.yLeft or axes == -1:
            self.selection[0][1].setY(pos.y())
            self.selection[1][1].setY(pos.y())
            
        self.emit(SIGNAL("selectionGeometryChanged()"))
        
    def getAxesAndPos(self, event):
        graph = self.parent()
        pos = graph.canvas().mapFrom(graph, event.pos())
        x = graph.invTransform(QwtPlot.xBottom, pos.x())
        y = graph.invTransform(QwtPlot.yLeft, pos.y())
        
        offset = 3
        dx = abs(graph.invTransform(QwtPlot.xBottom, pos.x() + offset) - x)
        dy = abs(graph.invTransform(QwtPlot.yLeft, pos.y() + offset) - y)
        
        x = abs(x)
        
        cx = self.selection[1][1].x()
        cy = self.selection[1][1].y()

        bottom = QRectF(QPointF(cx, cy), QPointF(graph.maxX, cy)).adjusted(-dx, dy, dx, -dy).normalized()
        left = QRectF(QPointF(cx, graph.maxY), QPointF(cx, cy)).adjusted(-dx, dy, dx, -dy).normalized()
        
        if bottom.contains(QPointF(x, y)) or bottom.contains(QPointF(-x, y)):
            axes = QwtPlot.yLeft
        elif left.contains(QPointF(x, y)) or left.contains(QPointF(-x, y)):
            axes = QwtPlot.xBottom
        else:
            axes = -1
        return axes, QPointF(x, y)
        
    def start(self, event):
        axes, pos = self.getAxesAndPos(event)
        self.updateAxes = axes
        self.updateSelection(axes, pos)
        
    def update(self, event):
        _, pos = self.getAxesAndPos(event)
        self.updateSelection(self.updateAxes, pos)
    
    def end(self, event):
        self.update(event)
        self.updateAxes = None
        
    def testSelection(self, data):
        if len(data) == 0:
            return []
        data = numpy.asarray(data)
        cutoffX = self.selection[1][1].x()
        cutoffY = self.selection[1][1].y()
        return (numpy.abs(data[:, 0]) >= cutoffX) & (data[:, 1] >= cutoffY)
    
class VulcanoGraph(OWGraph):
    def __init__(self, master, *args, **kwargs):
        OWGraph.__init__(self, *args, **kwargs)
        self.master = master
        self.cutoffX = 2.0
        self.cutoffY = 3.0
        self.maxX, self.maxY = 10, 10
        self.symbolSize = 5
        self.symetricSelections = True
        
        self.selectedCurve = self.addCurve("", brushColor=Qt.red)
        self.unselectedCurve = self.addCurve("", brushColor=Qt.blue)
        
        self.plotValues = {}

        self.setAxisAutoScale(QwtPlot.xBottom)
        self.setAxisAutoScale(QwtPlot.yLeft)
        
        self.reselect(replot=False)
        
    def setSelection(self, selection):
        self.selection = selection
        self.connect(self.selection, SIGNAL("selectionGeometryChanged()"), self.onSelectionChanged)
        if self.plotValues:
            self.updateSelectionArea()
            
    def onSelectionChanged(self):
        self.replot_()
        
    def splitSelected(self):
        test =  self.selection.testSelection(self.plotData)
        return (self.plotData[numpy.nonzero(test)], self.plotData[numpy.nonzero(~test)])
    
    def setPlotValues(self, values):
        self.plotValues = values
        self.plotData = numpy.array(values.values())
        self.replot_(setScale=True)

    def createSelectionRectCurve(self, p1, p2):
        curve = self.addCurve("selection", style=QwtPlotCurve.Lines, penColor=Qt.red, symbol=QwtSymbol.NoSymbol)
        curve.setData([p1.x(), p2.x(), p2.x(), p1.x(), p1.x()], [p1.y(), p1.y(), p2.y(), p2.y(), p1.y()])
        
    def items(self, title=None):
        for item in self.itemList():
            if str(item.title().text()) == title:
                yield item
        
    def updateSelectionArea(self):
        for c in self.items(title="selection"):
            c.detach()
        for p1, p2 in self.selection.selection:
            self.createSelectionRectCurve(p1, p2)

    def replot_(self, setScale=False):
        if self.plotValues:
            data = self.plotData
            self.maxX = numpy.max(numpy.abs(data[:,0]))
            self.maxY = numpy.max(data[:, 1])
            if setScale:
                self.setAxisScale(QwtPlot.xBottom, -self.maxX, self.maxX)
                self.setAxisScale(QwtPlot.yLeft, 0.0, self.maxY)
            
            selected, unselected = self.splitSelected()
            getData = lambda keys, dim: [self.plotValues[attr][dim] for attr in keys]
            
            self.selectedCurve.setData(selected[:,0], selected[:,1])
            self.selectedCurve.setBrush(QBrush(Qt.blue))
            self.unselectedCurve.setData(unselected[:, 0], unselected[:, 1])
            self.updateSelectionArea()
            self.master.infoLabel2.setText("%i selected genes" % len(selected))
        else:
            for curve in [self.selectedCurve, self.unselectedCurve]:
                curve.setData([],[])
            self.master.infoLabel2.setText("0 selected genes")
        self.replot()

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
            if event.button() == Qt.LeftButton:
                self.selection.start(event)
        else:
            OWGraph.mousePressEvent(self, event)

    def mouseMoveEvent(self, event):
        if self.state == SELECT:
            if event.buttons() & Qt.LeftButton:
                self.selection.update(event)
            if isinstance(self.selection, SymetricSelections):
                axes, pos = self.selection.getAxesAndPos(event)
                cursors = {QwtPlot.xBottom: Qt.SizeHorCursor,
                           QwtPlot.yLeft: Qt.SizeVerCursor}
                self.canvas().setCursor(cursors.get(axes, self._cursor))
        else:
            OWGraph.mouseMoveEvent(self, event)

    def mouseReleaseEvent(self, event):
        if self.state == SELECT:
            if event.button() == Qt.LeftButton:
                self.selection.end(event)
        else:
            OWGraph.mouseReleaseEvent(self, event)
            
    def reselect(self, replot=True):
        if self.symetricSelections:
            self.setSelection(SymetricSelections(self, x=self.maxX*0.80, y=self.maxY*0.80))
        else:
            self.setSelection(GraphSelections(self))
            self.canvas().setCursor(self._cursor)
        if replot:
            self.replot_()

    def updateSymbolSize(self):
        def setSize(curve, size):
            symbol = curve.symbol()
            symbol.setSize(size)
            if QWT_VERSION_STR >= "5.2":
                curve.setSymbol(symbol)
        setSize(self.selectedCurve, self.symbolSize)
        setSize(self.unselectedCurve, self.symbolSize)
        self.replot()

class OWVulcanoPlot(OWWidget):
    settingsList =["targetClass", "graph.cutoffX", "graph.cutoffY", "graph.symbolSize", "graph.symetricSelections", "showXTitle", "showYTitle"]
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
        
        toolbar = ZoomSelectToolbar(self, self.controlArea, self.graph, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])
        
        top_layout = toolbar.layout()
        top_layout.setDirection(QBoxLayout.TopToBottom)
        button_layotu = QHBoxLayout()
        top_layout.insertLayout(0, button_layotu)
        
        for i in range(1, top_layout.count()):
            item = top_layout.itemAt(1)
            top_layout.removeItem(item)
            button_layotu.addItem(item)
        
        
        OWGUI.checkBox(toolbar, self, "graph.symetricSelections", "Symetric selection", callback=self.graph.reselect)

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
                             
    def plot(self):
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
        self.graph.setAxisTitle(QwtPlot.xBottom, "log<sub>2</sub> (ratio)" if self.showXTitle else "")
        self.graph.setAxisTitle(QwtPlot.yLeft, "-log<sub>10</sub> (p_value)" if self.showYTitle else "")

    def updateTooltips(self):
        self.graph.tips.removeAll()
        for key, (logratio, logpval) in self.values.items():
            self.graph.tips.addToolTip(logratio, logpval, "<b>%s</b><hr>log<sub>2</sub>(ratio): %.5f<br>p-value: %.5f" \
                                       %(str(key) if self.genesInColumns else key.name, logratio, math.pow(10, -logpval)))

    def commit(self):
        if self.data and self.genesInColumns:
            items = sorted(self.values.items())
            test = self.graph.selection.testSelection([val for key, val in items])
            selected = [self.data[i] for t, (i, value) in zip(test, items) if t]
            if selected:
                data = orange.ExampleTable(self.data.domain, selected)
            else:
                data = None
        elif self.data:
            attrs = [(attr, self.values[attr])  for attr in self.data.domain.attributes if attr in self.values]
            test = self.graph.selection.testSelection([val for attr, val in attrs])
            selected = [attr for t, (attr, val) in zip(test, attrs) if t]
            newdomain = orange.Domain(selected + [self.data.domain.classVar])
            newdomain.addmetas(self.data.domain.getmetas())
            data = orange.ExampleTable(newdomain, self.data)
        else:
            data = None
        
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

                
        
