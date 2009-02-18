"""
<name>Vulcano Plot</name>
<description>Vulcano plot</description>
<priority>231</priority>
<contact>Ales Erjavec (ales.erjavec@fri.uni-lj.si)</contact>
"""

from OWWidget import *
from OWGraph import *
import OWGUI
import orange
from math import log
from statc import mean, ttest_ind

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

    def splitSelected(self):
        items = self.plotValues.items()
        return ([key for key, (x, y) in items if abs(x) >= self.cutoffX and y >= self.cutoffY],
                [key for key, (x, y) in items if abs(x) < self.cutoffX or y < self.cutoffY])
    
    def setPlotValues(self, values):
        self.plotValues = values
        self.replot_(setScale=True)

    def updateSelectionArea(self):
        x = numpy.array([self.maxX, self.cutoffX, self.cutoffX])
        y = numpy.array([self.cutoffY, self.cutoffY, self.maxY])
        self.leftSelectionCurve.setData(x, y)
        self.rightSelectionCurve.setData(-x, y)

    def replot_(self, setScale=False):
        if self.plotValues:
            data = numpy.array(self.plotValues.values())
            self.maxX = numpy.max(numpy.abs(data[:,0]))
            self.maxY = numpy.max(data[:, 1])
            if setScale:
                self.setAxisScale(QwtPlot.xBottom, -self.maxX, self.maxX)
                self.setAxisScale(QwtPlot.yLeft, 0.0, self.maxY)
            
            selected, unselected = self.splitSelected()
            getData = lambda keys, dim: [self.plotValues[attr][dim] for attr in keys]
            
            self.selectedCurve.setData(getData(selected, 0), getData(selected, 1))
            self.selectedCurve.setBrush(QBrush(Qt.blue))
            self.unselectedCurve.setData(getData(unselected, 0), getData(unselected, 1))
            self.updateSelectionArea()
        else:
            for curve in [self.selectedCurve, self.unselectedCurve, self.leftSelectionCurve, self.rightSelectionCurve]:
                curve.setData([],[])
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
        
        self.inputs = [("Example Table", ExampleTable, self.setData)]
        self.outputs =[("Example Table", ExampleTable)]

        self.targetClass = 0

        self.showXTitle = True
        self.showYTitle = True

        self.autoCommit = False        

        self.graph = VulcanoGraph(self)
        self.mainArea.layout().addWidget(self.graph)

        ## GUI
        self.infoLabel = OWGUI.label(OWGUI.widgetBox(self.controlArea, "Info"), self, "")
        self.infoLabel.setText("No data on input\n")

        self.targetClassCombo = OWGUI.comboBox(self.controlArea, self, "targetClass", "Target Class", callback=self.plot)

        box = OWGUI.widgetBox(self.controlArea, "Settings")
        OWGUI.hSlider(box, self, "graph.symbolSize", label="Symbol size:   ", minValue=2, maxValue=20, step=1, callback = self.graph.updateSymbolSize)
        OWGUI.checkBox(box, self, "showXTitle", "X axis title", callback=self.setAxesTitles)
        OWGUI.checkBox(box, self, "showYTitle", "Y axis title", callback=self.setAxesTitles)
        
        ZoomSelectToolbar(self, self.controlArea, self.graph, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])

        box = OWGUI.widgetBox(self.controlArea, "Commit")
        OWGUI.button(box, self, "Commit", callback=self.commit)
        OWGUI.checkBox(box, self, "autoCommit", "Commit automatically")

        self.connect(self.graphButton, SIGNAL("clicked()"), self.graph.saveToFile)
        
        OWGUI.rubber(self.controlArea)

        self.resize(600, 600)

    def setData(self, data=None):
        self.closeContext()
        self.data = data
        self.targetClassCombo.clear()
        self.targetClass = 0
        self.error(0)
        if data and data.domain.classVar:
            self.targetClassCombo.addItems([value for value in data.domain.classVar.values])
            self.infoLabel.setText("Genes: %i\nSamples: %i" %(len(data.domain.attributes), len(data)))
        elif data:
            self.infoLabel.setText("No data on input\n")
            self.error(0, "Class-labeled data set required.")
            self.data = None
        else:
            self.infoLabel.setText("No data on input\n")
        self.openContext("", data)
        self.plot()

    def plot(self):
##        self.graph.clear()
        self.values = {}
        if self.data:
            targetClass = self.data.domain.classVar(self.targetClass)
            self.progressBarInit()
            milestones = set(range(0, len(self.data.domain.attributes), max(len(self.data.domain.attributes)/100, 1)))
            for i, attr in enumerate(self.data.domain.attributes):
                sample1 = [float(ex[attr]) for ex in self.data if not ex[attr].isSpecial() and ex.getclass()==targetClass]
                sample2 = [float(ex[attr]) for ex in self.data if not ex[attr].isSpecial() and ex.getclass()!=targetClass]
                logratio = log(abs(mean(sample1)/mean(sample2)), 2)
                t, pval = ttest_ind(sample1, sample2)
                logpval = -log(pval, 10)
                self.values[attr] = (logratio, logpval)
                if i in milestones:
                    self.progressBarSet(100.0*i/len(self.data.domain.attributes))
            self.progressBarFinished()
        self.graph.setPlotValues(self.values)
        self.setAxesTitles()
        self.updateTooltips()

    def setAxesTitles(self):
        self.graph.setAxisTitle(QwtPlot.xBottom, "log2 (ratio)" if self.showXTitle else "")
        self.graph.setAxisTitle(QwtPlot.yLeft, "-log10 (p_value)" if self.showYTitle else "")

    def updateTooltips(self):
        self.graph.tips.removeAll()
        for attr, (logratio, logpval) in self.values.items():
            self.graph.tips.addToolTip(logratio, logpval, "<b>%s</b><hr>log2(ratio): %.5f<br>p-value: %.5f" \
                                       %(attr.name, logratio, math.pow(10, -logpval)))

    def commit(self):
        if self.data:
            check = lambda x,y:abs(x) >= self.graph.cutoffX and y >= self.graph.cutoffY
            selected = [attr for attr in self.data.domain.attributes if check(*self.values[attr])]
            newdomain = orange.Domain(selected + [self.data.domain.classVar])
            newdomain.addmetas(self.data.domain.getmetas())
            data = orange.ExampleTable(newdomain, self.data)
            self.send("Example Table", data)

    def commitIf(self):
        if self.autoCommit:
            self.commit()
        
if __name__ == "__main__":
    ap = QApplication(sys.argv)
    w = OWVulcanoPlot()
    d = orange.ExampleTable("E:\\affy(HD-CC)_GS_C2cpC5.tab")
##    d = orange.ExampleTable("../../orange/doc/datasets/brown-selected.tab")
    w.setData(d)
    w.show()
    ap.exec_()
    w.saveSettings()

                
        