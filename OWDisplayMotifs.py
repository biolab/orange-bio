"""
<name>Display Motifs</name>
<description>None.</description>
<icon>icons\DisplayMotifs.png</icon>
"""

import orange
from OWTools import *
from OWWidget import *
from OWGraph import *
from OWGUI import *

class subBarQwtCurve(QwtCurve):
    def __init__(self, parent = None, text = None):
        QwtCurve.__init__(self, parent, text)
        self.color = Qt.black

    def draw(self, p, xMap, yMap, f, t):
        p.setBackgroundMode(Qt.OpaqueMode)
        p.setBackgroundColor(self.color)
        p.setBrush(self.color)
##        p.setPen(Qt.black)
        p.setPen(self.color)
        if t < 0: t = self.dataSize() - 1
        if divmod(f, 2)[1] != 0: f -= 1
        if divmod(t, 2)[1] == 0:  t += 1
        for i in range(f, t+1, 2):
            px1 = xMap.transform(self.x(i))
            py1 = yMap.transform(self.y(i))
            px2 = xMap.transform(self.x(i+1))
            py2 = yMap.transform(self.y(i+1))
            p.drawRect(px1, py1, (px2 - px1), (py2 - py1))

class subBarQwtPlotCurve(QwtPlotCurve, subBarQwtCurve): # there must be a better way to do this
    def dummy():
        None

class OWDisplayMotifs(OWWidget):
    settingsList = []
    def __init__(self,parent=None):
        OWWidget.__init__(self, parent, "&Display Motifs", 1)

        #load settings
        self.loadSettings()

        # GUI
        self.box = QVBoxLayout(self.mainArea)
        self.graph = OWGraph(self.mainArea)
        self.graph.setYRlabels(None)
        self.graph.enableGridXB(0)
        self.graph.enableGridYL(0)
        self.graph.setAxisMaxMinor(QwtPlot.xBottom, 10)
        self.graph.setAxisMaxMajor(QwtPlot.xBottom, 10)
        self.graph.setAxisAutoScale(QwtPlot.xBottom)
        self.graph.setAxisScale(QwtPlot.xBottom, 0, 1020, 0)
##        self.graph.setCanvasColor(Qt.white)
        self.box.addWidget(self.graph)

        # inputs
        # data and graph temp variables
        self.inputs = [("Examples", ExampleTable, self.cdata, 0), ("Target", int, self.target), ("Genes", list, self.newGeneList), ("Motifs", list, self.newMotifList)]

        self.data = None
        self.motifLines = []
        self.colorByVariable = None
        self.visibleValues = []
        self.valueToCurveKey = {}
        self.outcomeValues = []
        self.geneList = []
        self.motifList = []

        self.callbackDeposit = []

        self.motifChiSquareThreshold = 0
        self.clusterPostProbThreshold = 0

        # GUI connections
        self.selValues = QVGroupBox(self.space)
        self.selcolorBy = QVGroupBox(self.controlArea)
        self.selcolorBy.setTitle("Color By")
        self.selValues.setTitle("Values")
        self.colorbyQCB = QComboBox(self.selcolorBy)
        self.valuesQLB = QListBox(self.selValues)
        self.unselectAllQLB = QPushButton("Unselect all", self.selValues)

##        self.motifChiSquareThresholdLWS = labelWithSpin_hb(self.selValues, self, "Motif chi-square threshold: ", 0, 100, "mcsth", 1)
        self.motifChiSquareThresholdQHB = QHBox(self.selValues)
        QLabel("Motif chi-square threshold: ", self.motifChiSquareThresholdQHB)
        self.motifChiSquareThresholdQSB = QSpinBox(0, 3000, 100, self.motifChiSquareThresholdQHB)
        self.motifChiSquareThresholdQSB.setValue(2000)
        self.connect(self.motifChiSquareThresholdQSB, SIGNAL("valueChanged(int)"), self.chiSquareChanged)

##        self.clusterThresholdLWS = labelWithSpin_hb(self.selValues, self.postProbChanged, "Cluster PostProb threshold: ", 0, 3000, "clusterThreshold", 100)
        self.clusterThresholdQHB = QHBox(self.selValues)
        QLabel("Cluster PostProb threshold: ", self.clusterThresholdQHB)
        self.clusterThresholdQSB = QSpinBox(0, 100, 1, self.clusterThresholdQHB)
        self.clusterThresholdQSB.setValue(50)
        self.connect(self.clusterThresholdQSB, SIGNAL("valueChanged(int)"), self.clusterPostProbChanged)

        self.valuesQLB.setSelectionMode(QListBox.Multi)
        #connect controls to appropriate functions
        self.connect(self.valuesQLB, SIGNAL("selectionChanged()"), self.valuesSelectionChange)
        self.connect(self.colorbyQCB, SIGNAL('activated (const QString &)'), self.selectColorBy)
        self.connect(self.unselectAllQLB, SIGNAL("clicked()"), self.unselAll)

    def unselAll(self):
        self.valuesQLB.clearSelection()

    def chiSquareChanged(self, value):
        print "cs:", value
        self.motifChiSquareThreshold = value
        self.calcMotifsGraph()
        self.updateMotifsGraph()

    def clusterPostProbChanged(self, value):
        print "cpp:", value
        self.clusterPostProbThreshold = value
        self.calcMotifsGraph()
        self.updateMotifsGraph()

    def valuesSelectionChange(self):
        visibleOutcomes = []
        for i in range(self.valuesQLB.numRows()):
            if self.valuesQLB.isSelected(i):
                visibleOutcomes.append(str(self.valuesQLB.item(i).text()))
        self.visibleValues = visibleOutcomes
        self.updateMotifsGraph()

    def setColorBy(self, list):
        self.colorbyQCB.clear()
        for i in list:
            self.colorbyQCB.insertItem(str(i))
        if len(list) > 0:
            self.colorbyQCB.setCurrentItem(0)

    def selectColorBy(self, variable):
        variable = str(variable)
        self.colorByVariable = variable
        values = self.data.domain[variable].values.native()
        values.sort()
        print values

        ## generate colors for all motifs (each color on its own curve)
        self.graph.removeCurves()
        self.valueToCurveKey = {}
        cn = 0
        allCn = len(values)
        for v in values:
            ## generate curve color
            newColor = QColor()
            newColor.setHsv(cn*255/allCn, 255, 255)
            ## create curve in graph
            curve = subBarQwtPlotCurve(self.graph)
            curve.color = newColor
            ckey = self.graph.insertCurve(curve)
            self.graph.setCurveStyle(ckey, QwtCurve.UserCurve)
            self.valueToCurveKey[v] = ckey
            cn += 1

        self.calcMotifsGraph()
        self.setValuesNames(values)

    def setValuesNames(self, values):
        self.valuesQLB.clear()
        for v in values:
            self.valuesQLB.insertItem(ColorPixmap(self.graph.curve(self.valueToCurveKey[v]).color), v)
        self.valuesQLB.selectAll(TRUE)
##        self.valuesSelectionChange()

    def target(self, targetValue):
        self.targetValue = targetValue

    def newGeneList(self, list):
        self.geneList = list
        self.calcMotifsGraph()
        self.updateMotifsGraph()
 
    def newMotifList(self, list):
        self.motifList = list
        self.calcMotifsGraph()
        self.updateMotifsGraph()

    def cdata(self, data, id):
        self.data = data
        self.motifLines = []

        if self.data == None:
            self.setColorBy( [] )
            self.graph.setYLlabels(None)
        else:
            potentialOutcomes = [v.name for v in self.data.domain if v.varType == orange.VarTypes.Discrete]
            self.setColorBy( potentialOutcomes )
            self.colorByVariable = potentialOutcomes[0]
            lineCn = 0
            self.outcomeValues = self.data.domain[potentialOutcomes[0]].values.native() 
            self.outcomeValues.sort()

            self.motifLines = {}
            for v in self.outcomeValues:
#                print v
                self.motifLines[v] = []

            self.graph.setYLlabels(self.outcomeValues)
            self.graph.setAxisScale(QwtPlot.yLeft, -0.5, len(self.outcomeValues) - 0.5, 1)

            for e in self.data:
                self.motifLines[str(e['sequenceID'].value)].append( e )
##            print self.motifLines.keys()
            self.selectColorBy(self.colorByVariable)

    def calcMotifsGraph(self):
        graphData = {}
        for (colorKey, curveKey) in self.valueToCurveKey.items():
            graphData[colorKey] = [[], []]

        lineCn = 0
        yVals = []
        for lineKey in self.outcomeValues:
            if (self.geneList <> None) and (self.geneList <> []) and (lineKey not in self.geneList): continue
            yVals.append( lineKey)
            for e in self.motifLines[lineKey]:
                motifName = e['motifName'].value
                if (self.motifList <> None) and (self.motifList <> []) and (motifName not in self.motifList): continue

                motifDistToEnd = int(e['distToEnd'].value)
                postProb = int(round(e['pvalue'].value * 100))

                if (1) or (postProb >= self.clusterPostProbThreshold): ## and (chiSquare >= self.motifChiSquareThreshold):
                    colorKey = str(e[self.colorByVariable].value)
                    # print motifNumber, colorKey
                    # set the point x values
                    graphData[colorKey][0].append(motifDistToEnd)
                    graphData[colorKey][0].append(motifDistToEnd + 3.0)
                    # set the point y values
                    graphData[colorKey][1].append(lineCn - 0.45)
                    graphData[colorKey][1].append(lineCn + 0.45)
            lineCn += 1

        self.graph.setYLlabels(yVals)
        self.graph.setAxisScale(QwtPlot.yLeft, -0.5, len(yVals) - 0.5, 1)

        for (colorKey, curveKey) in self.valueToCurveKey.items():
            self.graph.setCurveData(curveKey, graphData[colorKey][0], graphData[colorKey][1])

    def updateMotifsGraph(self):
        for (colorKey, curveKey) in self.valueToCurveKey.items():
            self.graph.curve(curveKey).setEnabled(colorKey in self.visibleValues)
        self.graph.update()

if __name__ == "__main__":
    a = QApplication(sys.argv)
    owdm = OWDisplayMotifs()
    a.setMainWidget(owdm)
    owdm.show()
    a.exec_loop()
    owdm.saveSettings()
