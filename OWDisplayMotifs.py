"""
<name>Display Motifs</name>
<description>None.</description>
<icon>icons\DisplayMotifs.png</icon>
"""

import orange
from OWTools import *
from OWWidget import *
from OWGraph import *
import OWGUI

from OWGraph import ColorPaletteHSV

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
    def __init__(self,parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, "&Display Motifs", 0)

        # set default settings
        self.colorBy = None
        self.pvalThresholdIndex = None
        self.pvalThreshold = None

        #load settings
        self.loadSettings()

        # GUI
        self.box = QVBoxLayout(self.mainArea)
        self.graph = OWGraph(self.mainArea)
        self.graph.setYRlabels(None)
        self.graph.enableGridXB(0)
        self.graph.enableGridYL(1)
        self.graph.setAxisMaxMinor(QwtPlot.xBottom, 10)
        self.graph.setAxisMaxMajor(QwtPlot.xBottom, 10)
        self.graph.setAxisAutoScale(QwtPlot.xBottom)
        self.graph.setAxisScale(QwtPlot.xBottom, -1020, 0, 0)
        self.box.addWidget(self.graph)

        # inputs
        # data and graph temp variables
        self.inputs = [("Examples", ExampleTable, self.cdata, 0), ("Genes", list, self.newGeneList), ("Motifs", list, self.newMotifList)]

        self.data = None
        self.motifLines = []
        self.visibleValues = []
        self.valueToCurveKey = {}
        self.allGenes = [] ## genes displayed always in same order
        self.geneList = [] ## selected genes
        self.motifList = [] ## selected motifs
        self.valuesPresentInData = []

        self.clusterPostProbThreshold = 0

        # GUI connections
        self.selValues = QVGroupBox(self.space)
        self.selcolorBy = QVGroupBox(self.controlArea)
        self.selcolorBy.setTitle("Color By")
        self.selValues.setTitle("Values")

        self.colorByCombo = OWGUI.comboBox(self.selcolorBy, self, "colorBy", items=[], callback=self.colorByChanged)
        self.pvalThresholdCombo = OWGUI.comboBox(self.selValues, self, "pvalThresholdIndex", items=[], callback=self.pvalThresholdChanged)
        
        self.valuesQLB = QListBox(self.selValues)
        self.unselectAllQLB = QPushButton("Unselect all", self.selValues)
        self.valuesQLB.setSelectionMode(QListBox.Multi)
        self.connect(self.valuesQLB, SIGNAL("selectionChanged()"), self.valuesSelectionChange)
        self.connect(self.unselectAllQLB, SIGNAL("clicked()"), self.unselAll)

    def unselAll(self):
        self.valuesQLB.clearSelection()

    def pvalThresholdChanged(self):
        print self.pvalThresholdIndex
        if self.pvalThresholdIndex:
            self.pvalThreshold = float(self.pvalThresholdCombo.text(self.pvalThresholdIndex))
        else:
            self.pvalThreshold = None
        self.calcMotifsGraph()
        self.updateMotifsGraph()

    def valuesSelectionChange(self):
        visibleOutcomes = []
        for i in range(self.valuesQLB.numRows()):
            if self.valuesQLB.isSelected(i):
                visibleOutcomes.append(str(self.valuesQLB.item(i).text()))
        self.visibleValues = visibleOutcomes
        self.updateMotifsGraph()

    def updateSelectionChange(self):
        ## make new colors only for those colorBy values present in data
        colors = ColorPaletteHSV(len(self.valuesPresentInData), 255)
        for (i, v) in enumerate(self.valuesPresentInData):
            ckey = self.valueToCurveKey[v]
            self.graph.curve(ckey).color = colors[i]

        self.setValuesNames(self.valuesPresentInData) 
        for i in range(self.valuesQLB.numRows()):
            if str(self.valuesQLB.item(i).text()) in self.visibleValues:
                self.valuesQLB.setSelected(i, 1)
            else:
                self.valuesQLB.setSelected(i, 0)

    def colorByChanged(self):
        self.graph.removeCurves()
        self.valueToCurveKey = {}

        if self.colorBy == None:
            self.setValuesNames(None)
            return
            
        values = self.data.domain[self.colorBy].values.native()
        values.sort()

        ## generate colors for each colorby value (each color on its own curve)
        colors = ColorPaletteHSV(len(values), 255)
        for (i, v) in enumerate(values):
            ## create curve in graph
            curve = subBarQwtPlotCurve(self.graph)
            curve.color = colors[i]
            ckey = self.graph.insertCurve(curve)
            self.graph.setCurveStyle(ckey, QwtCurve.UserCurve)
            self.valueToCurveKey[v] = ckey

        self.setValuesNames(values)
        self.calcMotifsGraph()

    def setValuesNames(self, values):
        self.valuesQLB.clear()
        if values == None:
            return
        for v in values:
            self.valuesQLB.insertItem(ColorPixmap(self.graph.curve(self.valueToCurveKey[v]).color), v)
        self.valuesQLB.selectAll(TRUE)

    ## signal processing
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
        self.colorByCombo.clear()

        if self.data == None:
            self.colorByChanged()
            self.graph.setYLlabels(None)
            potentialColorVariables = []
        else:
            potentialColorVariables = [v.name for v in self.data.domain if v.varType == orange.VarTypes.Discrete]
            for i in potentialColorVariables:
                self.colorByCombo.insertItem(str(i))

            ## break motif info according to gene (sequenceID)
            self.motifLines = {}
            for e in self.data:
                geneID = str(e['sequenceID'].value)
                tmpl = self.motifLines.get(geneID, [])
                tmpl.append( e)
                self.motifLines[geneID] = tmpl
            self.allGenes = self.motifLines.keys()
            self.allGenes.sort() ## genes displayed always in same order

            self.colorBy = 1;
            self.colorByChanged()

        if len(potentialColorVariables) > 0:
            self.colorByCombo.setCurrentItem(0)

    ## update graph
    def calcMotifsGraph(self):
        graphData = {}
        for (colorKey, curveKey) in self.valueToCurveKey.items():
            graphData[colorKey] = [[], []]

        lineCn = 0
        yVals = []
        self.colorByValuesPresentInData = []
        self.valuesPresentInData = []
        for lineKey in self.geneList:
##            if (self.geneList <> None) and (self.geneList <> []) and (lineKey not in self.geneList): continue
            yVals.append( lineKey)
            for e in self.motifLines[lineKey]:
                motifName = e['motifName'].value
                if (self.motifList <> None) and (self.motifList <> []) and (motifName not in self.motifList): continue

                motifDistToEnd = -int(e['distToEnd'].value)
                postProb = int(round(e['pvalue'].value * 100))
                dir = str(e['direction'])

                if (1) or (postProb >= self.clusterPostProbThreshold): ## and (chiSquare >= self.motifChiSquareThreshold):
                    colorKey = str(e[self.colorBy].value)
                    if colorKey not in self.valuesPresentInData:
                        self.valuesPresentInData.append( colorKey)
                    # print motifNumber, colorKey
                    # set the point x values
                    graphData[colorKey][0].append(motifDistToEnd)
                    graphData[colorKey][0].append(motifDistToEnd + 3.0)
                    # set the point y values
                    if dir == "1":
                        graphData[colorKey][1].append(lineCn - 0.45)
                        graphData[colorKey][1].append(lineCn + 0.0)
                    elif dir == "-1":
                        graphData[colorKey][1].append(lineCn + 0.0)
                        graphData[colorKey][1].append(lineCn + 0.45)
                    else:
                        graphData[colorKey][1].append(lineCn - 0.45)
                        graphData[colorKey][1].append(lineCn + 0.45)
            lineCn += 1

        self.graph.setYLlabels(yVals)
        self.graph.setAxisScale(QwtPlot.yLeft, -0.5, len(yVals) - 0.5, 1)

        for (colorKey, curveKey) in self.valueToCurveKey.items():
            self.graph.curve(curveKey).setEnabled(FALSE)
            self.graph.setCurveData(curveKey, graphData[colorKey][0], graphData[colorKey][1])
            self.graph.curve(curveKey).setEnabled(colorKey in self.visibleValues)

        self.valuesPresentInData.sort()
        self.visibleValues = [v for v in self.valuesPresentInData]
        self.updateSelectionChange()
        
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
