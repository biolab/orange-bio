"""
<name>Display Motifs</name>
<description>None.</description>
<contact>Tomaz Curk</contact>
<icon>icons\GenomeMap.png</icon>
<prototype>1</prototype>
"""

from collections import defaultdict

import orange
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWGraph import *
from Orange.OrangeWidgets.OWTools import *
from Orange.OrangeWidgets.OWWidget import *

class subBarQwtPlotCurve(QwtPlotCurve):
    def __init__(self, text = None):
        QwtPlotCurve.__init__(self, text or "")
        self.color = Qt.black

    def draw(self, p, xMap, yMap, f, t=-1):
        p.setBackgroundMode(Qt.OpaqueMode)
        p.setBrush(self.color)
        p.setPen(self.color)
        
        if t < 0:
            t = self.dataSize() - 1
            f = 0
        if divmod(f, 2)[1] != 0:
            f -= 1
        if divmod(t, 2)[1] == 0:
            t += 1
            
        for i in range(f, t+1, 2):
            px1 = xMap.transform(self.x(i))
            py1 = yMap.transform(self.y(i))
            px2 = xMap.transform(self.x(i+1))
            py2 = yMap.transform(self.y(i+1))
            p.drawRect(px1, py1, (px2 - px1), (py2 - py1))


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
        self.graph = OWGraph(self.mainArea)
        self.graph.setYRlabels(None)
        self.graph.enableGridXB(0)
        self.graph.enableGridYL(1)
        self.graph.setAxisMaxMinor(QwtPlot.xBottom, 10)
        self.graph.setAxisMaxMajor(QwtPlot.xBottom, 10)
        self.graph.setAxisAutoScale(QwtPlot.xBottom)
        self.graph.setAxisScale(QwtPlot.xBottom, -1020, 0, 0)
        self.mainArea.layout().addWidget(self.graph)
        
        # inputs
        # data and graph temp variables
        self.inputs = [("Examples", ExampleTable, self.cdata, Default), ("Genes", list, self.newGeneList, Default), ("Motifs", list, self.newMotifList, Default)]

        self.data = None
        self.motifLines = []
        self.visibleValues = []
        self.valueToCurve = {}
        self.allGenes = [] ## genes displayed always in same order
        self.geneList = [] ## selected genes
        self.motifList = [] ## selected motifs
        self.valuesPresentInData = []

        self.clusterPostProbThreshold = 0

        # GUI
        self.selValues = OWGUI.widgetBox(self.controlArea, "Values")
        self.selcolorBy = OWGUI.widgetBox(self.controlArea, "Color By")

        self.colorByCombo = OWGUI.comboBox(self.selcolorBy, self, "colorBy",
                                           items=[],
                                           callback=self.colorByChanged)
        
        self.pvalThresholdCombo = OWGUI.comboBox(self.selValues, self, "pvalThresholdIndex",
                                                 items=[],
                                                 callback=self.pvalThresholdChanged)
        
        self.valuesQLB = QListWidget(self.selValues)
        self.valuesQLB.setSelectionMode(QListWidget.MultiSelection)
        self.connect(self.valuesQLB, SIGNAL("itemSelectionChanged()"), self.valuesSelectionChange)
        self.selValues.layout().addWidget(self.valuesQLB)
        
        self.unselectAllQLB = OWGUI.button(self.selValues, self, "Unselect all",
                                           callback = self.unselAll)

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
        for i in range(self.valuesQLB.count()):
            item = self.valuesQLB.item(i)
            if item.isSelected():
                visibleOutcomes.append(str(item.text()))
        self.visibleValues = visibleOutcomes
        self.updateMotifsGraph()

    def updateSelectionChange(self):
        ## make new colors only for those colorBy values present in data
        colors = ColorPaletteHSV(len(self.valuesPresentInData))
        for (i, v) in enumerate(self.valuesPresentInData):
            curve = self.valueToCurve[v]
            curve.color = colors[i]

        self.setValuesNames(self.valuesPresentInData) 
        for i in range(self.valuesQLB.count()):
            item = self.valuesQLB.item(i)
            item.setSelected(str(item.text()) in self.visibleValues)
#            else:
#                self.valuesQLB.setSelected(i, 0)

    def colorByChanged(self):
#        self.graph.removeCurves()
        self.graph.removeDrawingCurves()
        self.valueToCurve = {}

        if self.colorBy == None:
            self.setValuesNames(None)
            return
            
        if self.potentialColorVariables:
            var = self.potentialColorVariables[self.colorBy]
            values = sorted(var.values)
#        values.sort()

        ## generate colors for each colorby value (each color on its own curve)
        colors = ColorPaletteHSV(len(values))
        for (i, v) in enumerate(values):
            ## create curve in graph
            curve = subBarQwtPlotCurve()
            curve.color = colors[i]
            curve.attach(self.graph)
#            ckey = self.graph.insertCurve(curve)
            curve.setStyle(QwtPlotCurve.UserCurve)
            self.valueToCurve[v] = curve

        self.setValuesNames(values)
        self.calcMotifsGraph()

    def setValuesNames(self, values):
        self.valuesQLB.clear()
        if values == None:
            return
        for v in values:
            self.valuesQLB.addItem(QListWidgetItem(QIcon(ColorPixmap(self.valueToCurve[v].color)), v))
            self.valuesQLB.item(self.valuesQLB.count() - 1).setSelected(True)
        

    ## signal processing
    def newGeneList(self, list):
        self.geneList = []
        
        if list is not None:
            self.geneList = list
            self.calcMotifsGraph()
            self.updateMotifsGraph()
 
    def newMotifList(self, list):
        self.motifList = []
        
        if list is not None:
            self.motifList = list
            self.calcMotifsGraph()
            self.updateMotifsGraph()
        
    def cdata(self, data):
        self.data = data
        self.motifLines = []
        self.colorByCombo.clear()
    
        if self.data == None:
            self.colorBy = None
            self.potentialColorVariables = []
            self.colorByChanged()
            self.graph.setYLlabels(None)
        else:
            self.potentialColorVariables = [v for v in self.data.domain.variables + self.data.domain.getmetas().values() \
                                            if v.varType == orange.VarTypes.Discrete]
            for var in self.potentialColorVariables:
                self.colorByCombo.addItem(str(var.name))
    
            ## break motif info according to gene (sequenceID)
            self.motifLines = defaultdict(list)
            for e in self.data:
                geneID = str(e['sequenceID'])
                self.motifLines[geneID].append(e)
                
            self.allGenes = sorted(self.motifLines.keys())
            if self.potentialColorVariables:
                self.colorBy = min(len(self.potentialColorVariables) - 1, self.colorBy or 0)
            else:
                self.colorBy = None
                
            self.colorByChanged()
    
        if len(self.potentialColorVariables) > 0:
            self.colorByCombo.setCurrentIndex(0)

    ## update graph
    def calcMotifsGraph(self):
        graphData = {}
        for (colorKey, curve) in self.valueToCurve.iteritems():
            graphData[colorKey] = [[], []]

        lineCn = 0
        yVals = []
        self.colorByValuesPresentInData = []
        self.valuesPresentInData = []
        for lineKey in self.geneList:
##            if (self.geneList <> None) and (self.geneList <> []) and (lineKey not in self.geneList): continue
            yVals.append( lineKey)
            for e in self.motifLines[lineKey]:
                motifName = str(e['motifName'])#.value
                if self.motifList is not None and self.motifList != [] and motifName not in self.motifList:
                    continue

                motifDistToEnd = -int(e['distToEnd']) #.value)
                postProb = int(round(float(e['pvalue']) * 100))
                dir = str(e['direction'])

                if 1 or (postProb >= self.clusterPostProbThreshold): ## and (chiSquare >= self.motifChiSquareThreshold):
                    colorAttr = self.potentialColorVariables[self.colorBy]
                    colorKey = str(e[colorAttr])
                    if colorKey not in self.valuesPresentInData:
                        self.valuesPresentInData.append(colorKey)
                    # print motifNumber, colorKey
                    # set the point x values
                    graphData[colorKey][0].append(motifDistToEnd)
                    graphData[colorKey][0].append(motifDistToEnd + 5.0)
                    # set the point y values
                    if dir == "1":
                        graphData[colorKey][1].append(lineCn + 0.0)
                        graphData[colorKey][1].append(lineCn + 0.30)
                    elif dir == "-1":
                        graphData[colorKey][1].append(lineCn - 0.30)
                        graphData[colorKey][1].append(lineCn + 0.0)
                    else:
                        graphData[colorKey][1].append(lineCn - 0.30)
                        graphData[colorKey][1].append(lineCn + 0.30)
            lineCn += 1

        vals = reduce(list.__add__, [val[0] for val in graphData.itervalues()], [])
        minX = min(vals or [0]) - 5.0
        maxX = max(vals or [0]) + 5.0
        
        self.graph.setYLlabels(yVals)
        self.graph.setAxisScale(QwtPlot.yLeft, -0.5, len(yVals) - 0.5, 1)
        self.graph.setAxisScale(QwtPlot.xBottom, minX, maxX)

        for (colorKey, curve) in self.valueToCurve.iteritems():
            curve.setData(graphData[colorKey][0], graphData[colorKey][1])
            curve.setVisible(colorKey in self.visibleValues)
        self.graph.replot()
        self.valuesPresentInData.sort()
        self.visibleValues = [v for v in self.valuesPresentInData]
        self.updateSelectionChange()
        
    def updateMotifsGraph(self):
        for (colorKey, curve) in self.valueToCurve.items():
            curve.setVisible(colorKey in self.visibleValues)
        self.graph.replot()

if __name__ == "__main__":
    a = QApplication(sys.argv)
    owdm = OWDisplayMotifs()
    owdm.cdata(orange.ExampleTable(os.path.expanduser("~/Desktop/motif.tab")))
    owdm.newGeneList(["1", "2"])
    owdm.show()
    a.exec_()
    owdm.saveSettings()
