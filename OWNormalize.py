"""
<name>Normalize Microarray Data</name>
<description>Normalization of custom cDNA microarray data.</description>
<icon>icons/Unknown.png</icon>
<priority>1150</priority>
<author>Peter Juvan (peter.juvan@fri.uni-lj.si)</author>
"""

import string, math
import Numeric, MA
import orange
from qttable import *
from OWWidget import *
import OWGUI
from OWGraph import *
from OWGraphTools import *      # color palletes, user defined curves, ...
import ColorPalette             # ColorButton


class Probe:
    def __init__(self, ID, ratio, color=QColor(0,0,255), symbol=1):
        """Probe name and its expected ratio;
        curve represents normalization factors where color and symbol is used to plot the curve on the graph.
        """
        self.ID = ID                # string
        self.ratio = ratio          # float
        self.color = color          # QColor
        self.symbol = symbol        # int (QwtSymbol.Style)
        self._signalSmpl = None     # Numeric.arrays
        self._signalRef = None
        self._bgSmpl = None
        self._bgRef = None
        self.curve = None


    def setData(self, data, varID, varSignalSmpl, varSignalRef, varBGSmpl, varBGRef):
        """Input: orange.ExampleTable, orange.Variable representing ID, signal/background, sample/reference,
        """
##        print "setData: varID, varSignalSmpl, varSignalRef, varBGSmpl, varBGRef", map(lambda x: x.name, [varID, varSignalSmpl, varSignalRef, varBGSmpl, varBGRef])
        pp = orange.Preprocessor_take()
        pp.values[varID] = self.ID
        etSub = pp(data)
        if len(etSub) == 0:
            self._signalSmpl = None
            self._signalRef = None
            self._bgSmpl = None
            self._bgRef = None
        else:
            ma = etSub.toMA("a")[0]
            self._signalSmpl = ma[:,data.domain.index(varSignalSmpl)]
            self._signalRef = ma[:,data.domain.index(varSignalRef)]
            self._bgSmpl = ma[:,data.domain.index(varBGSmpl)]
            self._bgRef = ma[:,data.domain.index(varBGRef)]


    def getSymbolPixmap(self, rect):
        """Input: QRect instance; output: QPixmap instance.
        """
        # init pixmap for table item
        symbolPixmap = QPixmap(rect.width(),rect.height())
        symbolPixmap.fill(QColor(255,255,255))
        painter = QPainter(symbolPixmap)
        symbol = QwtSymbol(self.symbol, QBrush(self.color, QBrush.SolidPattern), QPen(QColor(0,0,0),1), QSize(8,8))
        symbol.draw(painter, QPoint(rect.width()/2,rect.height()/2))
        painter.end()
        return symbolPixmap


    def getA_NF(self, subtrBG):
        print "TODO: use filters"
##            if filterName:
##                pp.values[n2v[filterName]] = "accept"
        if self._signalSmpl and self._signalRef and self._bgSmpl and self._bgRef:
            if subtrBG:
                netSmpl = self._signalSmpl - self._bgSmpl
                netRef =  self._signalRef  - self._bgRef
            else:
                netSmpl = self._signalSmpl
                netRef =  self._signalRef
            A = MA.sqrt(netSmpl*netRef)
            NF = self.ratio / (netSmpl/netRef)
            condition = Numeric.logical_not(Numeric.logical_or(MA.getmaskarray(A), MA.less_equal(A, 0).filled(1)))
##            print "condition:\n", condition
##            print "A:\n", A
##            print "NF:\n", NF
            return Numeric.asarray(MA.compress(condition,A)), Numeric.asarray(MA.compress(condition,NF))
        else:
            return Numeric.zeros((0,), Numeric.Float), Numeric.zeros((0,), Numeric.Float)


    def removeCurve(self, graph, replot):
        if self.curve:
            graph.removeCurve(self.curve)
            if replot:
                graph.replot()


    def replotCurve(self, graph, subtrBG, logAxisX, logAxisY, markerSize, replotProbeCurves, filters=None):
        print "TODO: use filters"
##            if filterName:
##                pp.values[n2v[filterName]] = "accept"
        if self.curve:
            graph.removeCurve(self.curve)
        self.curve = graph.insertCurve(self.ID)
        A, NF = self.getA_NF(subtrBG)
        if logAxisX:
            A = Numeric.log(A) / math.log(math.e)
        if logAxisY:
            NF = Numeric.log(NF) / math.log(math.e)
        graph.setCurveData(self.curve, A, NF)
        graph.setCurveStyle(self.curve, QwtCurve.NoCurve)
        qSymbol = QwtSymbol(self.symbol, QBrush(self.color, QBrush.SolidPattern), QPen(QColor(0,0,0),1), QSize(markerSize,markerSize))
        graph.setCurveSymbol(self.curve, qSymbol)
        if replotProbeCurves:
            graph.replot()


    def setCurveActive(self, graph, active):
        oldSymbol = graph.curveSymbol(self.curve)
        pen = oldSymbol.pen()
        if (active and pen.width() != 2) or (not active and pen.width() != 1):
            if active:
                pen.setWidth(2)
            else:
                pen.setWidth(1)
            newSymbol = QwtSymbol(oldSymbol.style(), oldSymbol.brush(), pen, oldSymbol.size())
            graph.setCurveSymbol(self.curve, newSymbol)
            graph.replot()
        
        


class OWNormalize(OWWidget):

    settingsList = [""]

    def __init__(self, parent = None, signalManager = None, name = "Normalize Microarray Data"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, wantStatusBar=False)  #initialize base class

        # set channels
        self.inputs = [("Examples", ExampleTable, self.onDataInput, Default), ("Probes", ExampleTable, self.onProbesInput)]
        self.outputs = [("Examples", ExampleTable), ("Probes", ExampleTable)]

        # constants
        self.stylesIndices = zip(["Circle","Rect","Diamond","Triangle","DTriangle","UTriangle","LTriangle","RTriangle","Cross","XCross"],
                                 range(1,10))
        self.sizeButtonColor = 18
        self.tblControlsColumn2Width = 43

        # general settings
        self.controlName = ""
        self.controlRatio = "1.0"
        self._probeSymbol = 0    # index for cmbProbeSymbol, equals to index+1 for QwtSymbol.Style
        self.probeColor = QColor(0,0,255)
        self.subtrBG = True
        self.useCV = True
        self.CV = "0.5"
        self.useMinIntensity = True
        self.minIntensity = "1.5"
        self.useMaxIntensity = True
        self.maxIntensity = "50000"
        self.normType = 0  #0: median, 1: LR, 2: LOESS
        # graph
        self.logAxisX = True
        self.logAxisY = False
        self.markerSize = 8
        self.mergeReplGraph = False
        # output
        self.mergeRepl = True
        self.mergeReplType = 0     #0: mean, 1: median
        self.mergeOtherType = 0     #0: use first, 1: concatenate
        self.outNumNonFiltered = True
        self.outNetSignal = True
        self.outNonNormLogRatio = True
        self.loadSettings()
        
        # context-specific settings
        self.data = None
        self.dataProbes = None
        self.varsFloat = []
        self.varsEnum = []
        self.varsOther = []
        # variable indices selected within combos
        self.varIdxID = None            # index in self.varsEnum
        self.varIdxSignalSmpl = None    # index in self.varsFloat
        self.varIdxSignalRef = None     # index in self.varsFloat
        self.varIdxBGSmpl = None        # index in self.varsFloat
        self.varIdxBGRef = None         # index in self.varsFloat
        # names of selected vars from listBox
        self.varOthersSelectedNames = []
        # controls/ratios
        self.inputContRatios = {}        # dict {controlName:ratio} filled from external exampleTable
        self.internContRatios = {}       # dict {controlName:ratio} updated together with self.tblControls
        self.probesExternal = {}        # dict {Probe.name: Probe} from external data
        self.probesUsed = {}            # dict {Probe.name: Probe} updated together with self.tblControls
        self.probeActive = None         # currently active probe

        # GUI
        self.controlArea.setFixedWidth(240)
        self.resize(800, 700)
        # control area: tabs
        self.tabs = QTabWidget(self.controlArea, 'tabWidget')
        # tab 1: vars
        boxVars = QVGroupBox(self)
        self.tabs.insertTab(boxVars, "Variables")
        self.cmbVarID = OWGUI.comboBox(boxVars, self, "varIdxID", "Probe ID", callback=self.varIDChange)
        self.cmbVarSignalSmpl = OWGUI.comboBox(boxVars, self, "varIdxSignalSmpl", "Sample signal", callback=lambda x="varSignalSmpl": self.varChange(x))
        self.cmbVarSignalRef = OWGUI.comboBox(boxVars, self, "varIdxSignalRef", "Reference signal", callback=lambda x="varSignalRef": self.varChange(x))
        self.cmbVarBGSmpl = OWGUI.comboBox(boxVars, self, "varIdxBGSmpl", "Sample background", callback=lambda x="varBGSmpl": self.varChange(x))
        self.cmbVarBGRef = OWGUI.comboBox(boxVars, self, "varIdxBGRef", "Reference background", callback=lambda x="varBGRef": self.varChange(x))
        OWGUI.button(boxVars, self, "Search for Default Names", callback=self.searchDefaultVarAssignmentClick)
        boxOtherVars = QVGroupBox(boxVars)
        boxOtherVars.setTitle('Other variables')
        self.lbVarOthers = QListBox(boxOtherVars)
        self.lbVarOthers.setSelectionMode(QListBox.Multi)
        self.connect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)

        # tab 2: table probe/ratio/marker
        boxControls = QVGroupBox(boxVars)
        self.tabs.insertTab(boxControls, "Controls")
        self.tblControls = QTable(boxControls)
        self.tblControls.setSelectionMode(QTable.Multi)
        self.tblControls.setNumCols(3)
        self.tblControls.setColumnWidth(2,self.tblControlsColumn2Width)
        self.connect(self.tblControls, SIGNAL("valueChanged(int,int)"), self.tblControlsChange)
        self.connect(self.tblControls, SIGNAL("currentChanged(int, int)"), self.tblControlsCurrentChanged)
        hheader=self.tblControls.horizontalHeader()
        hheader.setLabel(0, "Probe ID")
        hheader.setLabel(1, "Ratio")
        hheader.setLabel(2, "Marker")
        hheader.setClickEnabled(False)
        hheader.setMovingEnabled(False)
        hheader.setResizeEnabled(False)
        vheader=self.tblControls.verticalHeader()
        vheader.setClickEnabled(False)
        vheader.setMovingEnabled(False)
        vheader.setResizeEnabled(False)

        # tab 2: buttons
        boxBtns0 = QVGroupBox(boxControls)
        boxBtns0.setTitle("Select probes where ID contains")
        boxBtns00 = QHBox(boxBtns0)
        OWGUI.lineEdit(boxBtns00, self, "controlName")
        OWGUI.button(boxBtns00, self, "Select", callback=self.selectControlsClick)
        boxBtns01 = QHBox(boxBtns0)
        OWGUI.button(boxBtns01, self, "Select all", callback=self.selectControlsAllClick)
        OWGUI.button(boxBtns01, self, "Unselect all", callback=lambda repaint=True: self.tblControls.clearSelection(repaint))
        # empty pixmap to fill self.tblControls column 2
        self.tblControls.setNumRows(1)
        rect = self.tblControls.cellGeometry(0,2)
        self.symbolPixmapEmpty = QPixmap(rect.width(),rect.height())
        self.symbolPixmapEmpty.fill(QColor(255,255,255))
        self.tblControls.setNumRows(0)

        boxBtns1 = QVGroupBox(boxControls)
        boxBtns1.setTitle("Set ratio and marker for selected probes")
        boxBtns11 = QHBox(boxBtns1)
        OWGUI.lineEdit(boxBtns11, self, "controlRatio")
        self.btnProbeColor = OWGUI.button(boxBtns11, self, "", callback=self.probeColorClick)
        pxm = QPixmap(self.sizeButtonColor,self.sizeButtonColor)
        pxm.fill(self.probeColor)
        self.btnProbeColor.setPixmap(pxm)
        self.cmbProbeSymbol = OWGUI.comboBox(boxBtns11, self, "_probeSymbol")
        for styleName, styleIdx in self.stylesIndices:
            symbol = QwtSymbol(styleIdx, QBrush(QColor(255,255,255), QBrush.SolidPattern), QPen(QColor(0,0,0),1), QSize(self.markerSize,self.markerSize))
            pixmap = QPixmap(14,14)
            pixmap.fill(QColor(255,255,255))
            painter = QPainter(pixmap)
            symbol.draw(painter, QPoint(7,7))
            painter.end()
            self.cmbProbeSymbol.insertItem(pixmap, styleName)

        boxBtns12 = QHBox(boxBtns1)
        OWGUI.button(boxBtns12, self, "Set", callback=self.setRatiosClick)
        OWGUI.button(boxBtns12, self, "Clear", callback=self.clearRatiosClick)

        # tab 3: settings
        boxSettings = QVGroupBox(self)
        self.tabs.insertTab(boxSettings, "Settings")
        # tab 3: settings: filters
        boxFilters = QVGroupBox(boxSettings)
        boxFilters.setTitle('Filters')
        self.cbsubtrBG = OWGUI.checkBox(boxFilters, self, "subtrBG", "Subtract background", callback=self.settingsFiltersChange)
        # tab 3: settings: filters: CV
        boxCV = QHBox(boxFilters)
        self.cbCV = OWGUI.checkBox(boxCV, self, "useCV", "Coefficient of variation (CV):", callback=self.settingsFiltersChange)
        self.leCV = OWGUI.lineEdit(boxCV, self, "CV", callback=self.settingsFiltersChange)
        # tab 3: settings: filters: min.intensity
        boxMinIntensity = QHBox(boxFilters)
        self.cbMinIntensity = OWGUI.checkBox(boxMinIntensity, self, "useMinIntensity", "Min. intensity:", callback=self.settingsFiltersChange)
        self.leMinIntensity = OWGUI.lineEdit(boxMinIntensity, self, "minIntensity", callback=self.settingsFiltersChange)
        QLabel(" times above bg.", boxMinIntensity)
        # tab 3: settings: filters: max.intensity
        boxMaxIntensity = QHBox(boxFilters)
        self.cbMaxIntensity = OWGUI.checkBox(boxMaxIntensity, self, "useMaxIntensity", "Max. intensity:", callback=self.settingsFiltersChange)
        self.leMaxIntensity = OWGUI.lineEdit(boxMaxIntensity, self, "maxIntensity", callback=self.settingsFiltersChange)
        # tab 3: settings: normalization type
        rbgNormalization = OWGUI.radioButtonsInBox(boxSettings, self, value="normType", btnLabels=["Median (intensity independent)", "Linear regression", "Loess"], box="Normalization Type", callback=self.settingsNormTypeChange)
        # tab 3: settings: graph
        boxGraph = QVGroupBox(boxSettings)
        boxGraph.setTitle('Graph')
        self.cbLogAxisX = OWGUI.checkBox(boxGraph, self, "logAxisX", "Logarithmic X axis", callback=lambda ax=2: self.settingsGraphAxisChange(ax))
        self.cbLogAxisY = OWGUI.checkBox(boxGraph, self, "logAxisY", "Logarithmic Y axis", callback=lambda ax=0: self.settingsGraphAxisChange(ax))
        boxMSize = QHBox(boxGraph)
        QLabel("Marker size", boxMSize)
        self.cmbMarkerSize = OWGUI.comboBox(boxMSize, self, "markerSize", callback=self.settingsGraphChange, sendSelectedValue=1, valueType=int)
        for itemIdx, size in enumerate(range(1,16)):
            self.cmbMarkerSize.insertItem(str(size))
            if self.markerSize == size:
                self.cmbMarkerSize.setCurrentItem(itemIdx)
        self.cbMergeReplicasGraph = OWGUI.checkBox(boxGraph, self, value="mergeReplGraph", label="Merge replicas", callback=self.settingsGraphChange)
        # tab 3: settings: replication
        boxOutput = QVGroupBox(boxSettings)
        boxOutput.setTitle('Output')
##        boxRepl = QVButtonGroup("Replication", boxSettings)
        self.cbMergeReplicas = OWGUI.checkBox(boxOutput, self, value="mergeRepl", label="Merge replicas", callback=self.settingsReplicasChange)
        self.rbgMergeReplType = OWGUI.radioButtonsInBox(boxOutput, self, value="mergeReplType", btnLabels=["Mean", "Median"], box="Merge intensities", callback=self.settingsReplicasChange)
        self.rbgMergeOtherType = OWGUI.radioButtonsInBox(boxOutput, self, value="mergeOtherType", btnLabels=["Use first value", "Concatenate values"], box="Merge other variables", callback=self.settingsReplicasChange)
        # tab 3: settings: other output
        self.cbOutNumNonFiltered = OWGUI.checkBox(boxOutput, self, "outNumNonFiltered", "Number of replicas", callback=self.settingsOutputChange)
        self.cbOutNetSignal = OWGUI.checkBox(boxOutput, self, "outNetSignal", "Net intensity", callback=self.settingsOutputChange)
        self.cbOutNonNormLogRatio = OWGUI.checkBox(boxOutput, self, "outNonNormLogRatio", "Non-normalized log ratio", callback=self.settingsOutputChange)

        # main area: graph
        boxG = QVBox(self.mainArea)
        graphBoxLayout = QVBoxLayout(self.mainArea)
        graphBoxLayout.addWidget(boxG)
        self.graph = OWGraph(boxG)
        self.graph.setAutoLegend(False)
        self.graph.setAutoReplot(False)
##        self.graph.setAxisAutoScale(False)
##        self.graph.setAxisScale(0,0,7)
##        self.graph.setAxisScale(2,0,9)
        self.setGraphAxes(axes=[0,2])
        self.connect(self.graphButton, SIGNAL("clicked()"), self.graph.saveToFile)


    ###################################################################################
    ## EVENT HANDLERS: DATA IN / OUT

    def onDataInput(self, data):
        """Handles input of new data.
        """
        print "onDataInput"
        self.varsFloat = []
        self.varsEnum = []
        self.varsOther = []
        self.data = None
        if data <> None:
            if len(data.domain.getmetas())>0:
                # domain with all variables + metas as attributes
                domNoMeta = orange.Domain(data.domain.variables + data.domain.getmetas().values(), None)
                self.data = orange.ExampleTable(domNoMeta, data)
            else:
                self.data = data
            # divide vars to Enum, Float and String
            for var in self.data.domain.variables:
                if var.varType == orange.VarTypes.Continuous:
                    self.varsFloat.append(var)
                elif var.varType == orange.VarTypes.Discrete:
                    self.varsEnum.append(var)
                else:
                    self.varsOther.append(var)
        # fill combos and listBox with variables
        self.fillCmbVars()
        self.searchDefaultVarAssignmentClick()
        self.fillLbVarOthers()
        if self.data and self.varIdxID != None and self.dataProbes != None:
            self.processDataProbes(self.dataProbes, self.varsEnum[self.varIdxID].name)
        self.initProbes()
        self.sendData()
        self.sendProbes()


    def onProbesInput(self, dataProbes):
        """Handles input of probes data; appends Probe to self.probesExternal dict.
        """
        print "onProbesInput"
        self.dataProbes = dataProbes
        if self.data and self.varIdxID != None and dataProbes != None:
            self.processDataProbes(dataProbes, self.varsEnum[self.varIdxID].name)
            self.initProbes()
            self.sendData()
            self.sendProbes()


    def processDataProbes(self, dataProbes, varIDName):
        """Adds probes to self.probesExternal dict.
        """
        print "processDataProbes"
        varNames =  map(lambda var: var.name, dataProbes.domain.variables)
        if varIDName in varNames and "Ratio" in varNames:
            for ex in dataProbes:
                probe = Probe(ex[varIDName].native(), ex["Ratio"].native())
                try:
                    c = str(ex["ColorHSV"])
                    probe.color = QColor(int(c[0:2],16), int(c[2:4],16), int(c[4:6],16), QColor.Hsv)
                except TypeError:
                    pass
                try:
                    probe.symbol = int(ex["Symbol"].native())
                except TypeError:
                    pass
                self.probesExternal[probe.ID] = probe
        else:
            print "Warning: probe data is missing attribute named either " + varIDName + " or 'Ratio' or both"
            

    def sendData(self):
        """Compute norm. factors, plot them, normalize data and send out normalized data.
        """
##        self.send("Examples", self.getNormalizedData())


    def sendProbes(self):
        """Sends out example table with currently used probes.
        """
        if self.varsEnum:
            contRatDomain = orange.Domain([self.varsEnum[self.varIdxID], orange.FloatVariable("Ratio"), orange.StringVariable("ColorHSV"), orange.FloatVariable("Symbol")], None)
            contRatET = orange.ExampleTable(contRatDomain)
            pNames = self.probesUsed.keys()
            pNames.sort()
            for pName in pNames:
                probe = self.probesUsed[pName]
    ##            contRatET.append(orange.Example(contRatDomain, [pName, probe.ratio, str(probe.color.name())[1:], probe.symbol]))
                hsvStr = reduce(lambda a,b: a+b, map(lambda x: string.replace("%2s" % hex(x)[2:], " ", "0"), probe.color.hsv()), "")
                contRatET.append(orange.Example(contRatDomain, [pName, probe.ratio, hsvStr, probe.symbol]))
            self.send("Probes", contRatET)
        else:
            self.send("Probes", None)


    ###################################################################################
    ## EVENT HANDLERS: VARIABLE ASSIGNMENT

    def varIDChange(self):
        """Refresh listbox containing other variables and refill self.tblControls.
        """
        self.fillLbVarOthers()
        self.initProbes()
        self.sendData()
        self.sendProbes()

    def varChange(self, memberVarName):
        """Refresh listbox containing other variables.
        """
        self.fillLbVarOthers()
        self.sendData()

    def varOthersChange(self):
        """Updates list of selected other vars (lbVarOthers -> self.varOthersSelectedNames).
        """
        print "        varOthersChange               "
        self.varOthersSelectedNames = []
        for i in range(0, self.lbVarOthers.count()):
            if self.lbVarOthers.isSelected(i):
                self.varOthersSelectedNames.append(self.lbVarOthers.item(i).text())
        self.sendData()

    ###################################################################################
    ## EVENT HANDLERS: CONTROL ASSIGNMENT

    def tblControlsChange(self, row, col):
        """Handles direct changes to self.tblControls;
        updates self.internContRatios and sends out data and control/ratios.
        """
        resendData = True
        cName = str(self.tblControls.item(row, 0).text())
        try:
            ratio = float(eval(str(self.tblControls.item(row, 1).text())))
        except:
            ratio = None
        if ratio:   # not 0 or None
            if self.probesUsed.has_key(cName):
                probe = self.probesUsed[cName]
                probe.ratio = ratio
                probe.color = self.probeColor
                probe.symbol = self.getProbeSymbol()
            else:
                probe = Probe(cName, ratio, self.probeColor, self.getProbeSymbol())
                probe.setData(self.data, self.varsEnum[self.varIdxID], self.varsFloat[self.varIdxSignalSmpl],
                              self.varsFloat[self.varIdxSignalRef], self.varsFloat[self.varIdxBGSmpl],
                              self.varsFloat[self.varIdxBGRef])
                self.probesUsed[cName] = probe
            probe.replotCurve(self.graph, self.subtrBG, self.logAxisX, self.logAxisY, self.markerSize, True)
            self.updateProbeTable(row, probe)
        else:
            if self.probesUsed.has_key(cName):
                probe = self.probesUsed.pop(cName)
                probe.removeCurve(self.graph, True)
            else:
                print "tblControlsChange: not resending data !!!!!!!"
                resendData=False
            self.updateProbeTable(row, None)
        if resendData:
            self.sendData()
            self.sendProbes()


    def tblControlsCurrentChanged(self, row, col):
        """Handles changes of currently selected cell in self.tblControls;
        makes markers of the selected probe thicker, sets self.probeActive.
        """
        print "tblControlsCurrentChanged"
        pName = str(self.tblControls.item(row, 0).text())
        if self.probesUsed.has_key(pName):
            probe = self.probesUsed[pName]
            self.controlRatio = str(probe.ratio)
            self.setProbeSymbol(probe.symbol)
            self.probeColor = probe.color
            # update color of button
            self.btnProbeColor.pixmap().fill(self.probeColor)
            self.btnProbeColor.repaint()
            # activate currently selected probe markers
            self.activateProbe(probe)
        else:
            self.activateProbe(None)
        


    def selectControlsClick(self):
        """Select probes where ID contains self.controlName.
        """
        for idx in range(self.tblControls.numRows()):
            if string.lower(self.controlName) in string.lower(self.tblControls.item(idx, 0).text()):
                sel = QTableSelection()
                sel.init(idx,0)
                sel.expandTo(idx,0)
                self.tblControls.addSelection(sel)

    def selectControlsAllClick(self):
        """Clears all selections and selects all probes.
        """
        self.tblControls.clearSelection(False)
        sel = QTableSelection()
        sel.init(0,0)
        sel.expandTo(self.tblControls.numRows()-1,0)
        self.tblControls.addSelection(sel)


    def probeColorClick(self):
        probeColor = QColorDialog.getColor(self.probeColor, self)
        if probeColor.isValid():
            self.probeColor = probeColor
        self.btnProbeColor.pixmap().fill(self.probeColor)
        self.btnProbeColor.repaint()

        
    def setRatiosClick(self):
        """Sets ratios for the selected controls
        """
        self.setSelectedProbes(self.controlRatio, self.probeColor, self.getProbeSymbol())
        self.sendData()
##        self.sendControlRatios()
        self.sendProbes()


    def clearRatiosClick(self):
        """Clears ratios for the selected controls
        """
        self.setSelectedProbes("", None, None)
        self.sendData()
##        self.sendControlRatios()
        self.sendProbes()


    ###################################################################################
    ## EVENT HANDLERS: SETTINGS

    def settingsFiltersChange(self):
        """Handles changes of filter settings, which affects graph curves and ouput data.
        """
        self.replotProbeCurves()
        self.sendData()


    def settingsNormTypeChange(self):
        """Handles changes of normalization type, which affects normalization curve and output data.
        """
        self.replotNormCurve()
        self.sendData()


    def settingsReplicasChange(self):
        """Handles changes of replicas settings, which affects output data.
        """
        self.rbgMergeReplType.setEnabled(self.mergeRepl)
        self.rbgMergeOtherType.setEnabled(self.mergeRepl)
        self.cbOutNumNonFiltered.setEnabled(self.mergeRepl)
        self.sendData()


    def settingsGraphAxisChange(self, axis):
        """Handles changes of graph axis settings; replot axes and curves.
        """
        self.setGraphAxes([axis])
        self.replotProbeCurves()


    def settingsGraphChange(self):
        """Handles changes of graph settings; replot curves.
        """
        self.replotProbeCurves()


    def settingsOutputChange(self):    
        """Handles changes of output settings; send out data.
        """
        self.sendData()

    ###################################################################
    ## UTILITY FUNCTIONS: GENERAL

    def getProbeSymbol(self):
        return self._probeSymbol + 1


    def setProbeSymbol(self, idx):
        self._probeSymbol = idx - 1


    ###################################################################
    ## UTILITY FUNCTIONS: VARIABLE ASSIGNMENT

    def searchDefaultVarAssignmentClick(self):
        """Select default variables in combos based on their names.
        """
        if self.data:
            # smart variable assignment: ID
            self.varIdxID = None
            for idx,var in enumerate(self.varsEnum):
                if "ID" in var.name:
                    self.varIdxID = idx
                    break
            if self.varIdxID == None:
                self.varIdxID = 0
            # smart variable assignment: signal, background (smpl & ref)
            colNames = ["Raw intensity (med) {Smpl}", "Raw intensity (med) {Ref}", "Background (med) {Smpl}", "Background (med) {Ref}"]
            vars = [None, None, None, None]
            varsFloat = list(self.varsFloat)
            for cIdx, cName in enumerate(colNames):
                for var in varsFloat:
                    if cName in var.name:
                        vars[cIdx] = var
                        varsFloat.remove(var)
                        print cName, "break1"
                        break
                if vars[cIdx] == None:
                    vars[cIdx] = varsFloat[0]
                    varsFloat.pop(0)
                    print cName, "break2"
            # select vars in combos, update listbox with other variables
            self.varIdxSignalSmpl = self.varsFloat.index(vars[0])
            self.varIdxSignalRef = self.varsFloat.index(vars[1])
            self.varIdxBGSmpl = self.varsFloat.index(vars[2])
            self.varIdxBGRef = self.varsFloat.index(vars[3])
            self.fillLbVarOthers()


    def fillCmbVars(self):
        """ Fills combos with variables
        """
        cmbNames = ["cmbVarID", "cmbVarSignalSmpl", "cmbVarSignalRef", "cmbVarBGSmpl", "cmbVarBGRef"]
        varsFloatNames = map(lambda var: var.name, self.varsFloat)
        varsNames = [map(lambda var: var.name, self.varsEnum), varsFloatNames, varsFloatNames, varsFloatNames, varsFloatNames]
        if self.data:
            for cmbName, vNames in zip(cmbNames, varsNames):
                self.__dict__[cmbName].clear()
                self.__dict__[cmbName].insertStrList(vNames)
        else:
            for cmbName in cmbNames:
                self.__dict__[cmbName].clear()


    def fillLbVarOthers(self):
        """Fills listBox with variables not selected by combos
        """
        self.lbVarOthers.clear()
        if self.data:
            varNames = map(lambda var: var.name, self.data.domain.variables)
            if self.varsEnum[self.varIdxID].name in varNames:
                varNames.remove(self.varsEnum[self.varIdxID].name)
            if self.varsFloat[self.varIdxSignalSmpl].name in varNames:
                varNames.remove(self.varsFloat[self.varIdxSignalSmpl].name)
            if self.varsFloat[self.varIdxSignalRef].name in varNames:
                varNames.remove(self.varsFloat[self.varIdxSignalRef].name)
            if self.varsFloat[self.varIdxBGSmpl].name in varNames:
                varNames.remove(self.varsFloat[self.varIdxBGSmpl].name)
            if self.varsFloat[self.varIdxBGRef].name in varNames:
                varNames.remove(self.varsFloat[self.varIdxBGRef].name)
            self.lbVarOthers.insertStrList(varNames)
            # select items (self.lbVarOthers <- self.varOthersSelectedNames)
            for idx in range(self.lbVarOthers.count()):
##                print self.lbVarOthers.item(idx).text(), self.varOthersSelectedNames, self.lbVarOthers.item(idx).text() in self.varOthersSelectedNames
                self.disconnect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
                self.lbVarOthers.setSelected(idx, self.lbVarOthers.item(idx).text() in self.varOthersSelectedNames)
                self.connect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)


    ###################################################################
    ## UTILITY FUNCTIONS: CONTROLS ASSIGNMENT

    def initProbes(self):
        """Clears self.tblControls and self.probesUsed;
        copies ratios from self.probesExternal to self.probesUsed (only those that are in self.tblControls);
        initializes self.tblControls with values of the selected attribute representing probe IDs;
        fills probe ratios and plots curves
        """
        print "initProbes"
        # clear self.probesUsed and curves
        for probe in self.probesUsed.values():
            probe.removeCurve(self.graph, False)
        self.graph.replot()
        self.probesUsed = {}
##        self.probeActive = None
        if self.data and self.varIdxID != None:
            idList = list(self.varsEnum[self.varIdxID].values)
            idList.sort()
            # self.tblControls, set horizontal header labels
            self.tblControls.horizontalHeader().setLabel(0, self.varsEnum[self.varIdxID].name)
            self.tblControls.setNumRows(len(idList))
            firstFilledRow = -1
            for row, id in enumerate(idList):
                OWGUI.tableItem(self.tblControls, row, 0, id, editType=QTableItem.Never)#, background=QColor(160,160,160))
##                self.tblControls.setItem(row, 0, QTableItem(self.tblControls, QTableItem.Never, id, self.symbolPixmapEmpty))
                probeExt = self.probesExternal.get(id, None)
                if probeExt:
                    probeExt.setData(self.data, self.varsEnum[self.varIdxID], self.varsFloat[self.varIdxSignalSmpl],
                                     self.varsFloat[self.varIdxSignalRef], self.varsFloat[self.varIdxBGSmpl], self.varsFloat[self.varIdxBGRef])
                    self.probesUsed[probeExt.ID] = probeExt
                    # update self.tblControls
                    self.updateProbeTable(row, probeExt)
                    lastFilledRow = row
                    if firstFilledRow == -1:
                        firstFilledRow = row
                    # plot curve
                    probeExt.replotCurve(self.graph, self.subtrBG, self.logAxisX, self.logAxisY, self.markerSize, False)
                else:
                    self.updateProbeTable(row, None)
            # adjust table width and select a cell, replot graph
            self.tblControls.adjustColumn(0)
            self.tblControls.setColumnWidth(1, max(35, self.tblControls.visibleWidth() - self.tblControls.columnWidth(0) - self.tblControls.columnWidth(2)))
            if firstFilledRow != -1:
##                print "setcurrcell start"
                self.tblControls.setCurrentCell(firstFilledRow, 1)
                self.tblControlsCurrentChanged(firstFilledRow, 1)
##                print "setcurrcell end"
##                self.activeProbe = self.probesUsed[str(self.tblControls.item(firstFilledRow, 0).text())]
##                self.activeProbe.setCurveActive(self.graph, True)
                self.graph.replot()
            elif len(idList) > 0:
                self.tblControls.setCurrentCell(0, 1)
        else:
            self.tblControls.horizontalHeader().setLabel(0, "Probe ID")
            self.tblControls.setNumRows(0)
            

    def setSelectedProbes(self, ratioStr, color, symbol):
        """Updates self.probesUsed of the selected probes with a given ratio string, color and symbol;
        calls updateProbeTable(row, probe), replots/removes curves
        """
        for selNum in range(self.tblControls.numSelections()):
            sel = self.tblControls.selection(selNum)
##            print "  sel:", sel.topRow(), sel.bottomRow()+1
            for row in range(sel.topRow(), sel.bottomRow()+1):
##                OWGUI.tableItem(self.tblControls, row, 1, ratioStr, editType=QTableItem.Always, background=QColor(160,160,160))
##                self.tblControls.setText(row, 1, ratioStr)
                cName = str(self.tblControls.item(row, 0).text())
                try:
                    ratio = float(eval(ratioStr))
                except:
                    ratio = None
                if ratio:
                    if self.probesUsed.has_key(cName):
                        probe = self.probesUsed[cName]
                        probe.ratio = ratio
                        probe.color = color
                        probe.symbol = symbol
                    else:
                        probe = Probe(cName, ratio, color, symbol)
                        self.probesUsed[cName] = probe
                        probe.setData(self.data, self.varsEnum[self.varIdxID], self.varsFloat[self.varIdxSignalSmpl],
                                                   self.varsFloat[self.varIdxSignalRef], self.varsFloat[self.varIdxBGSmpl],
                                                   self.varsFloat[self.varIdxBGRef])
                    probe.replotCurve(self.graph, self.subtrBG, self.logAxisX, self.logAxisY, self.markerSize, False)
                    self.updateProbeTable(row, self.probesUsed[cName])
                else:
                    print "setSelectedProbes: remove probe & curve"
                    if self.probesUsed.has_key(cName):
                        probe = self.probesUsed.pop(cName)
                        probe.removeCurve(self.graph, False)
                        if self.probeActive == probe:
                            self.activateProbe(None)
                    self.updateProbeTable(row, None)
        # activate probes (if there was no probe in the current row before)
        cName = str(self.tblControls.item(self.tblControls.currentRow(), 0).text())
        if self.probesUsed.has_key(cName):
            self.activateProbe(self.probesUsed[cName])
        # replot if there were any changes
        if self.tblControls.numSelections() > 0:
            self.graph.replot()


    def updateProbeTable(self, row, probe):
        """Updates self.tblControls with data from the probe (columns 1 and 2).
        """
        if probe:
            txt = str(probe.ratio)
            pxm = probe.getSymbolPixmap(self.tblControls.cellGeometry(row,2))
####            self.tblControls.setText(row, 1, str(probe.ratio))
##            item1 = QTableItem(self.tblControls, QTableItem.Always, str(probe.ratio))
##            self.tblControls.setItem(row, 1, item1)
##            item2 = QTableItem(self.tblControls, QTableItem.Never, "", probe.getSymbolPixmap(self.tblControls.cellGeometry(row,2)))
##            self.tblControls.setItem(row, 2, item2)
        else:
            txt = ""
            pxm = self.symbolPixmapEmpty
        self.tblControls.setItem(row, 1, QTableItem(self.tblControls, QTableItem.OnTyping, txt))
##            self.tblControls.setText(row, 1, "")
        self.tblControls.setItem(row, 2, QTableItem(self.tblControls, QTableItem.Never, "", pxm))


    def activateProbe(self, probe):
        """Deactive currently activated probe;
        activate new probe or None.
        """
##        if self.probeActive != probe:
        if self.probeActive:
            self.probeActive.setCurveActive(self.graph, False)
        if probe:
            probe.setCurveActive(self.graph, True)
        self.probeActive = probe

    ###################################################################
    ## UTILITY FUNCTIONS: NORMALIZATION


    ###################################################################
    ## UTILITY FUNCTIONS: PLOT GRAPH

    def setGraphAxes(self, axes=None):
        """According to selected scaling sets up axis labels and scales;
        axis: 0: vertical left, 1: vertical right, 2: horizontal bottom, 3: horizontal top
        """
        titles = {False: ["Normalization factor"]*2 + ["Average intensity"]*2, True: ["log2(normalization factor)"]*2 + ["log2(average intensity)"]*2}
        useLog = [self.logAxisY]*2 + [self.logAxisX]*2
        if axes==None: axes = [0,2]
        for axis in axes:
            self.graph.setAxisTitle(axis, titles[useLog[axis]][axis])
    

    def replotProbeCurves(self):
        """Replots all probe curves.
        """
        print "replotProbeCurves"
        for pID, probe in self.probesUsed.items():
            probe.replotCurve(self.graph, self.subtrBG, self.logAxisX, self.logAxisY, self.markerSize, False)
##        self.activateProbe(self.probeActive)
        if self.probeActive:
            self.probeActive.setCurveActive(self.graph, True)
        self.graph.replot()


    def replotNormCurve(self):
        """Replots normalization curve.
        """
        print "TODO"
        

##    lr = MultLinReg(Numeric.reshape(Numeric.asarray(As), (len(As),1)), SFs)
##    if verbose:
##        # loess
##        loessAs = Numeric.arange(min(As), max(As), (max(As)-min(As))/100.).tolist()
##        loessSFs = Numeric.asarray(statc.loess(zip(As, SFs), loessAs, windowSize))[:,1]
##        plotedItems.append(P.plot(loessAs, loessSFs, c="black", linewidth=2))
##        # lin. reg.
##        plotedItems.append(P.plot(As, lr.y_hat, c="blue", linewidth=2))
##        # median
##        plotedItems.append(P.plot(As, Numeric.resize(NumExtn.median(SFs), (len(As),)), c="red", linewidth=2))
##        # legend
##        P.figlegend(plotedItems, contUsed + ["loess", "lin.reg", "median"], 'upper right', numpoints=2, markerscale=1)
        
    


                
# test widget appearance
if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWNormalize()
    a.setMainWidget(ow)
    ow.show()
    ow.onDataInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\STEROLTALK\array-pro\experiments\Tadeja 2nd image analysis\10vs10mg raw data\chol 0560.tab"))
    ow.onProbesInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\Orange\OWNormalize\steroltalk v0 controlRatios.tab"))
    a.exec_loop()
    #save settings 
    ow.saveSettings()
