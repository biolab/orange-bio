"""
<name>Normalize Expression Array Data</name>
<description>Normalization of expression array data.</description>
<icon>icons/Normalize.png</icon>
<priority>1150</priority>
<author>Peter Juvan (peter.juvan@fri.uni-lj.si)</author>
"""

"""
TODO: settingsList
"""

import string, math
import Numeric, MA, MLab
import LinearAlgebra
import NumExtn
import orange
from qttable import *
from OWWidget import *
import OWGUI, OWToolbars
from OWGraph import *
from OWGraphTools import *      # color palletes, user defined curves, ...
import ColorPalette             # ColorButton
import types

import chipstat



# global debugging variables
D1 = False
D2 = False
D3 = False
D4 = False
D5 = False
D6 = False # Probes

class OWNormalize(OWWidget):

    settingsList = ["varNameA", "varNameB", "varNameSignalSmpl", "varNameSignalRef", "varNameBGSmpl", "varNameBGRef", "varNameBGSmplSD", "varNameBGRefSD",
                    "defNameA", "defNameB", "defNameSmpl1", "defNameSmpl2", "defNameRef1", "defNameRef2", "defNameForeground", "defNameBackground", "defNameMean", "defNameSD",
                    "displayVarAAliases", "commitOnChange"]
    # constants
    stylesIndices = zip(["<none>","Circle","Rect","Diamond","Triangle","DTriangle","UTriangle","LTriangle","RTriangle","Cross","XCross"], range(11))
    sizeButtonColor = 18
    # self.tblControls column indices
    tcPKey = 0      # hidden column for probe keys
    tcMarker = 1
    tcRatio = 2
    tcNPAll = 3         # number of all probes
    tcNPAccepted = 4    # number of non-filtered probes
    tcVarA = 5
    tcVarAAlias = 6      # hidden if not self.displayVarAAliases
    tcVarB = 7          # hidden if probe varB is <none>
    # merge level
    MergeLevelNone = 0
    MergeLevelPerVarsAB = 1
    MergeLevelPerVarA = 2
    # mergeOtherTypes
    MergeOtherTypeMean = 0
    MergeOtherTypeMedian = 1
    MergeOtherTypeConc = 2
    # normalization curve styles
    normCurveStyles = {0: QwtCurve.Lines, 1: QwtCurve.Spline}


    def __init__(self, parent = None, signalManager = None, name = "Normalize Microarray Data"):
        if D1 or D6: print "OWNormalize.__init__"
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, wantStatusBar=False)  #initialize base class
        # store original caption title for extending it with expression and probe data file names
        self.captionTitleBase = self.captionTitle

        # set channels
        self.inputs = [("Expression Data", ExampleTable, self.onDataInput), ("Probe Data", ExampleTable, self.onProbesInput)]
        self.outputs = [("Expression Data", ExampleTable), ("Probe Data", ExampleTable)]

        # defaults: filters
        self._def_subtrBG = 0
        self._def_useCV = True
        self._def_maxCV = 0.5
        self._def_useMinIntensity = True
        self._def_minIntensityRatio = 1.5
        self._def_useMaxFGIntensity = True
        self._def_maxFGIntensity = 60000
        self._def_useMaxBGIntensity = False
        self._def_maxBGIntensity = 400
        # defaults: normalization
        self._def_arrayType = 0  # 0: whole-genome array, 1: custom array with control probes
        self._def_normRange = 2  # 0: global, 1: local, per var B values, 2: combined
        self._def_minNumControlProbes = 2
        self._def_approxFunction = 2   # 0: median, 1: LR, 2: LOESS
        self._def_loessWindow = 60
        self._def_loessWeight = 0.0

        # general settings
        self.controlName = ""
        self.ratioStr = ""
        self.probeSymbolIdx = 0    # index of selected cmbProbeSymbol item
        self.probeColor = ProbeSet.NoColor
        # general settings: filters
        self.subtrBG = self._def_subtrBG
        self.useCV = self._def_useCV
        self.maxCV = self._def_maxCV
        self.useMinIntensity = self._def_useMinIntensity
        self.minIntensityRatio = self._def_minIntensityRatio
        self.useMaxFGIntensity = self._def_useMaxFGIntensity
        self.maxFGIntensity = self._def_maxFGIntensity
        self.useMaxBGIntensity = self._def_useMaxBGIntensity
        self.maxBGIntensity = self._def_maxBGIntensity
        # general settings: normalization
        self.arrayType = self._def_arrayType
        self.normRange = self._def_normRange
        self.minNumControlProbes = self._def_minNumControlProbes
        self.approxFunction = self._def_approxFunction
        self.loessWindow = self._def_loessWindow
        self.loessWeight = self._def_loessWeight

        # settings
        self.logAxisX = True
        self.logAxisY = True
        self.markerSize = 9
        self.mergeReplGraph = False
        self.showLegend = False
        self.tracking = True
        self.normCurveStyleIdx = 0
        self.displayVarAAliases = True
        self.commitOnChange = True

        # output
        self.mergeLevel = OWNormalize.MergeLevelPerVarA             #0: none, 1: per vars A & B (ID & Type), 2: per var A (ID)
        self.mergeIntensitiesType = 1   #0: mean, 1: median
        self.mergeOtherType = OWNormalize.MergeOtherTypeMedian        # Mean, Median, Concatenate
        self.mergeOtherRemoveDupl = False                           # whether or not to remove duplicates from comma separated list of values
        
        self.outVarAAliases = True
        self.outNumProbes = True
        self.outNetSignal = True
        self.outA = True
        self.outMRaw = True
        self.outMCentered = True
##        self.outNormLogRatioCentered = False
        self.autoSendSelection = 1

        # context-specific settings
        self.data = None
        self.dataProbes = None
        self.varsAll = {}
        self.varsFloat = {}
        self.varsEnumStr = {}
        # variable indices selected within combos
        self.varNameA = None
        self.varNameB = "<none>"
        self.varNameSignalSmpl = None
        self.varNameSignalRef = None
        self.varNameBGSmpl = None
        self.varNameBGRef = None
        self.varNameBGSmplSD = "<none>"
        self.varNameBGRefSD = "<none>"
        # default var names
        self.defNameA = "id"
        self.defNameB = ""
        self.defNameSmpl1 = "smpl"
        self.defNameSmpl2 = ""
        self.defNameRef1 = "ref"
        self.defNameRef2 = ""
        self.defNameForeground = "raw"
        self.defNameBackground = "background"
        self.defNameMean = "med"
        self.defNameSD = "st.dev"
        # names of selected vars from listBox
        self.varsOtherSelected = {}

        # load previous settings
        self.loadSettings()

        # GUI
        self.controlArea.setFixedWidth(265)
        self.controlArea.setFixedWidth(275)
        # main area: graph
        boxG = QVBox(self.mainArea)
        graphBoxLayout = QVBoxLayout(self.mainArea)
        graphBoxLayout.addWidget(boxG)
        self.graph = OWGraphMA(boxG)
        self.graph.setAutoReplot(False)
        self.setGraphAxes(axes=[0,2])
        self.settingsProbeTrackingChange()  # connect events to self.graph
        self.connect(self.graph, SIGNAL("legendClicked(long)"), self.onLegendClicked)
        self.graph.enableGraphLegend(self.showLegend)
        # save graph button
        self.connect(self.graphButton, SIGNAL("clicked()"), self.graph.saveToFile)
        
        # control area: tabs
        self.tabs = QTabWidget(self.controlArea, 'tabWidget')
        # tab 1: vars
        boxVars = QVGroupBox(self)
        self.tabs.insertTab(boxVars, "Var")
        boxGroupBy = QVGroupBox('Group probes by', boxVars)
        self.cmbVarA = OWGUI.comboBox(boxGroupBy, self, "varNameA", label="ID", labelWidth=33, orientation="horizontal", callback=self.varABChange, sendSelectedValue=1, valueType=str)
        ### Var B, Type, Other, Additional, Optional, Opt.var, Extra, Alternative var, Alt.var 
        self.cmbVarB = OWGUI.comboBox(boxGroupBy, self, "varNameB", label="Alt.var", labelWidth=33, orientation="horizontal", callback=self.varABChange, sendSelectedValue=1, valueType=str)
        boxFGI = QVGroupBox('Average foreground intensity', boxVars)
        self.cmbVarSignalSmpl = OWGUI.comboBox(boxFGI, self, "varNameSignalSmpl", label="Smpl", labelWidth=33, orientation="horizontal", callback=self.varFBChange, sendSelectedValue=1, valueType=str)
        self.cmbVarSignalRef = OWGUI.comboBox(boxFGI, self, "varNameSignalRef", label="Ref", labelWidth=33, orientation="horizontal", callback=self.varFBChange, sendSelectedValue=1, valueType=str)
        boxBGI = QVGroupBox('Average background intensity', boxVars)
        self.cmbVarBGSmpl = OWGUI.comboBox(boxBGI, self, "varNameBGSmpl", label="Smpl", labelWidth=33, orientation="horizontal", callback=self.varFBChange, sendSelectedValue=1, valueType=str)
        self.cmbVarBGRef = OWGUI.comboBox(boxBGI, self, "varNameBGRef", label="Ref", labelWidth=33, orientation="horizontal", callback=self.varFBChange, sendSelectedValue=1, valueType=str)
        boxBGSD = QVGroupBox('Background std. deviation', boxVars)
        self.cmbVarBGSmplSD = OWGUI.comboBox(boxBGSD, self, "varNameBGSmplSD", label="Smpl", labelWidth=33, orientation="horizontal", callback=self.varFBChange, sendSelectedValue=1, valueType=str)
        self.cmbVarBGRefSD = OWGUI.comboBox(boxBGSD, self, "varNameBGRefSD", label="Ref", labelWidth=33, orientation="horizontal", callback=self.varFBChange, sendSelectedValue=1, valueType=str)
        # tab 1: default var names
        boxDefaultNames = QVGroupBox('Default variable names', boxVars)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameA", label="ID", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameB", label="Alt.var", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameSmpl1", label="Smpl 1", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameSmpl2", label="Smpl 2", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameRef1", label="Ref 1", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameRef2", label="Ref 2", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameForeground", label="Foreground", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameBackground", label="Background", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameMean", label="Average", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.lineEdit(boxDefaultNames, self, "defNameSD", label="Std. deviation", labelWidth=70, orientation='horizontal', box=None, tooltip=None)
        OWGUI.button(boxDefaultNames, self, "Search for Default Variables", callback=self.defaultVarAssignmentClick)

        # tab 4: normalization
        boxNorm = QVGroupBox(self)
        self.tabs.insertTab(boxNorm, "Norm")
        # tab 4: normalization: array type
        self.boxArrayType = OWGUI.radioButtonsInBox(boxNorm, self, value="arrayType", box='Expression array type', btnLabels=["Genome array", "Boutique array with control probes"], callback=self.settingsArrayTypeChange)
        # tab 4: normalization: range, type
        self.boxNormRange = OWGUI.radioButtonsInBox(boxNorm, self, value="normRange", box='Normalization range', btnLabels=["Global, entire microarray", "Per type of probe", "Combined (w.r.t. num. of control probes)"], callback=self.settingsNormalizationChange)
        sldMinNumControlProbes = OWGUI.qwtHSlider(self.boxNormRange, self, "minNumControlProbes", box="Min. number of control probes", minValue=2, maxValue=300, step=1, precision=0, callback=None, logarithmic=0, ticks=0, maxWidth=110)
        self.connect(sldMinNumControlProbes, SIGNAL("sliderReleased()"), self.settingsNormalizationChange)
        boxapproxFunction = OWGUI.radioButtonsInBox(boxNorm, self, value="approxFunction", box='Approx. function', btnLabels=["Median (intensity independent)", "Linear regression", "Loess"], callback=self.settingsNormalizationChange)
        # tab 4: normalization type: loess settings
        sldLoessWindow = OWGUI.qwtHSlider(boxapproxFunction, self, "loessWindow", box="Window size (% of points)", minValue=1, maxValue=99, step=1, precision=0, logarithmic=0, ticks=0, maxWidth=110)
        self.connect(sldLoessWindow, SIGNAL("sliderReleased()"), self.settingsNormalizationChange)
        boxSldLoessWeight = QHBox(boxapproxFunction)
        sldLoessWeight = OWGUI.qwtHSlider(boxSldLoessWeight, self, "loessWeight", box="Weight of non-control probes [0,1]", minValue=0, maxValue=1, step=0.01, precision=2, logarithmic=0, ticks=0, maxWidth=110)
        self.connect(sldLoessWeight, SIGNAL("sliderReleased()"), self.settingsNormalizationChange)
        boxSldLoessWeight.setEnabled(False)
        # tab 4: default button
        OWGUI.button(boxNorm, self, "Set &Default Values", callback=self.normalizationAllChange)

        # tab 3: filters
        boxFilters = QVGroupBox(self)
        self.tabs.insertTab(boxFilters, "Filter")
        # tab 3: filters: subtract BG
        self.cbSubtrBG = OWGUI.checkBox(boxFilters, self, "subtrBG", "Subtract background", callback=self.settingsSubstrBGChange)
        # tab 3: filters: CV
        self.boxMaxCV = QVGroupBox('Max. coeff. of variation (CV)', boxFilters)
        OWGUI.checkBox(self.boxMaxCV, self, "useCV", "Enabled", callback=self.settingsFilterMaxCVChange)
        sldMaxCV = OWGUI.qwtHSlider(self.boxMaxCV, self, "maxCV", minValue=0, maxValue=2, step=0.01, precision=2, callback=None, logarithmic=0, ticks=0, maxWidth=110)
        self.connect(sldMaxCV, SIGNAL("sliderReleased()"), self.settingsFilterMaxCVChange)
        self.lblInfoFilterMaxCV = QLabel("\n", self.boxMaxCV)
        # tab 3: filters: minIntensityRatio
        boxMinIntRatio = QVGroupBox('Min. signal to background ratio', boxFilters)
        OWGUI.checkBox(boxMinIntRatio, self, "useMinIntensity", "Enabled", callback=self.settingsFilterMinIntRatioChange)
        sldMinInt = OWGUI.qwtHSlider(boxMinIntRatio, self, "minIntensityRatio", minValue=0, maxValue=5, step=0.01, precision=2, callback=None, logarithmic=0, ticks=0, maxWidth=110)
        self.connect(sldMinInt, SIGNAL("sliderReleased()"), self.settingsFilterMinIntRatioChange)
        self.lblInfoFilterMinIntRatio = QLabel("\n", boxMinIntRatio)
        # tab 3: filters: maxFGIntensity
        boxMaxFGIntensity = QVGroupBox('Max. foreground intensity', boxFilters)
        OWGUI.checkBox(boxMaxFGIntensity, self, "useMaxFGIntensity", "Enabled", callback=self.settingsFilterMaxFGIntChange)
        sldMaxFGInt = OWGUI.qwtHSlider(boxMaxFGIntensity, self, "maxFGIntensity", minValue=0, maxValue=65536, step=1, precision=0, callback=None, logarithmic=0, ticks=0, maxWidth=110)
        self.connect(sldMaxFGInt, SIGNAL("sliderReleased()"), self.settingsFilterMaxFGIntChange)
        self.lblInfoFilterMaxFGInt = QLabel("\n", boxMaxFGIntensity)
        # tab 3: filters: maxBGIntensity
        boxMaxBGIntensity = QVGroupBox('Max. background intensity', boxFilters)
        OWGUI.checkBox(boxMaxBGIntensity, self, "useMaxBGIntensity", "Enabled", callback=self.settingsFilterMaxBGIntChange)
        sldMaxBGInt = OWGUI.qwtHSlider(boxMaxBGIntensity, self, "maxBGIntensity", minValue=0, maxValue=4096, step=1, precision=0, callback=None, logarithmic=0, ticks=0, maxWidth=110)
        self.connect(sldMaxBGInt, SIGNAL("sliderReleased()"), self.settingsFilterMaxBGIntChange)
        self.lblInfoFilterMaxBGInt = QLabel("\n", boxMaxBGIntensity)
        # tab 3: default button
        OWGUI.button(boxFilters, self, "Set &Default Values", callback=self.filtersAllChange)

        # tab 2: table probe/ratio/marker
        boxProbes = QVGroupBox(boxVars)
        self.tabs.insertTab(boxProbes, "Probe")
        self.tblControls = QTable(boxProbes)
        self.tblControls.setNumCols(8)
        self.tblControls.setColumnWidth(OWNormalize.tcMarker, 20)
        self.tblControls.setColumnWidth(OWNormalize.tcRatio, 15)
        self.tblControls.setColumnWidth(OWNormalize.tcNPAll, 10)
        self.tblControls.setColumnWidth(OWNormalize.tcNPAccepted, 10)
        self.tblControls.setColumnWidth(OWNormalize.tcVarA, 15)
        self.tblControls.setColumnWidth(OWNormalize.tcVarAAlias, 0)
        self.tblControls.setColumnWidth(OWNormalize.tcVarB, 0)
        self.connect(self.tblControls, SIGNAL("valueChanged(int,int)"), self.tblControlsValueChange)
        self.connect(self.tblControls, SIGNAL("currentChanged(int, int)"), self.tblControlsCurrentChanged)
        self.connect(self.tblControls , SIGNAL('selectionChanged()'), self.tblControlsSelectionChanged)
        self.connect(self.tblControls, SIGNAL("doubleClicked(int, int, int, const QPoint &)"), self.tblControlsDoubleClicked)
        hheader=self.tblControls.horizontalHeader()
        self.connect(hheader,SIGNAL("clicked(int)"), self.tblControlsHHeaderClicked)
        self.sortby = 0
        hheader.setLabel(OWNormalize.tcMarker, "")
        hheader.setLabel(OWNormalize.tcRatio, "Ratio")
        hheader.setLabel(OWNormalize.tcNPAll, "##")
        hheader.setLabel(OWNormalize.tcNPAccepted, "#")
        hheader.setLabel(OWNormalize.tcVarA, "ID")
        hheader.setLabel(OWNormalize.tcVarAAlias, "Alias")
        hheader.setLabel(OWNormalize.tcVarB, "Alt.var")
        hheader.setMovingEnabled(False)
        # hide vertical header and columns pKey, name
        self.tblControls.setLeftMargin(0)
        self.tblControls.verticalHeader().hide()
        self.tblControls.hideColumn(OWNormalize.tcPKey)
        self.tblControls.hideColumn(OWNormalize.tcVarAAlias)
        self.tblControls.hideColumn(OWNormalize.tcVarB)
        # tab 2: buttons
        boxBtns0 = QVGroupBox("Select probes where ID contains", boxProbes)
        boxBtns00 = QHBox(boxBtns0)
        OWGUI.lineEdit(boxBtns00, self, "controlName")
        OWGUI.button(boxBtns00, self, "Select", callback=self.btnSelectControlsClick)
        boxBtns01 = QHBox(boxBtns0)
        OWGUI.button(boxBtns01, self, "Select all", callback=self.btnSelectControlsAllClick)
        OWGUI.button(boxBtns01, self, "Unselect all", callback=self.btnUnselectControlsAllClick)
        boxBtns1 = QVGroupBox("Set ratio and marker for selected probes", boxProbes)
        boxBtns11 = QHBox(boxBtns1)
        pxm = QPixmap(OWNormalize.sizeButtonColor,OWNormalize.sizeButtonColor)
        pxm.fill(self.probeColor)
        self.cmbProbeSymbol = OWGUI.comboBox(boxBtns11, self, "probeSymbolIdx")#, callback=self.cmbProbeSymbolActivated)
        for styleName, styleIdx in OWNormalize.stylesIndices:
            symbol = QwtSymbol(styleIdx, QBrush(QColor(255,255,255), QBrush.SolidPattern), QPen(QColor(0,0,0),1), QSize(self.markerSize,self.markerSize))
            pixmap = QPixmap(14,14)
            pixmap.fill(QColor(255,255,255))
            painter = QPainter(pixmap)
            symbol.draw(painter, QPoint(7,7))
            painter.end()
            self.cmbProbeSymbol.insertItem(pixmap, styleName)
        self.btnProbeColor = OWGUI.button(boxBtns11, self, "", callback=self.probeColorClick)
        self.btnProbeColor.setPixmap(pxm)
        leRatio = OWGUI.lineEdit(boxBtns11, self, "ratioStr")
        self.connect(leRatio, SIGNAL("returnPressed()"), self.leRatioReturnPressed)
        boxBtns12 = QHBox(boxBtns1)
        OWGUI.button(boxBtns12, self, "Set", callback=self.btnSetProbesClick)
        OWGUI.button(boxBtns12, self, "Clear", callback=self.btnClearProbesClick)

        # tab 6: output
        boxOutput = QVGroupBox(self)
        self.tabs.insertTab(boxOutput, "Out")
        # tab 6: output: merge replicas
        boxMerge = QVGroupBox('Merge replicas', boxOutput)
        OWGUI.radioButtonsInBox(boxMerge, self, value="mergeLevel", btnLabels=["None", "ID &  Alt.var", "ID"], box="Group probes by matching variable(s)", callback=self.settingsOutputReplicasChange)
        self.rbgMergeIntensitiesType = OWGUI.radioButtonsInBox(boxMerge, self, value="mergeIntensitiesType", btnLabels=["Mean", "Median"], box="Average calculation", callback=self.settingsOutputReplicasChange)
        # tab 6: output: additional info
        boxAdditional = QVGroupBox('Additional info', boxOutput)
        self.cbOutVarAAliases = OWGUI.checkBox(boxAdditional, self, "outVarAAliases", "ID alias", callback=self.settingsOutputChange)
        self.cbOutNumProbes = OWGUI.checkBox(boxAdditional, self, "outNumProbes", "Number of probes", callback=self.settingsOutputChange)
        OWGUI.checkBox(boxAdditional, self, "outNetSignal", "Net intensities", callback=self.settingsOutputChange)
        OWGUI.checkBox(boxAdditional, self, "outA", "A (log2 average intensity)", callback=self.settingsOutputChange)
        OWGUI.checkBox(boxAdditional, self, "outMRaw", "M raw", callback=self.settingsOutputChange)
        OWGUI.checkBox(boxAdditional, self, "outMCentered", "M centered", callback=self.settingsOutputChange)
##        OWGUI.checkBox(boxAdditional, self, "outNormLogRatioCentered", "M normalized centered", callback=self.settingsOutputChange)
        # tab 6: output: other variables
        boxOtherVars = QVGroupBox('Other variables', boxOutput)
        self.lbVarOthers = QListBox(boxOtherVars)
        self.lbVarOthers.setSelectionMode(QListBox.Multi)
        self.connect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
        self.boxMergeOtherType = QVGroupBox("Merge", boxOutput)
        self.boxMergeOtherType.setEnabled(self.mergeLevel and len(self.varsOtherSelected) > 0)
        rbgMergeOtherType = OWGUI.radioButtonsInBox(self.boxMergeOtherType, self, value="mergeOtherType", btnLabels=["Mean", "Median", "Concatenate values"], box="Continuous variables", callback=self.settingsOutputOtherChange)
        boxMergeOtherTypeD = QVGroupBox('Non-continuous variables', self.boxMergeOtherType)
        QLabel("Values are concatenated by default.", boxMergeOtherTypeD)
        self.cbMergeOtherRemoveDupl = OWGUI.checkBox(boxMergeOtherTypeD, self, "mergeOtherRemoveDupl", "Remove duplicate values", callback=self.settingsOutputOtherChange)

        # tab 5: settings
        boxSettings = QVGroupBox(self)
        self.tabs.insertTab(boxSettings, "Settings")
        # tab 5: settings: graph
        boxGraph = QVGroupBox('Graph', boxSettings)
        boxMSize = QHBox(boxGraph)
        QLabel("Marker size", boxMSize)
        cmbMarkerSize = OWGUI.comboBox(boxMSize, self, "markerSize", callback=self.settingsGraphChange, sendSelectedValue=1, valueType=int)
        for itemIdx, size in enumerate(range(3,16)):
            cmbMarkerSize.insertItem(str(size))
            if self.markerSize == size:
                cmbMarkerSize.setCurrentItem(itemIdx)
        OWGUI.checkBox(boxGraph, self, "logAxisX", "Logarithmic X axis", callback=lambda ax=2: self.settingsGraphAxisChange(ax))
        OWGUI.checkBox(boxGraph, self, "logAxisY", "Logarithmic Y axis", callback=lambda ax=0: self.settingsGraphAxisChange(ax))
        cbMergeReplicas = OWGUI.checkBox(boxGraph, self, value="mergeReplGraph", label="Merge replicas", callback=self.settingsGraphChange)
        cbMergeReplicas.setEnabled(False)
        OWGUI.checkBox(boxGraph, self, value="showLegend", label="Show legend", callback=self.settingsShowLegendChange)
        OWGUI.checkBox(boxGraph, self, value="tracking", label="Tracking", callback=self.settingsProbeTrackingChange)
        boxNormCurveStyle = OWGUI.radioButtonsInBox(boxGraph, self, box='Curve style', value="normCurveStyleIdx", btnLabels=["Line", "Spline"], callback=self.settingsGraphChange)
        self.zoomSelectToolbar = OWToolbars.ZoomSelectToolbar(self, boxGraph, self.graph, self.autoSendSelection)
        # tab 5: settings: Probes
        boxProbes = QVGroupBox("Probes", boxSettings)
        OWGUI.checkBox(boxProbes, self, 'displayVarAAliases', 'Display ID aliases', callback=self.adjustProbeTableColumns)
        # tab 5: settings: commit
        boxCommit = QVGroupBox("Output", boxSettings)
        OWGUI.checkBox(boxCommit, self, 'commitOnChange', 'Commit data on selection change', callback=self.commitChange)
        self.btnCommit = OWGUI.button(boxCommit, self, "Commit", callback=self.commitClicked, disabled=1)
        
        # control area: info
        boxProbeInfo = QVGroupBox("Info", self.controlArea)
        self.lblProbeInfo = QLabel("\n\n", boxProbeInfo)

        self.resize(1000, 752)

        # INITIALIZATION: controls/ratios, probe info, filters, filter info
        self.probes = Probes(self.graph, self.subtrBG, self.logAxisX, self.logAxisY, self.markerSize, OWNormalize.normCurveStyles[self.normCurveStyleIdx])
        self.setInfoProbes()
        self.setProbeFilters()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()
        self.probes.setNormalizationParameters(self.normRange, self.minNumControlProbes, self.approxFunction, self.loessWindow, self.loessWeight)


    def updateCaptionTitle(self):
        """updates caption title accorting to the expression and probe data filenames
        """
        ct = self.captionTitleBase
        ctMid = ""
        ctEnd = ""
        exprDataName = ""
        if self.data:
            exprDataName = self.data.name
        probeDataName = ""
        if self.dataProbes:
            probeDataName = self.dataProbes.name
        if exprDataName or probeDataName:
            ct += " ["
            ctEnd = "]"
        if exprDataName and probeDataName:
            ctMid = ", "
        if exprDataName:
            ct += "E: " + exprDataName
        ct += ctMid
        if probeDataName:
            ct += "P: " + probeDataName
        self.setCaptionTitle(ct + ctEnd)
        



    ###################################################################################
    ## EXPRESSION DATA IN / OUT, SEARCH DEFAULT VARS CLICK
    ###################################################################################

    def onDataInput(self, data):
        """Handles input of new data:
            - put metas among variables
            - fill combos, select default variables
            - initialize self.probes
        """
        if D1 or D2: print "OWNormalize.onDataInput"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        # move metas among normal variables, remove string variables and variables with duplicate names, store to self.data
        self.varsAll = {}
        self.varsFloat = {}
        self.varsEnumStr = {}
        self.data = None
        if data <> None:
            # remove string variables and variables with duplicate names
            newVarList = []
            for var in data.domain.variables + data.domain.getmetas().values():
                if self.varsAll.has_key(var.name):
                    print "Warning: domain contains two variables with the same name: %s; the first one will be used"  % var.name
                    continue
                self.varsAll[var.name] = var
                if var.varType == orange.VarTypes.Continuous:
                    self.varsFloat[var.name] = var
                else:
                    self.varsEnumStr[var.name] = var
                newVarList.append(var)
            domNoMeta = orange.Domain(newVarList, None)
            self.data = orange.ExampleTable(domNoMeta, data)
            self.data.name = data.name
        # update caption title
        self.updateCaptionTitle()
        # fill combos and listBox with variables
        self.fillCmbVars()
        self.fillLbVarOthers()
        self.updateVarAssignment()
        self.initProbes()
        # send data
        if self.commitOnChange:
            self.sendData()
            self.sendProbes()
        qApp.restoreOverrideCursor()       


    def defaultVarAssignmentClick(self):
        """Select default variables in combos based on their names.
        """
        if D1: print "OWNormalize.defaultVarAssignmentClick"
        if self.data:
            qApp.restoreOverrideCursor()
            qApp.setOverrideCursor(QWidget.waitCursor)
            self.setDefaultVarAssignment()
            self.initProbes()
            # send data
            if self.commitOnChange:
                self.sendData()
                self.sendProbes()
            qApp.restoreOverrideCursor()
            
        
    def fillCmbVars(self):
        """ Fills combos with variables; add "<none>":None for 2nd Probe combo.
        """
        if D1: print "OWNormalize.fillCmbVars"
        cmbNamesAllVars = ["cmbVarA", "cmbVarB"]
        cmbNamesFloatVars = ["cmbVarSignalSmpl", "cmbVarSignalRef", "cmbVarBGSmpl", "cmbVarBGRef", "cmbVarBGSmplSD", "cmbVarBGRefSD"]
        if self.data:
            varNamesAllVars = zip(map(lambda x: x.lower(), self.varsAll.keys()), self.varsAll.keys())
            varNamesAllVars.sort()
            varNamesFloatVars = zip(map(lambda x: x.lower(), self.varsFloat.keys()), self.varsFloat.keys())
            varNamesFloatVars.sort()
            # fill VarA & VarB combos
            for cmbName in cmbNamesAllVars:
                self.__dict__[cmbName].clear()
                self.__dict__[cmbName].insertStrList(map(lambda x: x[1], varNamesAllVars))
            # fill intensity combos
            for cmbName in cmbNamesFloatVars:
                self.__dict__[cmbName].clear()
                self.__dict__[cmbName].insertStrList(map(lambda x: x[1], varNamesFloatVars))
        else:
            for cmbName in cmbNamesAllVars + cmbNamesFloatVars:
                self.__dict__[cmbName].clear()            
        # add "<none>"
        if self.varsAll.has_key("<none>"):
            print "Warning: doman consists of discrete variable named '<none>'; this name is reserved and should not be used"
        self.cmbVarB.insertItem("<none>", 0)
        self.cmbVarBGSmplSD.insertItem("<none>", 0)
        self.cmbVarBGRefSD.insertItem("<none>", 0)


    def updateVarAssignment(self):
        """check that the current var names are present in the data
        if not, set to the first variable;
        othervise rewrite the value of the member variable so that the right item in the combo is shown
        """
        if D1 or D2: print "OWNormalize.updateVarAssignment"
        if self.data and len(self.varsAll) > 0:
            varsAllNames = self.varsAll.keys()
            varsAllNames.sort()
            varsEnumStrNames = self.varsEnumStr.keys()
            varsEnumStrNames.sort()
            varsFloatNames = self.varsFloat.keys()
            varsFloatNames.sort()

            # test variable names if exist in the data; rewrite if exist, reset if not
            if not self.varNameA in varsAllNames:
                self._setDefaultVarAssignment_varA(varsEnumStrNames + varsFloatNames)
            else:
                self.varNameA = self.varNameA
            if not self.varNameSignalSmpl in varsAllNames:
                self._setDefaultVarAssignment_varSignalSmpl(varsFloatNames)
            else:
                self.varNameSignalSmpl = self.varNameSignalSmpl
            if not self.varNameSignalRef in varsAllNames:
                self._setDefaultVarAssignment_varSignalRef(varsFloatNames)
            else:
                self.varNameSignalRef = self.varNameSignalRef
            if not self.varNameBGSmpl in varsAllNames:
                self._setDefaultVarAssignment_varBGSmpl(varsFloatNames)
            else:
                self.varNameBGSmpl = self.varNameBGSmpl
            if not self.varNameBGRef in varsAllNames:
                self._setDefaultVarAssignment_varBGRef(varsFloatNames)
            else:
                self.varNameBGRef = self.varNameBGRef
            if not self.varNameB in varsAllNames:
                self._setDefaultVarAssignment_varB(varsEnumStrNames + varsFloatNames)
            else:
                self.varNameB = self.varNameB
            if not self.varNameBGSmplSD in varsAllNames:
                self._setDefaultVarAssignment_varBGSmplSD(varsFloatNames)
            else:
                self.varNameBGSmplSD = self.varNameBGSmplSD
            if not self.varNameBGRefSD in varsAllNames:
                self._setDefaultVarAssignment_varBGRefSD(varsFloatNames)
            else:
                self.varNameBGRefSD = self.varNameBGRefSD
        

    def setDefaultVarAssignment(self):
        """Select default variables in combos based on their names.
        """
        if D1 or D2: print "OWNormalize.setDefaultVarAssignment"
        if self.data and len(self.varsAll) > 0:
            varsAllNames = self.varsAll.keys()
            varsAllNames.sort()
            varsEnumStrNames = self.varsEnumStr.keys()
            varsEnumStrNames.sort()
            varsFloatNames = self.varsFloat.keys()
            varsFloatNames.sort()
            self._setDefaultVarAssignment_varA(varsEnumStrNames + varsFloatNames)
            self._setDefaultVarAssignment_varB(varsEnumStrNames + varsFloatNames)
            self._setDefaultVarAssignment_varSignalSmpl(varsFloatNames)
            self._setDefaultVarAssignment_varSignalRef(varsFloatNames)
            self._setDefaultVarAssignment_varBGSmpl(varsFloatNames)
            self._setDefaultVarAssignment_varBGRef(varsFloatNames)
            self._setDefaultVarAssignment_varBGSmplSD(varsFloatNames)
            self._setDefaultVarAssignment_varBGRefSD(varsFloatNames)


    def _setDefaultVarAssignment_varA(self, dataVarNames):
        self.varNameA = None
        if self.data and len(dataVarNames) > 0:
            # first search variables with equal name, then by substrings, then all variables
            if self.defNameA:
                for name in dataVarNames:
                    if self.defNameA == name:
                        self.varNameA = name
                        break
                if self.varNameA == None:
                    candVarsLenName = []    # list of tuples (len, varName)
                    for name in dataVarNames:
                        if self.defNameA in name:
                            candVarsLenName.append((len(name),name))
                    if len(candVarsLenName) > 0:
                        candVarsLenName.sort()
                        self.varNameA = candVarsLenName[0][1]
            if self.varNameA == None:
                self.varNameA = dataVarNames[0]


    def _setDefaultVarAssignment_varB(self, dataVarNames):
        self.varNameB = "<none>"
        if self.data and len(dataVarNames) > 0:
            # first search variables with equal name, then by substrings, then all variables
            if self.defNameB:
                for name in dataVarNames:
                    if self.defNameB == name:
                        self.varNameB = name
                        break
                if self.varNameB == "<none>":
                    candVarsLenName = []    # list of tuples (len, varName)
                    for name in dataVarNames:
                        if self.defNameB in name:
                            candVarsLenName.append((len(name),name))
                    if len(candVarsLenName) > 0:
                        candVarsLenName.sort()
                        self.varNameB = candVarsLenName[0][1]
            if self.varNameB <> "<none>":
                # select self.varNameB among the other variables (tab Output)
                self.disconnect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
                for idx in range(self.lbVarOthers.count()):
                    if self.varNameB == self.lbVarOthers.item(idx).text():
                        self.lbVarOthers.setSelected(idx, True)
                self.fillVarsOtherSelected()
                self.connect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
            # enable/disable self.boxMergeOtherType
            self.boxMergeOtherType.setEnabled(self.mergeLevel and len(self.varsOtherSelected) > 0)
        

    def _setDefaultVarAssignment_varSignalSmpl(self, dataVarNames):
        self.varNameSignalSmpl = None
        if self.data and len(dataVarNames) > 0:
            candVarsLenName = []    # list of tuples (len, varName)
            for vName in dataVarNames:
                if vName and self.defNameSmpl1 in vName and self.defNameSmpl2 in vName and self.defNameForeground in vName and self.defNameMean in vName:
                    candVarsLenName.append((len(vName),vName))
            if len(candVarsLenName) > 0:
                candVarsLenName.sort()
                self.varNameSignalSmpl = candVarsLenName[0][1]
            else:
                self.varNameSignalSmpl = dataVarNames[0]

    def _setDefaultVarAssignment_varSignalRef(self, dataVarNames):
        self.varNameSignalRef = None
        if self.data and len(dataVarNames) > 0:
            candVarsLenName = []    # list of tuples (len, varName)
            for vName in dataVarNames:
                if vName and self.defNameRef1 in vName and self.defNameRef2 in vName and self.defNameForeground in vName and self.defNameMean in vName:
                    candVarsLenName.append((len(vName),vName))
            if len(candVarsLenName) > 0:
                candVarsLenName.sort()
                self.varNameSignalRef = candVarsLenName[0][1]
            else:
                self.varNameSignalRef = dataVarNames[0]

    def _setDefaultVarAssignment_varBGSmpl(self, dataVarNames):
        self.varNameBGSmpl = None
        if self.data and len(dataVarNames) > 0:
            candVarsLenName = []    # list of tuples (len, varName)
            for vName in dataVarNames:
                if vName and self.defNameSmpl1 in vName and self.defNameSmpl2 in vName and self.defNameBackground in vName and self.defNameMean in vName:
                    candVarsLenName.append((len(vName),vName))
            if len(candVarsLenName) > 0:
                candVarsLenName.sort()
                self.varNameBGSmpl = candVarsLenName[0][1]
            else:
                self.varNameBGSmpl = dataVarNames[0]

    def _setDefaultVarAssignment_varBGRef(self, dataVarNames):
        self.varNameBGRef = None
        if self.data and len(dataVarNames) > 0:
            candVarsLenName = []    # list of tuples (len, varName)
            for vName in dataVarNames:
                if vName and self.defNameRef1 in vName and self.defNameRef2 in vName and self.defNameBackground in vName and self.defNameMean in vName:
                    candVarsLenName.append((len(vName),vName))
            if len(candVarsLenName) > 0:
                candVarsLenName.sort()
                self.varNameBGRef = candVarsLenName[0][1]
            else:
                self.varNameBGRef = dataVarNames[0]

    def _setDefaultVarAssignment_varBGSmplSD(self, dataVarNames):
        self.varNameBGSmplSD = "<none>"
        if self.data and len(dataVarNames) > 0:
            candVarsLenName = []    # list of tuples (len, varName)
            for vName in dataVarNames:
##                if vName and self.defNameSmpl1 in vName and self.defNameSmpl2 in vName and self.defNameBackground in vName and self.defNameSD in vName and self.defNameSmpl and self.defNameBackground and self.defNameSD:
                if vName and self.defNameSmpl1 in vName and self.defNameSmpl2 in vName and self.defNameBackground in vName and self.defNameSD in vName:
                    candVarsLenName.append((len(vName),vName))
            if len(candVarsLenName) > 0:
                candVarsLenName.sort()
                self.varNameBGSmplSD = candVarsLenName[0][1]
        # enable/disable Max. CV slider
        self.boxMaxCV.setEnabled(self.varNameBGSmplSD != "<none>" and self.varNameBGRefSD != "<none>")

    def _setDefaultVarAssignment_varBGRefSD(self, dataVarNames):
        self.varNameBGRefSD = "<none>"
        if self.data and len(dataVarNames) > 0:
            candVarsLenName = []    # list of tuples (len, varName)
            for vName in dataVarNames:
##                if vName and self.defNameRef1 in vName and self.defNameRef2 in vName and self.defNameBackground in vName and self.defNameSD in vName and self.defNameRef and self.defNameBackground and self.defNameSD:
                if vName and self.defNameRef1 in vName and self.defNameRef2 in vName and self.defNameBackground in vName and self.defNameSD in vName:
                    candVarsLenName.append((len(vName),vName))
            if len(candVarsLenName) > 0:
                candVarsLenName.sort()
                self.varNameBGRefSD = candVarsLenName[0][1]
        # enable/disable Max. CV slider
        self.boxMaxCV.setEnabled(self.varNameBGSmplSD != "<none>" and self.varNameBGRefSD != "<none>")
                    

    def fillLbVarOthers(self):
        """Fills listBox with variables not selected by combos
        """
        if D1: print "OWNormalize.fillLbVarOthers"
        self.lbVarOthers.clear()
        if self.data:
            varNamesAllVars = zip(map(lambda x: x.lower(), self.varsAll.keys()), self.varsAll.keys())
            varNamesAllVars.sort()
            self.lbVarOthers.insertStrList(map(lambda x: x[1], varNamesAllVars))
            # select items (self.lbVarOthers <- self.varsOtherSelected)
            self.disconnect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
            for idx in range(self.lbVarOthers.count()):
                self.lbVarOthers.setSelected(idx, self.lbVarOthers.item(idx).text() in self.varsOtherSelected.keys())
            self.connect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)


    def initProbes(self):
        """Init self.probes:
            - reload probe data
            - fill self.tblControls
            - update self.graph
            - update infos
        """
        if D1 or D2 or D6: print "OWNormalize.initProbes"
        if self.data:
            pbStep = 100./len(self.data)
            self.progressBarInit()
            self.probes.initProbes(self.data, self.varNameA, self.varNameB, self.varNameSignalSmpl, self.varNameSignalRef, self.varNameBGSmpl, self.varNameBGRef,
                                   self.varNameBGSmplSD, self.varNameBGRefSD, callback=lambda: self.progressBarAdvance(pbStep))
            self.progressBarFinished()
        else:
            self.probes.clear(False)
        # process external probe data
        self.processDataProbes()
        # fill / update probe table & probe info
        self.fillProbeTable()
        self.setInfoProbes()
        # filter info
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()


    def sendData(self):
        """Compute norm. factors, plot them, normalize data and send out normalized data.
        """
        if D1 or D2 or D6: print "OWNormalize.sendData"
        if self.data:
            # etNum, varListNum: normalized log2 ratios, num. probes, net intensities, A, non-norm M;
            varListNum = [orange.FloatVariable("M normalized")]
            self.progressBarInit()
            l2r = self.probes.getLog2Ratio_norm_masked(self.mergeLevel, self.mergeIntensitiesType, self.progressBarAdvance, False)
            self.progressBarFinished()
            maData = MA.reshape(l2r, (l2r.shape[0], 1))
            # control ratios
            self.progressBarInit()
            varListNum.append(orange.FloatVariable("Control ratio"))
            maData = MA.concatenate([maData, MA.reshape(self.probes.getControlRatios(self.mergeLevel, self.mergeIntensitiesType, self.progressBarAdvance), (maData.shape[0], 1))], 1)
            self.progressBarFinished()
            if self.outNumProbes:
                self.progressBarInit()
                varListNum += [orange.FloatVariable("Num. probes"), orange.FloatVariable("Num. accepted probes")]
                maData = MA.concatenate([maData, self.probes.getNumReplicas_nonFiltered(self.mergeLevel, self.progressBarAdvance)], 1)
                self.progressBarFinished()
            if self.outNetSignal:
                self.progressBarInit()
                varListNum += [orange.FloatVariable("Net intensity (Smpl)"), orange.FloatVariable("Net intensity (Ref)")]
                maData = MA.concatenate([maData, self.probes.getNetIntensity_smpl_ref(self.mergeLevel, self.mergeIntensitiesType, self.progressBarAdvance)], 1)
                self.progressBarFinished()
            if self.outA:
                self.progressBarInit()
                varListNum.append(orange.FloatVariable("A"))
                maData = MA.concatenate([maData, MA.reshape(self.probes.getA_masked(self.mergeLevel, self.mergeIntensitiesType, self.progressBarAdvance), (maData.shape[0], 1))], 1)
                self.progressBarFinished()
            if self.outMRaw:
                self.progressBarInit()
                varListNum.append(orange.FloatVariable("M raw"))
                maData = MA.concatenate([maData, MA.reshape(self.probes.getLog2Ratio_raw_masked(self.mergeLevel, self.mergeIntensitiesType, self.progressBarAdvance, False), (maData.shape[0], 1))], 1)
                self.progressBarFinished()
            if self.outMRaw and self.outMCentered:
                self.progressBarInit()
                varListNum.append(orange.FloatVariable("M raw centered"))
                maData = MA.concatenate([maData, MA.reshape(self.probes.getLog2Ratio_raw_masked(self.mergeLevel, self.mergeIntensitiesType, self.progressBarAdvance, True), (maData.shape[0], 1))], 1)
                self.progressBarFinished()
            if self.outMCentered:
                self.progressBarInit()
                varListNum.append(orange.FloatVariable("M normalized centered"))
                maData = MA.concatenate([maData, MA.reshape(self.probes.getLog2Ratio_norm_masked(self.mergeLevel, self.mergeIntensitiesType, self.progressBarAdvance, True), (maData.shape[0], 1))], 1)
                self.progressBarFinished()
            etNum = chipstat.ma2orng(maData, orange.Domain(varListNum, None))

            # valListList_byAttr: list of lists of values of individual attributes from varListNames + varListOtherCSV; needs to be transposed before converted to ExampleTable
            # varListNames: data table with varA (cloned and renamed to "ID"), varA aliases ("ID alias") and varB
            varAclone = self.data.domain[self.varNameA].clone()
            varAclone.name = "ID"
            varAclone.getValueFrom = lambda e,r: e[self.varNameA]
            varListNames = [varAclone]
            valListList_byAttr = [self.probes.getValsA(self.mergeLevel)]
            if self.outVarAAliases:
                varListNames.append(orange.StringVariable("ID alias"))
                valListList_byAttr.append(self.probes.getValsAAlias(self.mergeLevel))
            if self.varNameB <> "<none>" and self.mergeLevel <> OWNormalize.MergeLevelPerVarA:
                # cannot get unique varB values if self.mergeLevel == OWNormalize.MergeLevelPerVarA
                varListNames.append(self.data.domain[self.varNameB])
                valListList_byAttr.append(self.probes.getValsB(self.mergeLevel))

            # varListOtherCSV: subset of examples/attributes from self.data where non-continuous vars are replaced by String vars named as var.name + " list"
            varListOtherCSV = []
            if len(self.varsOtherSelected) > 0:
                domainOtherVars = orange.Domain(self.varsOtherSelected.values(), None)
                lstLstOtherVars = list(orange.ExampleTable(domainOtherVars, self.data))
                for varIdx, var in enumerate(domainOtherVars):
                    vals = map(lambda x: x[varIdx].native(), lstLstOtherVars)
                    if var.varType == orange.VarTypes.Continuous and self.mergeOtherType <> OWNormalize.MergeOtherTypeConc:
                        valListList_byAttr.append(self.probes._mergeFunc[self.mergeLevel](MA.asarray(vals, MA.Float), Probes.mergeTypes[self.mergeOtherType], lambda x: None))
                        varListOtherCSV.append(var)
                        var.numberOfDecimals = 3    # by default it is set to zero because usually all values are integers; after merging we need higher precision (3 is the default value)
                    else:
                        valListList_byAttr.append(self.probes._concatFunc[self.mergeLevel](vals, removeDupl=self.mergeOtherRemoveDupl))
                        varListOtherCSV.append(orange.StringVariable(var.name + " list"))

            # etNamesOtherCSV: ExampleTable with variables from varListNames + varListOtherCSV, data from valListList_byAttr
            etNamesOtherCSV = orange.ExampleTable(orange.Domain(varListNames+varListOtherCSV, None), Numeric.transpose(Numeric.asarray(valListList_byAttr, Numeric.PyObject)).tolist())

            # final example table: Log2 ratio first, then metas with names, then CSV of other selected attributes
            domOut = orange.Domain([varListNum[0]], None)
            for var in varListNames + varListNum[1:] + varListOtherCSV:
                domOut.addmeta(orange.newmetaid(), var)
            etOut = orange.ExampleTable(domOut, orange.ExampleTable([etNum, etNamesOtherCSV]))
            etOut.name = (self.data.name + " normalized").strip()

        else:
            etOut = None
        self.send("Expression Data", etOut)


    ###################################################################################
    ## PROBE DATA IN / OUT
    ###################################################################################

    def onProbesInput(self, dataProbes):
        """Handles input of probes data.
        """
        if D1 or D2: print "OWNormalize.onProbesInput"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.dataProbes = dataProbes
        # update caption title
        self.updateCaptionTitle()
        # process external probe data
        self.processDataProbes()
        # fill / update probe table & info
        self.fillProbeTable()
        self.setInfoProbes()
        # filters info
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()
        # send data and probes data
        if self.commitOnChange:
            self.sendData()
            self.sendProbes()
        qApp.restoreOverrideCursor()


    def processDataProbes(self):
        """Copy data from orange.ExampleTable to self.probes.
        """
        if D1 or D2 or D6: print "OWNormalize.processDataProbes"
        if self.dataProbes:
            varNames = map(lambda var: var.name, self.dataProbes.domain.variables)
            if "Var A" in varNames and "ColorRGB" in varNames and "Ratio" in varNames and "Symbol" in varNames:
                # list of probe keys
                if self.varNameB <> "<none>" and "Var B" in varNames:
                    pKeysIDsNames = map(lambda e: (str(e["Var A"].native()) + str(e["Var B"].native()), str(e["Var A"].native()), str(e["Var B"].native())), self.dataProbes)
                else:
                    pKeysIDsNames = map(lambda e: (str(e["Var A"].native()), str(e["Var A"].native()), ""), self.dataProbes)
                # get color, symbol and ratio of probes
                for e, (pKey, pID, pName) in zip(self.dataProbes, pKeysIDsNames):
                    c = str(e["ColorRGB"])
                    try:
                        color = QColor(int(c[0:2],16), int(c[2:4],16), int(c[4:6],16), QColor.Rgb)
                    except ValueError:
                        # this test is here because Excel changes "000000" to "0"
                        if c=="0":
                            color = QColor(0,0,0)
                        else:
                            color = ProbeSet.NoColor
                    try:
                        symbol = int(e["Symbol"].native())
                    except TypeError:
                        symbol = QwtSymbol.None
                    ratio = str(e["Ratio"])
                    # read optional name for varA and varB values
                    try:
                        valAAlias = str(e["Var A alias"].native())
                    except TypeError:
                        valAAlias = None
                    ###########################################################
                    # if pKey exists
                    probe = self.probes.get(pKey)
                    if probe:
                        self.probes.setMarker(probe, symbol, color, refresh=False)
                        self.probes.setRatio(probe, ratio, recalc=False, refresh=False)
                        if valAAlias:
                            probe.valAAlias = valAAlias
                    # if pKey does not exist
                    else:
                        for pk in self.probes.getKeysFromValAValB(pID, pName):
                            probe = self.probes.get(pk)
                            self.probes.setMarker(probe, symbol, color, refresh=False)
                            self.probes.setRatio(probe, ratio, recalc=False, refresh=False)
                            if valAAlias:
                                probe.valAAlias = valAAlias
                if len(self.dataProbes) > 0:
                    self.probes.replotProbeCurves(refresh=False)
                    self.probes.updateReplotNormCurves(refresh=True)
            else:
                print "Warning: probe data must consist of attributes 'Var A', 'Ratio', 'ColorRGB' and 'Symbol'; optional attributes: 'Var B' and 'Var A alias'"


    def sendProbes(self):
        """Sends out example table with currently used probes.
        """
        if D1 or D2 or D6: print "OWNormalize.sendProbes"
        if self.data:
            vals = []
            pKeys = self.probes.keys()
            pKeys.sort()
            for pKey in pKeys:
                probe = self.probes[pKey]
                vals.append([probe.valA, probe.valAAlias, probe.valB, probe.ratioExpr, probe.getColorRgbStr(), probe.symbol])
            # create domain, clone and rename variables to "Var A" and "Var B"
            if self.varNameB <> "<none>":
                varB = self.varsAll[self.varNameB].clone()
                varB.name = "Var B"
            else:
                varB = orange.StringVariable("Var B")
            varAAlias = orange.StringVariable("Var A alias")
            domain = orange.Domain([self.varsAll[self.varNameA].clone(), varAAlias, varB, orange.StringVariable("Ratio"), orange.StringVariable("ColorRGB"), orange.FloatVariable("Symbol")], None)
            domain[0].name = "Var A"
            # copy values into example table
            et = orange.ExampleTable(domain, vals)
            # remove 'Var B' if selv.varNameB == "<none>"
            if self.varNameB == "<none>":
                et = orange.ExampleTable(orange.Domain(et.domain[0:2] + et.domain[3:], None), et)
            if self.dataProbes:
                et.name = self.dataProbes.name
            else:
                et.name = "Probe Data"
            self.send("Probe Data", et)
        else:
            self.send("Probe Data", None)


    ###################################################################################
    ## SELECTION OF VARIABLES
    ###################################################################################

    def varABChange(self):
        """Refresh listbox containing other variables and refill self.tblControls.
        """
        if D1: print "OWNormalize.varABChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.initProbes()
        # send data
        if self.commitOnChange:
            self.sendData()
            self.sendProbes()
        qApp.restoreOverrideCursor()


    def varFBChange(self):
        """Refresh listbox containing other variables (tab Out);
        enables/disables Max. CV slider (tab Filter).
        """
        if D1 or D6: print "OWNormalize.varFBChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        # enable/disable Max. CV slider
        self.boxMaxCV.setEnabled(self.varNameBGSmplSD != "<none>" and self.varNameBGRefSD != "<none>")
        # update data
        self.probes.updateProbeData(self.data, self.varNameSignalSmpl, self.varNameSignalRef, self.varNameBGSmpl, self.varNameBGRef, self.varNameBGSmplSD, self.varNameBGRefSD, recalc=True)
        # update info
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def varOthersChange(self):
        """Updates list of selected other vars (lbVarOthers -> self.varsOtherSelected);
        enables merge options for other variables.
        """
        self.fillVarsOtherSelected()
        self.boxMergeOtherType.setEnabled(self.mergeLevel and len(self.varsOtherSelected) > 0)
        if self.commitOnChange:
            self.sendData()


    def fillVarsOtherSelected(self):        
        if D1: print "OWNormalize.varOtherChange"
        self.varsOtherSelected = {}
        if self.data:
            for i in range(0, self.lbVarOthers.count()):
                if self.lbVarOthers.isSelected(i):
                    varName = str(self.lbVarOthers.item(i).text())
                    self.varsOtherSelected[varName] = self.data.domain[varName]


    ###################################################################################
    ## PROBE TABLE (self.tblControls)
    ###################################################################################

    def fillProbeTable(self):
        # init self.tblControls
        if D1 or D2 or D6: print "OWNormalize.fillProbeTable"
        if self.probes:
            # table row items and their number, header labels, len of sorting keys where spaces are added in front
            allProbes = self.probes.values()
            numProbes = len(allProbes)
            self.tblControls.setNumRows(numProbes)
            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcVarA, "ID")
            sortKeyLen = int(math.log(numProbes, 10))+1
            pbStep = 100./len(self.probes)
            self.progressBarInit()
            # generate sorting keys: varA, varAAlias
            valAIdx = zip(map(lambda pr: (pr.valA, pr.valAAlias, pr.valB), allProbes), range(numProbes))
            valAIdx.sort()
            valARank = dict(zip([x[1] for x in valAIdx], range(numProbes)))
            valAAliasIdx = zip(map(lambda pr: (pr.valAAlias, pr.valA, pr.valB), allProbes), range(numProbes))
            valAAliasIdx.sort()
            valAAliasRank = dict(zip([x[1] for x in valAAliasIdx], range(numProbes)))
            if self.varNameB <> "<none>":
                # generate sorting keys: varB)
                valBIdx = zip(map(lambda pr: (pr.valB, pr.valA, pr.valAAlias), allProbes), range(numProbes))
                valBIdx.sort()
                valBRank = dict(zip([x[1] for x in valBIdx], range(numProbes)))
                # fill rows
                for row, (key,pr) in enumerate(self.probes.items()):
                    # store the table row index of the probe
                    pr.tblRowIdx = row
                    numProbesAll = pr.getNumProbes()
                    numProbesAcc = self.probes.getNumProbes_nonFiltered(pr)
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcNPAll, str(numProbesAll), editType=QTableItem.Never, sortingKey="%05d%05d"%(numProbesAll,numProbesAcc))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcNPAccepted, str(numProbesAcc), editType=QTableItem.Never, sortingKey="%05d%05d"%(numProbesAcc,numProbesAll))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcVarA, str(pr.valA), editType=QTableItem.Never, sortingKey=self.sortingKey(valARank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcVarAAlias, str(pr.valAAlias), editType=QTableItem.OnTyping, sortingKey=self.sortingKey(valAAliasRank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcVarB, str(pr.valB), editType=QTableItem.Never, sortingKey=self.sortingKey(valBRank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcPKey, key)
                    self.fillProbeTableItemMarkerRatio(row, pr)
                    self.progressBarAdvance(pbStep)
                self.tblControls.horizontalHeader().setLabel(OWNormalize.tcVarB, self.varNameB)
            else:
                # fill rows
                for row, (key,pr) in enumerate(self.probes.items()):
                    # store the table row index of the probe
                    pr.tblRowIdx = row
                    numProbesAll = pr.getNumProbes()
                    numProbesAcc = self.probes.getNumProbes_nonFiltered(pr)
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcNPAll, str(numProbesAll), editType=QTableItem.Never, sortingKey="%05d%05d"%(numProbesAll,numProbesAcc))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcNPAccepted, str(numProbesAcc), editType=QTableItem.Never, sortingKey="%05d%05d"%(numProbesAcc,numProbesAll))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcVarA, str(pr.valA), editType=QTableItem.Never, sortingKey=self.sortingKey(valARank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcVarAAlias, str(pr.valAAlias), editType=QTableItem.OnTyping, sortingKey=self.sortingKey(valAAliasRank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcPKey, key)
                    self.fillProbeTableItemMarkerRatio(row, pr)
                    self.progressBarAdvance(pbStep)
            self.progressBarFinished()
        else:
            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcVarA, "ID")
            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcVarB, "Alt.var")
            self.tblControls.setNumRows(0)
        self.adjustProbeTableColumns()


    def adjustProbeTableColumns(self):
        """shows/hides appropriate columns, adjusts their widths
        """
        # show/hide columns with aliases and varB values and adjust their widths
        lastColumn = OWNormalize.tcVarA
        if self.displayVarAAliases:
            self.tblControls.showColumn(OWNormalize.tcVarAAlias)
            self.tblControls.adjustColumn(OWNormalize.tcVarAAlias)
            lastColumn = OWNormalize.tcVarAAlias
        else:
            self.tblControls.hideColumn(OWNormalize.tcVarAAlias)
            self.tblControls.setColumnWidth(OWNormalize.tcVarAAlias, 0)

        if self.varNameB <> "<none>":
            self.tblControls.showColumn(OWNormalize.tcVarB)
            self.tblControls.adjustColumn(OWNormalize.tcVarB)
            lastColumn = OWNormalize.tcVarB
        else:
            self.tblControls.hideColumn(OWNormalize.tcVarB)
            self.tblControls.setColumnWidth(OWNormalize.tcVarB, 0)
        # adjust other columns' width
        self.tblControls.adjustColumn(OWNormalize.tcRatio)
        self.tblControls.adjustColumn(OWNormalize.tcNPAll)
        self.tblControls.adjustColumn(OWNormalize.tcNPAccepted)
        self.tblControls.adjustColumn(OWNormalize.tcVarA)
        # calculate the width of the last column
        colWidthSum = self.tblControls.columnWidth(OWNormalize.tcMarker) + self.tblControls.columnWidth(OWNormalize.tcRatio) + self.tblControls.columnWidth(OWNormalize.tcNPAll) + self.tblControls.columnWidth(OWNormalize.tcNPAccepted) + self.tblControls.columnWidth(OWNormalize.tcVarA) + self.tblControls.columnWidth(OWNormalize.tcVarAAlias) + self.tblControls.columnWidth(OWNormalize.tcVarB)
        lastColWidth = max(self.tblControls.columnWidth(lastColumn), self.tblControls.visibleWidth() - colWidthSum + self.tblControls.columnWidth(lastColumn))
        self.tblControls.setColumnWidth(lastColumn, lastColWidth)


    def fillProbeTableItemMarkerRatio(self, row, probe):
        """Updates a given row of self.tblControls: marker & ratio.
        """
        pxm = self.probes.getSymbolPixmap(probe, self.tblControls.cellGeometry(row, OWNormalize.tcMarker))
        ratioSortKey = self.probes.getRatioSortingKey(probe)
        markerSortKey = probe.getMarkerSortKey()
        OWGUI.tableItem(self.tblControls, row, OWNormalize.tcRatio, self.probes.getRatioStr(probe), editType=QTableItem.OnTyping, sortingKey=ratioSortKey+markerSortKey) #, background=QColor(160,160,160))
        OWGUI.tableItem(self.tblControls, row, OWNormalize.tcMarker, "", editType=QTableItem.Never, sortingKey=markerSortKey+ratioSortKey, pixmap=pxm)#, background=QColor(160,160,160))


    def sortingKey(self, val, len):
        """Returns a string with leading spaces followed by str(val), whose length is at least len.
        """
        s = "%" + str(int(len)) + "." + str(int(len/2)) + "f"
        return s % val


    def updateProbeTableNumAcceptedProbes(self):
        # updates column tcNPAccepted after filters change
        if D1 or D2 or D6: print "OWNormalize.updateProbeTableNumAcceptedProbes"
        if self.probes:
            # table row items and their number, header labels, len of sorting keys where spaces are added in front
            for row in range(self.tblControls.numRows()):
                pr = self.probes.get(str(self.tblControls.text(row, OWNormalize.tcPKey)))
                numProbesAll = pr.getNumProbes()
                numProbesAcc = self.probes.getNumProbes_nonFiltered(pr)
                OWGUI.tableItem(self.tblControls, row, OWNormalize.tcNPAccepted, str(numProbesAcc), editType=QTableItem.Never, sortingKey="%05d%05d"%(numProbesAcc,numProbesAll))
 

    def tblControlsHHeaderClicked(self, col):
        """Sorts the table by column col; sets probe.tblRowIdx accordingly
        """
        if D1: print "OWNormalize.tblControlsSort"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        if col == self.sortby-1:
            self.sortby = - self.sortby
        else:
            self.sortby = col+1
        self.tblControls.sortColumn(col, self.sortby>=0, TRUE)
        self.tblControls.horizontalHeader().setSortIndicator(col, self.sortby<0)
        qApp.restoreOverrideCursor()
        # update probe.tblRowIdx
        for row in range(self.tblControls.numRows()):
            pKey = str(self.tblControls.item(row, OWNormalize.tcPKey).text())
            self.probes[pKey].tblRowIdx = row


    def tblControlsValueChange(self, row, col):
        """Handles direct changes to self.tblControls;
        if col == tcRatio: updates self.internContRatios and sends out expression and probe data;
        if col == tcVarAAlias: updates probe.valAAlias and sends out probe data;
        """
        if D1 or D3: print "OWNormalize.tblControlsValueChange"
        if col == OWNormalize.tcRatio:
            pKey = str(self.tblControls.item(row, OWNormalize.tcPKey).text())
            ratio = str(self.tblControls.item(row, OWNormalize.tcRatio).text())
            if ratio <> self.probes[pKey].ratioExpr:
                # set ratio
                probe = self.probes[pKey]
                newRatio = self.probes.setRatio(probe, ratio, recalc=True)
                # update table & info
                self.fillProbeTableItemMarkerRatio(row, probe)
                self.setInfoProbes()
                self.setInfoFilterMaxCV()
                self.setInfoFilterMinRatio()
                self.setInfoFilterMaxFGInt()
                self.setInfoFilterMaxBGInt()
                if self.commitOnChange:
                    self.sendData()
                    self.sendProbes()
        if col == OWNormalize.tcVarAAlias:
            pKey = str(self.tblControls.item(row, OWNormalize.tcPKey).text())
            alias = str(self.tblControls.item(row, OWNormalize.tcVarAAlias).text())
            if alias <> self.probes[pKey].valAAlias:
                self.probes[pKey].valAAlias = alias
                if self.commitOnChange:
                    if self.outVarAAliases:
                        self.sendData()
                    self.sendProbes()
                
                


    def tblControlsCurrentChanged(self, row, col):
        """Handles changes of currently selected cell in self.tblControls;
        activate the current probe.
        """
        if D1: print "OWNormalize.tblControlsCurrentChanged"
        if row >= 0 and col >= 0:
            pKey = str(self.tblControls.item(row, OWNormalize.tcPKey).text())
            self.probes.setCurveActive(pKey, True, refresh=True)


    def tblControlsSelectionChanged(self):
        if D1: print "OWNormalize.tblControlsSelectionChanged"
        self.activateSelectedProbes()


    def tblControlsDoubleClicked(self, row, col, button, mousePos):
        """Adjust values of self.ratioStr, self.btnProbeColor and self.cmbProbeSymbol.
        """
        if D1: print "OWNormalize.tblControlsSelectionChanged"
        if row>=0 and col>=0:
            pKey = str(self.tblControls.item(row, OWNormalize.tcPKey).text())
            probe = self.probes[pKey]
            # update GUI: ratio lineEdit, probeSymbol, probeColor
            self.probeSymbolIdx = probe.symbol
            self.ratioStr = probe.ratioExpr
            if probe.color <> ProbeSet.NoColor:
                self.probeColor = probe.color
                # update color of button
                self.btnProbeColor.pixmap().fill(self.probeColor)
                self.btnProbeColor.repaint()


    ###################################################################################
    ## ASSIGNMENT OF RATIOS / COLORS / SYMBOLS
    ###################################################################################

    def btnSelectControlsClick(self):
        """Select probes where ID contains self.controlName.
        """
        if D1: print "OWNormalize.btnSelectControlsClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.tblControls.setCurrentCell(-1,-1)
        numSelectionsOld = self.tblControls.numSelections()
        for idx in range(self.tblControls.numRows()):
            if string.lower(self.controlName) in string.lower(self.tblControls.item(idx, OWNormalize.tcVarA).text()):
                sel = QTableSelection()
                sel.init(idx,OWNormalize.tcVarA)
                sel.expandTo(idx,OWNormalize.tcVarA)
                self.tblControls.addSelection(sel)
        if self.tblControls.numSelections() <> numSelectionsOld:
            self.activateSelectedProbes()
        qApp.restoreOverrideCursor()


    def btnSelectControlsAllClick(self):
        """Clears all selections and selects all probes.
        """
        if D1: print "OWNormalize.btnSelectControlsAllClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.tblControls.clearSelection(False)
        sel = QTableSelection()
        sel.init(0,OWNormalize.tcMarker)
        sel.expandTo(self.tblControls.numRows()-1,OWNormalize.tcVarB)
        self.tblControls.addSelection(sel)
        self.activateSelectedProbes()
        qApp.restoreOverrideCursor()


    def btnUnselectControlsAllClick(self):
        """Clears all selections and selects all probes.
        """
        if D1: print "OWNormalize.btnUnselectControlsAllClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.tblControls.setCurrentCell(-1,-1)
        self.tblControls.clearSelection(True)
        qApp.restoreOverrideCursor()


    def leRatioReturnPressed(self):
        """If new ratio is entered and return pressed, selected controls get updated.
        """
        if D1: print "OWNormalize.leRatioReturnPressed"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.updateSelectedProbes(self.ratioStr, self.probeColor, self.probeSymbolIdx, True, False, False)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()
        if self.commitOnChange:
            self.sendData()
            self.sendProbes()
        qApp.restoreOverrideCursor()


    def probeColorClick(self):
        """Show color-selection window, update button color and color of the symbols of the selected probes.
        """
        if D1: print "OWNormalize.probeColorCick"
        probeColor = QColorDialog.getColor(self.probeColor, self)
        if probeColor.isValid():
            self.probeColor = probeColor
        self.btnProbeColor.pixmap().fill(self.probeColor)
        self.btnProbeColor.repaint()

        
    def btnSetProbesClick(self):
        """Sets ratio, color and symbol for the selected controls.
        """
        if D1: print "OWNormalize.btnSetProbesClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.updateSelectedProbes(self.ratioStr, self.probeColor, self.probeSymbolIdx, True, True, True)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()
        if self.commitOnChange:
            self.sendData()
            self.sendProbes()
        qApp.restoreOverrideCursor()


    def btnClearProbesClick(self):
        """Clears ratios for the selected controls
        """
        if D1: print "OWNormalize.btnClearProbesClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.updateSelectedProbes("", ProbeSet.NoColor, QwtSymbol.None, True, True, True)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()
        if self.commitOnChange:
            self.sendData()
            self.sendProbes()
        qApp.restoreOverrideCursor()


    def updateSelectedProbes(self, ratioStr, color, symbol, updateRatio, updateColor, updateSymbol):
        """Updates selected probes with a given ratio string, color and symbol;
        updates only those values where the correcponding update value == True;
        calls fillProbeTableItemMarkerRatio(row, probe), replots/removes curves.
        """
        if D1 or D2: print "OWNormalize.updateSelectedProbes"
        numSelections = self.tblControls.numSelections()
        for selNum in range(self.tblControls.numSelections()):
            sel = self.tblControls.selection(selNum)
            for row in range(sel.topRow(), sel.bottomRow()+1):
                # update probe
                pKey = str(self.tblControls.item(row, OWNormalize.tcPKey).text())
                probe = self.probes[pKey]
                if updateRatio:
                    self.probes.setRatio(probe, ratioStr, recalc=False, refresh=False)
                if updateColor or updateSymbol:
                    if updateColor:
                        newColor = color
                    else:
                        newColor = probe.color
                    if updateSymbol:
                        newSymbol = symbol
                    else:
                        newSymbol = probe.symbol
                    self.probes.setMarker(probe, newSymbol, newColor, refresh=False)
                # update table
                self.fillProbeTableItemMarkerRatio(row, probe)
        if self.tblControls.numSelections() > 0:
            self.probes.updateReplotNormCurves(refresh=True)


    def activateSelectedProbes(self):
        """Activate currently selected probes;
        """
        if D1 or D2 or D6: print "OWNormalize.activateSelectedProbes"
        activePKeys = {}    # dict pKey:pKey (to remove duplicates)
        for selNum in range(self.tblControls.numSelections()):
            sel = self.tblControls.selection(selNum)
            for row in range(sel.topRow(), sel.bottomRow()+1):
                item = self.tblControls.item(row, OWNormalize.tcPKey)
                if item:
                    aPKey = str(item.text())
                    activePKeys[aPKey] = aPKey
        self.probes.setCurveActiveList(activePKeys.keys(), refresh=True)


    def setInfoProbes(self):
        if D4: print "OWNormalize.setInfoProbes"
        if self.probes:
            self.lblProbeInfo.setText("Probes in total:\t%d control,\t%d other.\nAccepted:\t%d control,\t%d other.\nPlotted:\t\t%d control,\t%d other." %
                                      (self.probes.getNumProbesControls(), self.probes.getNumProbesOthers(),
                                       self.probes.getNumProbesControls_nonFiltered(), self.probes.getNumProbesOthers_nonFiltered(),
                                       self.probes.getNumProbesControls_nonFiltered_plotted(), self.probes.getNumProbesOthers_nonFiltered_plotted()))
        else:
            self.lblProbeInfo.setText("No data on input.\n\n")


    ###################################################################################
    ## GRAPH
    ###################################################################################

    def setGraphAxes(self, axes=None):
        """According to selected scaling sets up axis labels and scales;
        axis: 0: vertical left, 1: vertical right, 2: horizontal bottom, 3: horizontal top
        """
        if D1: print "OWNormalize.setGraphAxes"
        titles = {False: ["Centered ratio"]*2 + ["Average intensity"]*2, True: ["M: centered log2 ratio"]*2 + ["A: log2 average intensity"]*2}
        useLog = [self.logAxisY]*2 + [self.logAxisX]*2
        if axes==None: axes = [0,2]
        for axis in axes:
            self.graph.setAxisTitle(axis, titles[useLog[axis]][axis])
    

    def onLegendClicked(self, key):
        """Change active probe curve
        """
        if D1: print "OWNormalize.onLegendClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        pKey = self.graph.curve(key).key
        if pKey <> None:
            self.probes.switchCurveActive(pKey, refresh=True)
        qApp.restoreOverrideCursor()


    def onMouseMoved(self, e):
        """Find closest curve (if close enough), activate corresponding probe, print tooltip.
        """
        if D1 or D6: print "OWNormalize.onMouseMoved"
        (curveKey, distPoints, x, y, pointKey) = self.graph.closestCurve(e.x(), e.y())
        curve = self.graph.curve(curveKey)
        # if we have a curve and the curve has a "key" (member variable specific for our curves)
        if curve and curve.__dict__.has_key("key"):
            # if we are close enough to the curve
            if distPoints<=self.markerSize/2:
                # activata probe and/or normalization curve (if not already active)
                self.probes.setCurveActiveList([curve.key], refresh=True)
                # show tooltip
                self.probes.showDataTooltip(curve.key, x, y)
            else:
                # deactivate all curves
                self.probes.setCurveActiveList([], refresh=True)


    def onMousePressed(self, e):
        """If left button is pressed, the active control is made current in self.tblControls;
        first, OWGraph.onMousePressed(e) is executed (by default), followed by this code.
        """
        if D1: print "OWNormalize.onMousePressed"
        actCrvKeys = self.probes.getActiveCurveKeys()
        if e.button() == Qt.LeftButton and len(actCrvKeys) > 0:
            # clear selections from tblControls (without fireing events)
            self.disconnect(self.tblControls , SIGNAL('selectionChanged()'), self.tblControlsSelectionChanged)
            self.tblControls.clearSelection(True)
            selectedRows = []
            for crvKey in actCrvKeys:
                probe = self.probes.get(crvKey)
                if probe:
                    selectedRows.append(probe.tblRowIdx)
                    # select appropriate table cells
                    sel = QTableSelection()
                    sel.init(probe.tblRowIdx,OWNormalize.tcVarA)
                    sel.expandTo(probe.tblRowIdx,OWNormalize.tcVarA)
                    self.tblControls.addSelection(sel)
            # set current cell
            selectedRows.sort()
            if selectedRows:
                self.tblControls.setCurrentCell(selectedRows[0], OWNormalize.tcMarker)
            self.connect(self.tblControls , SIGNAL('selectionChanged()'), self.tblControlsSelectionChanged)

    ###################################################################################
    ## SETTINGS
    ###################################################################################

    def settingsSubstrBGChange(self):
        if D1: print "OWNormalize.settingsSubstrBGChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.setSubtrBG(self.subtrBG)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxFGInt()
        self.setInfoFilterMaxBGInt()
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def settingsFilterMaxCVChange(self):
        """Handles changes of filter settings, which affects graph curves and ouput data.
        """
        if D1: print "OWNormalize.settingsFiltersChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.setProbeFilters()
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.updateProbeTableNumAcceptedProbes()
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def settingsFilterMinIntRatioChange(self):
        """Handles changes of filter settings, which affects graph curves and ouput data.
        """
        if D1: print "OWNormalize.settingsFiltersChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.setProbeFilters()
        self.setInfoProbes()
        self.setInfoFilterMinRatio()
        self.updateProbeTableNumAcceptedProbes()
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def settingsFilterMaxFGIntChange(self):
        """Handles changes of filter settings, which affects graph curves and ouput data.
        """
        if D1: print "OWNormalize.settingsFilterMaxFGIntChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.setProbeFilters()
        self.setInfoProbes()
        self.setInfoFilterMaxFGInt()
        self.updateProbeTableNumAcceptedProbes()
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def settingsFilterMaxBGIntChange(self):
        """Handles changes of filter settings, which affects graph curves and ouput data.
        """
        if D1: print "OWNormalize.settingsFilterMaxBGIntChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.setProbeFilters()
        self.setInfoProbes()
        self.setInfoFilterMaxBGInt()
        self.updateProbeTableNumAcceptedProbes()
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def setProbeFilters(self):      
        if D1: print "OWNormalize.setProbeFilters"
        if self.useCV:
            maxCV = self.maxCV
        else:
            maxCV = Probes.bigVal
        if self.useMinIntensity:  
            minIntensityRatio = self.minIntensityRatio
        else:
            minIntensityRatio = 0
        if self.useMaxFGIntensity:
            maxFGIntensity = self.maxFGIntensity
        else:
            maxFGIntensity = Probes.bigVal
        if self.useMaxBGIntensity:
            maxBGIntensity = self.maxBGIntensity
        else:
            maxBGIntensity = Probes.bigVal
        self.probes.setFilterParameters(maxCV, minIntensityRatio, maxFGIntensity, maxBGIntensity)


    def setInfoFilterMaxCV(self):
        if D4: print "OWNormalize.setInfoFilterMaxCV"
        if self.probes:
            if self.probes.getNumProbesControls() > 0:
                ratioControls = 100. * self.probes.getNumFilteredControlsMaxCV() / self.probes.getNumProbesControls()
            else:
                ratioControls = 0
            if self.probes.getNumProbesOthers() > 0:
                ratioOthers = 100. * self.probes.getNumFilteredOthersMaxCV() / self.probes.getNumProbesOthers()
            else:
                ratioOthers = 0
            self.lblInfoFilterMaxCV.setText("%d (%.2f%s) control probes removed.\n%d (%.2f%s) other probes removed." % (self.probes.getNumFilteredControlsMaxCV(), ratioControls, "%", self.probes.getNumFilteredOthersMaxCV(), ratioOthers, "%"))
        else:
            self.lblInfoFilterMaxCV.setText("No data on input.\n")


    def setInfoFilterMinRatio(self):
        if D4: print "OWNormalize.setInfoFilterMinIntRatio"
        if self.probes:
            if self.probes.getNumProbesControls() > 0:
                ratioControls = 100. * self.probes.getNumFilteredControlsMinRatio() / self.probes.getNumProbesControls()
            else:
                ratioControls = 0
            if self.probes.getNumProbesOthers() > 0:
                ratioOthers = 100. * self.probes.getNumFilteredOthersMinRatio() / self.probes.getNumProbesOthers()
            else:
                ratioOthers = 0
            self.lblInfoFilterMinIntRatio.setText("%d (%.2f%s) control probes removed.\n%d (%.2f%s) other probes removed." % (self.probes.getNumFilteredControlsMinRatio(), ratioControls, "%", self.probes.getNumFilteredOthersMinRatio(), ratioOthers, "%"))
        else:
            self.lblInfoFilterMinIntRatio.setText("No data on input.\n")


    def setInfoFilterMaxFGInt(self):
        if D4: print "OWNormalize.setInfoFilterMaxFGInt"
        if self.probes:
            if self.probes.getNumProbesControls() > 0:
                ratioControls = 100. * self.probes.getNumFilteredControlsMaxFGInt() / self.probes.getNumProbesControls()
            else:
                ratioControls = 0
            if self.probes.getNumProbesOthers() > 0:
                ratioOthers = 100. * self.probes.getNumFilteredOthersMaxFGInt() / self.probes.getNumProbesOthers()
            else:
                ratioOthers = 0
            self.lblInfoFilterMaxFGInt.setText("%d (%.2f%s) control probes removed.\n%d (%.2f%s) other probes removed." % (self.probes.getNumFilteredControlsMaxFGInt(), ratioControls, "%", self.probes.getNumFilteredOthersMaxFGInt(), ratioOthers, "%"))
        else:
            self.lblInfoFilterMaxFGInt.setText("No data on input.\n")


    def setInfoFilterMaxBGInt(self):
        if D4: print "OWNormalize.setInfoFilterMaxBGInt"
        if self.probes:
            if self.probes.getNumProbesControls() > 0:
                ratioControls = 100. * self.probes.getNumFilteredControlsMaxBGInt() / self.probes.getNumProbesControls()
            else:
                ratioControls = 0
            if self.probes.getNumProbesOthers() > 0:
                ratioOthers = 100. * self.probes.getNumFilteredOthersMaxBGInt() / self.probes.getNumProbesOthers()
            else:
                ratioOthers = 0
            self.lblInfoFilterMaxBGInt.setText("%d (%.2f%s) control probes removed.\n%d (%.2f%s) other probes removed." % (self.probes.getNumFilteredControlsMaxBGInt(), ratioControls, "%", self.probes.getNumFilteredOthersMaxBGInt(), ratioOthers, "%"))
        else:
            self.lblInfoFilterMaxBGInt.setText("No data on input.\n")


    def settingsArrayTypeChange(self):
        """Handles changes of array type (0: genome array, 1: boutique array)
        """
##TODO!!!        
        pass


    def settingsNormalizationChange(self):
        """Handles changes of normalization type, which affects normalization curve and output data.
        """
        if D1: print "OWNormalize.settingsNormalizationChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.setNormalizationParameters(self.normRange, self.minNumControlProbes, self.approxFunction, self.loessWindow, self.loessWeight)
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def filtersAllChange(self):
        """handles changes caused by Set Default Values button"""
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        chngF = False
        # subtract background
        if self.subtrBG <> self._def_subtrBG:
            self.subtrBG = self._def_subtrBG
            self.probes.setSubtrBG(self.subtrBG)
            chngF = True
        # CV            
        if self.useCV <> self._def_useCV:
            self.useCV = self._def_useCV
            chngF = True
        if self.maxCV <> self._def_maxCV:
            self.maxCV = self._def_maxCV
            chngF = True
        # min. intensity ratio
        if self.useMinIntensity <> self._def_useMinIntensity:
            self.useMinIntensity = self._def_useMinIntensity
            chngF = True
        if self.minIntensityRatio <> self._def_minIntensityRatio:
            self.minIntensityRatio = self._def_minIntensityRatio
            chngF = True
        # max FG intensity
        if self.useMaxFGIntensity <> self._def_useMaxFGIntensity:
            self.useMaxFGIntensity = self._def_useMaxFGIntensity
            chngF = True
        if self.maxFGIntensity <> self._def_maxFGIntensity:
            self.maxFGIntensity = self._def_maxFGIntensity
            chngF = True
        # max BG intensity
        if self.useMaxBGIntensity <> self._def_useMaxBGIntensity:
            self.useMaxBGIntensity = self._def_useMaxBGIntensity
            chngF = True
        if self.maxBGIntensity <> self._def_maxBGIntensity:
            self.maxBGIntensity = self._def_maxBGIntensity
            chngF = True
        # refresh
        if chngF:
            self.setProbeFilters()
            self.setInfoProbes()
            self.setInfoFilterMaxCV()
            self.setInfoFilterMinRatio()
            self.setInfoFilterMaxFGInt()
            self.setInfoFilterMaxBGInt()
            if self.commitOnChange:
                self.sendData()
        qApp.restoreOverrideCursor()


    def normalizationAllChange(self):
        """handles changes caused by Set Default Values button"""
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        chngN = False
        if self.normRange <> self._def_normRange:
            self.normRange = self._def_normRange
            chngN = True
        if self.minNumControlProbes <> self._def_minNumControlProbes:
            self.minNumControlProbes = self._def_minNumControlProbes
            chngN = True
        if self.approxFunction <> self._def_approxFunction:
            self.approxFunction = self._def_approxFunction
            chngN = True
        if self.loessWindow <> self._def_loessWindow:
            self.loessWindow = self._def_loessWindow
            chngN = True
        if self.loessWeight <> self._def_loessWeight:
            self.loessWeight = self._def_loessWeight
            chngN = True
        # refresh
        if chngN:
            self.probes.setNormalizationParameters(self.normRange, self.minNumControlProbes, self.approxFunction, self.loessWindow, self.loessWeight)
            if self.commitOnChange:
                self.sendData()
        qApp.restoreOverrideCursor()


    ###################################################################################
    ## OUTPUT SETTINGS
    ###################################################################################


    def settingsOutputReplicasChange(self):
        """Handles changes of replicas settings, which affect output data;
        enables / disables some output options.
        """
        if D1: print "OWNormalize.settingsOutputReplicasChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.rbgMergeIntensitiesType.setEnabled(self.mergeLevel)
        self.boxMergeOtherType.setEnabled(self.mergeLevel and len(self.varsOtherSelected) > 0)
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()


    def settingsOutputOtherChange(self):        
        """Handles changes of output settings, which affect output data only if other variables are selected.
        """
        if D1: print "OWNormalize.settingsOutputChange"
        if len(self.varsOtherSelected) > 0:
            qApp.restoreOverrideCursor()
            qApp.setOverrideCursor(QWidget.waitCursor)
            if self.commitOnChange:
                self.sendData()
            qApp.restoreOverrideCursor()


    def settingsOutputChange(self): 
        """Handles changes of output settings; send out data.
        """
        if D1: print "OWNormalize.settingsOutputChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        if self.commitOnChange:
            self.sendData()
        qApp.restoreOverrideCursor()
        
    def settingsGraphAxisChange(self, axis):
        """Handles changes of graph axis settings; replot axes and curves.
        """
        if D1: print "OWNormalize.settingsGraphAxisChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.setPlotParameters(self.logAxisX, self.logAxisY, self.markerSize, OWNormalize.normCurveStyles[self.normCurveStyleIdx], refresh=True)
        self.setGraphAxes([axis])
        qApp.restoreOverrideCursor()


    def settingsGraphChange(self):
        """Handles changes of graph settings; replot curves.
        """
        if D1: print "OWNormalize.settingsGraphChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.setPlotParameters(self.logAxisX, self.logAxisY, self.markerSize, OWNormalize.normCurveStyles[self.normCurveStyleIdx], refresh=True)
        qApp.restoreOverrideCursor()


    def settingsShowLegendChange(self):
        """Enables/disables legend.
        """
        if D1: print "OWNormalize.settingsShowLegendChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.graph.enableGraphLegend(self.showLegend)
        qApp.restoreOverrideCursor()


    def settingsProbeTrackingChange(self):
        if D1: print "OWNormalize.settingsProbeTrackingChange"
        signal1 = SIGNAL("plotMouseMoved(const QMouseEvent &)")
        signal2 = SIGNAL('plotMousePressed(const QMouseEvent&)')
        if self.tracking:
            self.connect(self.graph, signal1, self.onMouseMoved)
            self.connect(self.graph, signal2, self.onMousePressed)
        else:
            self.disconnect(self.graph, signal1, self.onMouseMoved)
            self.disconnect(self.graph, signal2, self.onMousePressed)
            

    def commitChange(self):
        self.btnCommit.setEnabled(not self.commitOnChange)
        if self.commitOnChange:
            self.commitClicked()


    def commitClicked(self):
        """Handles Commit click, sends out examples and probes.
        """
        if D1: print "OWNormalize.commitClicked"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.sendData()
        self.sendProbes()
        qApp.restoreOverrideCursor()



class QwtPlotKeyCurve(QwtPlotCurve):
    """QwtPlotCurve with additional member variable: key.
    """

    def __init__(self, parent, name, key):
        self.key = key
        QwtPlotCurve.__init__(self, parent, name)



class OWGraphMA(OWGraph):
    """OWGraph for MA plot with curves of type QwtPlotKeyCurve;
    additionally handles zooming and mouse clicks.
    """

    def insertCurve(self, title, key, xAxis=QwtPlot.xBottom, yAxis=QwtPlot.yLeft):
        curve = QwtPlotKeyCurve(self, title, key)
        curve.setAxis(xAxis, yAxis)
        return OWGraph.insertCurve(self, curve)

    
    def zoomOut(self):
        """Overridden in order to fix for autoscaling when there is no zoom.
        """
        if len(self.zoomStack):
            (xmin, xmax, ymin, ymax) = self.zoomStack.pop()
            if len(self.zoomStack):
                self.setAxisScale(QwtPlot.xBottom, xmin, xmax)
                self.setAxisScale(QwtPlot.yLeft, ymin, ymax)
            else:
                self.setAxisAutoScale(QwtPlot.xBottom)
                self.setAxisAutoScale(QwtPlot.xTop)
                self.setAxisAutoScale(QwtPlot.yLeft)
                self.setAxisAutoScale(QwtPlot.yRight)
            self.replot()
            return 1
        return 0


    def onMousePressed(self, e):
        self.mouseCurrentlyPressed = 1
        self.mouseCurrentButton = e.button()
        self.xpos = e.x()
        self.ypos = e.y()

        # ####
        # ZOOM
        if e.button() == Qt.LeftButton and self.state == ZOOMING:
            self.tempSelectionCurve = SelectionCurve(self, pen = Qt.DashLine)
            self.zoomKey = OWGraph.insertCurve(self, self.tempSelectionCurve)

        # ####
        # SELECT RECTANGLE
        elif e.button() == Qt.LeftButton and self.state == SELECT_RECTANGLE:
            self.tempSelectionCurve = SelectionCurve(self)
            key = OWGraph.insertCurve(self, self.tempSelectionCurve)
            self.selectionCurveKeyList.append(key)

        # ####
        # SELECT POLYGON
        elif e.button() == Qt.LeftButton and self.state == SELECT_POLYGON:
            if self.tempSelectionCurve == None:
                self.tempSelectionCurve = SelectionCurve(self)
                key = OWGraph.insertCurve(self, self.tempSelectionCurve)
                self.selectionCurveKeyList.append(key)
                self.tempSelectionCurve.addPoint(self.invTransform(QwtPlot.xBottom, self.xpos), self.invTransform(QwtPlot.yLeft, self.ypos))
            self.tempSelectionCurve.addPoint(self.invTransform(QwtPlot.xBottom, self.xpos), self.invTransform(QwtPlot.yLeft, self.ypos))

            if self.tempSelectionCurve.closed():    # did we intersect an existing line. if yes then close the curve and finish appending lines
                self.tempSelectionCurve = None
                self.replot()
                if self.autoSendSelectionCallback: self.autoSendSelectionCallback() # do we want to send new selection
        self.event(e)


class NormCurveDataDict(dict):
    """Dictionary of NormCurveData objects.
    """
    class NormCurveData:
        """Class for storing data about normalization curve.
        """
        def __init__(self, key, computInd, normInd, probeList):
            self._key = key                     # a key for NormCurveDataDict, equal to QwtPlotKeyCurve.key
            self.computInd = computInd          # data indices for computing the normalization curve
            self.normInd = normInd              # data indices for plotting tick marks; they indicate the probes which are normalized by this curve
            self.probeList = probeList          # list of ProbeSet
            self.curveList = []                 # list of long key returned by OWGraphMA.insertCurve, i.e. QwtPlot.insertCurve

        def getKey(self):
            return self._key

        def extendComputInd(self, computInd):
            tmpDict = dict(zip(self.computInd, self.computInd))
            tmpDict.update(dict(zip(computInd, computInd)))
            self.computInd = tmpDict.keys()

        def extendNormInd(self, normInd):
            tmpDict = dict(zip(self.normInd, self.normInd))
            tmpDict.update(dict(zip(normInd, normInd)))
            self.normInd = tmpDict.keys()

        def extendProbeList(self, probeList):
            tmpDict = dict(zip(self.probeList, self.probeList))
            tmpDict.update(dict(zip(probeList, probeList)))
            self.probeList = tmpDict.keys()


    def __setitem__(self, key, ncd):
        if ncd.__class__.__name__ <> "NormCurveData":
            raise ValueError, "instance of class NormCurveData expected, got %s" % ncd.__class__.__name__
        if ncd.getKey() <> key:
            raise ValueError, "key (%s) and ncd.key (%s) do not match" % (str(key), str(ncd.getKey()))
        dict.__setitem__(self, key, ncd)

    def add(self, key, computInd, normInd, probeList):
        if self.has_key(key):
            self.__getitem__(key).extendComputInd(computInd)
            self.__getitem__(key).extendNormInd(normInd)
            self.__getitem__(key).extendProbeList(probeList)
        else:
            self.__setitem__(key, NormCurveDataDict.NormCurveData(key, computInd, normInd, probeList))


class ProbeSet:
    """A set of equivalent probes, their intensities, and their expected ratio, display symbol & color;
    curve represents normalization factors where color and symbol is used to plot the curve on the graph.
    """
    PenWidthActiveProbe = 3
    PenWidthActiveCurve = 3
    PenWidthInactiveProbe = 1
    PenWidthInactiveCurve = 1
    PenWidths = {False:PenWidthInactiveCurve,
                 True: PenWidthActiveCurve}
    NoColor = QColor(0,0,0)
    
    def __init__(self, valA, valB, pKey):
        if D1: print "ProbeSet.__init__"
        self.valA = valA            # string
        self.valAAlias = valA       # string that is shown in the probe table
        self.valB = valB            # string
        self.pKey = pKey            # ID / ID+name
        # ratio
        self.ratioExpr = ""         # string from which ratio is evaluated
        # marker
        self.color = ProbeSet.NoColor       # QColor
        self.symbol = QwtSymbol.None        # int (QwtSymbol.Style)
        # table, gaph
        self.tblRowIdx = None       # row index in OWNormalize.tblControls
        self.curve = None           # curve in OWNormalize.graph
        # data indices
        self._dataIndices = []


    ############################################
    # DATA
    ############################################

    def getDataIndices(self):
        if D1: print "ProbeSet.getDataIndices"
        return self._dataIndices


    def addProbeIdx(self, dataIdx):
        if D1: print "ProbeSet.addProbeIdx"
        self._dataIndices.append(dataIdx)


    def getNumProbes(self):
        return len(self._dataIndices)


    ############################################
    # TABLE
    ############################################

    def generateSymbolPixmap(self, rect):
        """Input: QRect instance; output: QPixmap instance.
        """
        if D1: print "ProbeSet.generateSymbolPixmap"
        # init pixmap for table item
        symbolPixmap = QPixmap(rect.width(),rect.height())
        symbolPixmap.fill(QColor(255,255,255))
        if self.symbol <> QwtSymbol.None:
            painter = QPainter(symbolPixmap)
            symbol = QwtSymbol(self.symbol, QBrush(self.color, QBrush.SolidPattern), QPen(QColor(0,0,0),1), QSize(8,8))
            symbol.draw(painter, QPoint(rect.width()/2,rect.height()/2))
            painter.end()
        return symbolPixmap


    def getColorRgbStr(self):
        """returns hexadecimal color string, e.g. "00ff10
        """
        return string.replace("%2s%2s%2s" % (hex(self.color.red())[2:], hex(self.color.green())[2:], hex(self.color.blue())[2:]), " ", "0")


    def getMarkerSortKey(self):
        """string for sorting merkers and storing pixmaps to dict
        """
        ch = self.color.hsv()
        return "%02d%03d%03d%03d" % (self.symbol, ch[0], ch[1], ch[2])



class Probes(dict):
    """Dictionary of ProbeSet items; key: pKey (probe varA + optionally probe varB, item: ProbeSet
    Technical issues:
        recalc: recalculate l2r-s
        refresh: self.graph.replot()
    TODO:
        - setFilterParameters: divide into three functions, one for each parameter
    """

    mergeTypes = {0:MA.average, 1:NumExtn.medianMA}
    bigVal = 1e20
    midVal = 1e10

    NormRangeGlobal = 0
    NormRangeLocal = 1
    NormRangeCombined = 2

    # global normalization curve key
    # local norm. curve names equal varB values
    NormCurveNameGlobal = "NormCurveNameGlobal"

    def __init__(self, graph, subtrBG, logAxisX, logAxisY, markerSize, normCurveStyle):
        if D1 or D2 or D6: print "Probes.__init__"
        self._active = {}   # currently active probeSet and normalization curves; key & value: curve.key
        self.graph = graph
        # plot parameters
        self.subtrBG = subtrBG
        self.logAxisX = logAxisX
        self.logAxisY = logAxisY
        self.markerSize = markerSize
        self.normCurveStyle = normCurveStyle
        # var names
        self.varNameA = None
        self.varNameB = "<none>"
        # var A values, var A aliases, var B values
        self._valAList = []         # consecutive list of varA values (replicated)
        self._valBList = []         # consecutive list of varB values (replicated)
        # data: Numeric arrays
        self.__sigSmpl = MA.zeros((0,), MA.Float)
        self.__sigRef = MA.zeros((0,), MA.Float)
        self.__bgSmpl = MA.zeros((0,), MA.Float)
        self.__bgRef = MA.zeros((0,), MA.Float)
        self.__bgSmplSD = MA.zeros((0,), MA.Float)
        self.__bgRefSD = MA.zeros((0,), MA.Float)
        # net intensity function dict (indexed by self.subtrBG)
        self.__netSmpl_masked_func = {0: lambda cond: self._sigSmpl_masked(cond), 1: lambda cond: self._sigSmpl_masked(cond) - self._bgSmpl_masked(cond)}
        self.__netRef_masked_func =  {0: lambda cond: self._sigRef_masked(cond),  1: lambda cond: self._sigRef_masked(cond) -  self._bgRef_masked(cond) }
        # ratios
        self.__ratio = None     # MA array with control ratios and masked non-controls
        # Numeric array: 0: OK, 1: filtered out
        self.__filterMaxCV = None
        self.__filterMinRatio = None
        self.__filterMaxFGInt = None
        self.__filterMaxBGInt = None
##        self._control = Numeric.zeros((0,), Numeric.Int)
        self.__plotted = Numeric.zeros((0,), Numeric.Int)
        # normalization functions
        self._approxFunctionDict = {0:self._getNormCurveMedian, 1:self._getNormCurveLinReg, 2:self._getNormCurveLoess}
        self._approxFunction = self._approxFunctionDict[2]
        # normalization curves
        self._ncdd = NormCurveDataDict()
        # default parameters
        self._normRange = None  # 0: NormRangeGlobal, 1: NormRangeLocal (per var B values); 2: NormRangeCombined
        self._minNumControlProbes = 0
        self.maxCV = Probes.bigVal
        self.minIntensityRatio = 0
        self.maxFGIntensity = Probes.bigVal
        self.maxBGIntensity = Probes.bigVal
        self.loessWindow = 60
        self.loessWeight = 0
        # for normalization per varB values
        self._valB2ind = {}     # {valB1:[indices1], valB2:[indices2]}; key: probe valB; value: list of data indices
        self._valB2probes = {}  # {valB1:[probes1],  valB2:[probes2]};  key: probe valB; value: list of ProbeSets
        self._valA2ind = {}       # {valA1:[indices1], valA2:[indices2]}; key: probe valA; value: list of data indices
        # normalized log2 ratio
        self._Mnorm = MA.zeros((0,), MA.Float)
        # different levels of merging
        self._mergeFunc = {OWNormalize.MergeLevelNone:self.__mergeReplicasNone,
                           OWNormalize.MergeLevelPerVarsAB:self.__mergeReplicasPerVarsAB,
                           OWNormalize.MergeLevelPerVarA:self.__mergeReplicasPerVarA}
        self._concatFunc = {OWNormalize.MergeLevelNone:self.__concatReplicasNone,
                            OWNormalize.MergeLevelPerVarsAB:self.__concatReplicasPerVarsAB,
                            OWNormalize.MergeLevelPerVarA:self.__concatReplicasPerVarA}
        # dictionary of marker pixmaps
        self.__markerPixmapDict = {}    # key: ProbeSet.getMarkerSortKey(), item: pixmap


    def clear(self, refresh=True):
        if D1 or D2 or D6: print "Probes.clear"
        self._clearNormCurves(False)
        self._removeAllProbeCurves(refresh)
        dict.clear(self)
        self._active = {}
        self._valAList = []
        self._valBList = []
        self.__sigSmpl = MA.zeros((0,), MA.Float)
        self.__sigRef = MA.zeros((0,), MA.Float)
        self.__bgSmpl = MA.zeros((0,), MA.Float)
        self.__bgRef = MA.zeros((0,), MA.Float)
        self.__bgSmplSD = MA.zeros((0,), MA.Float)
        self.__bgRefSD = MA.zeros((0,), MA.Float)
        # ratio
        self.__ratio = None
        # Numeric array: 0: OK, 1: filtered out
        self.__filterMaxCV = None
        self.__filterMinRatio = None
        self.__filterMaxFGInt = None
        self.__filterMaxBGInt = None
##        self._control = Numeric.zeros((0,), Numeric.Int)
        self.__plotted = Numeric.zeros((0,), Numeric.Int)
        # for normalization per varB values
        self._valB2ind = {}
        self._valB2probes = {}
        self._valA2ind = {}
        # normalized log2 ratio
        self._Mnorm = MA.zeros((0,), MA.Float)
        # dictionary of marker pixmaps
        self.__markerPixmapDict = {}


    def setFilterParameters(self, maxCV, minIntensityRatio, maxFGIntensity, maxBGIntensity):
        if D1 or D2 or D6: print "Probes.setFilterParameters"
        if maxCV <> self.maxCV or minIntensityRatio <> self.minIntensityRatio or maxFGIntensity <> self.maxFGIntensity or maxBGIntensity <> self.maxBGIntensity:
            self.maxCV = maxCV                        
            self.minIntensityRatio = minIntensityRatio
            self.maxFGIntensity = maxFGIntensity
            self.maxBGIntensity = maxBGIntensity
            self.__filterMaxCV = None
            self.__filterMinRatio = None
            self.__filterMaxFGInt = None
            self.__filterMaxBGInt = None
            self.replotProbeCurves(False)
            self.updateReplotNormCurves(True)


    def setNormalizationParameters(self, normRange, minNumControlProbes, approxFunction, loessWindow, loessWeight):
        """approxFunction: 0: median, 1: LR, 2: LOESS
        """
        if D1 or D2 or D6: print "Probes.setNormalizationParameters"
        change = False
        if self._minNumControlProbes <> minNumControlProbes or normRange <> self._normRange:
            self._normRange = normRange
            self._minNumControlProbes = minNumControlProbes
            change = True
        if self._approxFunctionDict[approxFunction] <> self._approxFunction or loessWindow <> self.loessWindow or loessWeight <> self.loessWeight:
            self._normRange = normRange
            self._approxFunction = self._approxFunctionDict[approxFunction]
            self.loessWindow = loessWindow
            self.loessWeight = loessWeight
            change = True
        if change:
            # refresh=True is needed in order to avoid error upon moving a mouse over an old curve while self._ncdd is being updated after changing normRange
            self._clearNormCurves(True)
            self.updateReplotNormCurves(True)


    def setSubtrBG(self, subtrBG):
        if D1 or D2 or D6: print "Probes.setSubtrBG"
        if self.subtrBG <> subtrBG:
            self.subtrBG = subtrBG
            self.replotProbeCurves(False)
            self.updateReplotNormCurves(True)
        

    def setPlotParameters(self, logAxisX, logAxisY, markerSize, normCurveStyle, refresh=True):
        if D1 or D2 or D6: print "Probes.setPlotParameters"
        if self.logAxisX <> logAxisX or  self.logAxisY <> logAxisY or self.markerSize <> markerSize or self.normCurveStyle <> normCurveStyle:
            self.logAxisX = logAxisX
            self.logAxisY = logAxisY
            self.markerSize = markerSize
            self.normCurveStyle = normCurveStyle
            self.replotProbeCurves(False)
            self.updateReplotNormCurves(refresh)


    def initProbes(self, data, varNameA, varNameB, varNameSignalSmpl, varNameSignalRef, varNameBGSmpl, varNameBGRef, varNameBGSmplSD, varNameBGRefSD, callback):
        """Input: orange.ExampleTable, orange.Variable.names;
        stores self.varNameA and self.varNameB.
        """
        if D1 or D2 or D6: print "Probes.initProbes"
        self.clear(refresh=True)
        self.varNameA = varNameA
        self.varNameB = varNameB
##        self._control = Numeric.zeros((len(data),), Numeric.Int)
        self.__plotted = Numeric.zeros((len(data),), Numeric.Int)
        self.__ratio = MA.ones((len(data),), MA.Float) * MA.masked
        # update data and probe data indices
        self.updateProbeData(data, varNameSignalSmpl, varNameSignalRef, varNameBGSmpl, varNameBGRef, varNameBGSmplSD, varNameBGRefSD, recalc=False)
        # set varPKey and add varNameB (if necessary)
        domainNew = orange.Domain(data.domain)
        # add var 'pKey(varNameAvarNameB)' with string values from varNameA (+ varNameB)
        varPKey = orange.StringVariable("pKey(" + self.varNameA + self.varNameB + ")")
        domainNew.addmeta(orange.newmetaid(), varPKey)
        if self.varNameB <> "<none>":
            varPKey.getValueFrom = lambda example, returnWhat: orange.Value(varPKey, str(example[self.varNameA].native()) + str(example[self.varNameB].native()))
        else:
            varPKey.getValueFrom = lambda example, returnWhat: orange.Value(varPKey, str(example[self.varNameA].native()))
            # add varNameB with unknown values
            domainNew.addmeta(orange.newmetaid(), orange.StringVariable("<none>"))
        data = orange.ExampleTable(domainNew, data)
        # create / update ProbeSets, update self._valAList, self._valBList, self._valB2ind, self._valB2probes, self._valA2ind
        for eIdx, e in enumerate(data):
            pKey = e[varPKey].native()
            valA = e[varNameA].native()
            valB = e[varNameB].native()
            # update self._valAList, self._valBList
            self._valAList.append(valA)
            self._valBList.append(valB)
            # add/update ProbeSet
            if dict.has_key(self, pKey):
                ps = dict.get(self, pKey)
            else:
                ps = ProbeSet(valA, valB, pKey)
                dict.__setitem__(self, pKey, ps)
            ps.addProbeIdx(eIdx)
            # store probes/data indices for normalization per valB values
            if self._valB2ind.has_key(valB):
                self._valB2ind[valB].append(eIdx)
                self._valB2probes[valB].append(ps)
            else:
                self._valB2ind[valB] = [eIdx]
                self._valB2probes[valB] = [ps]
            if self._valA2ind.has_key(valA):
                self._valA2ind[valA].append(eIdx)
            else:
                self._valA2ind[valA] = [eIdx]
            # progressbar callback
            callback()
        self.replotProbeCurves(False)
        self.updateReplotNormCurves(True)
        

    def updateProbeData(self, data, varNameSignalSmpl, varNameSignalRef, varNameBGSmpl, varNameBGRef, varNameBGSmplSD, varNameBGRefSD, recalc=True):
        """Update signal and background of the selected probes, construct new filter.
        """
        if D1 or D2 or D6: print "Probes.updateProbeData"
        if data:
            # keep only the variables that you need
            domVarList = [data.domain[varNameSignalSmpl], data.domain[varNameSignalRef], data.domain[varNameBGSmpl], data.domain[varNameBGRef]]
            if varNameBGSmplSD <> "<none>":
                domVarList.append(data.domain[varNameBGSmplSD])
            if varNameBGRefSD <> "<none>":
                domVarList.append(data.domain[varNameBGRefSD])
            dom = orange.Domain(domVarList, None)
            data = orange.ExampleTable(dom, data)
            # convert to MA
            dataMA = data.toMA("a")[0]
            self.__sigSmpl = dataMA[:,data.domain.index(varNameSignalSmpl)]
            self.__sigRef = dataMA[:,data.domain.index(varNameSignalRef)]
            self.__bgSmpl = dataMA[:,data.domain.index(varNameBGSmpl)]
            self.__bgRef = dataMA[:,data.domain.index(varNameBGRef)]
            if varNameBGSmplSD <> "<none>":
                self.__bgSmplSD = dataMA[:,data.domain.index(varNameBGSmplSD)]
            else:
                self.__bgSmplSD = None
            if varNameBGRefSD <> "<none>":
                self.__bgRefSD = dataMA[:,data.domain.index(varNameBGRefSD)]
            else:
                self.__bgRefSD = None
            self.__filterMaxCV = None
            self.__filterMinRatio = None
            self.__filterMaxFGInt = None
            self.__filterMaxBGInt = None
            if recalc:
                self.replotProbeCurves(False)
                self.updateReplotNormCurves(True)

    ###########################################################################################################

    def getSymbolPixmap(self, probe, rect):
        markerSortKey = probe.getMarkerSortKey()
        pxm = self.__markerPixmapDict.get(markerSortKey)
        if not pxm:
            pxm = probe.generateSymbolPixmap(rect)
            self.__markerPixmapDict[markerSortKey] = pxm
        return pxm



    def getKeysFromValAValB(self, valASub, valBSub):
        """Return a list of probe keys which contain a given part of valA and valB (optionally);
        """
        if D1: print "Probes.getKeysFromValAValB"
        probesValAValB = {}   # probes where both valA and valB matches
        probesValA = {}       # probes where only valA matches
        valASub = str(valASub)
        valBSub = str(valBSub)
        for pKey, probe in self.items():
            if valASub and (valASub in str(probe.valA)):
                probesValA[pKey] = probe
                if valBSub and (valBSub in str(probe.valB)):
                    probesValAValB[pKey] = probe
        if len(probesValAValB) > 0:
            return probesValAValB.keys()
        else:
            return probesValA.keys()


    ############################################
    # NUMBER OF PROBES
    ############################################

    def _isControlArr(self):
        """returns array of len(data) where 1 for controls and 0 for others
        """
        return Numeric.asarray(Numeric.logical_not(MA.getmaskarray(self.__ratio)), Numeric.Int)

    def _isOtherArr(self):
        """returns array of len(data) where 0 for controls and 1 for others
        """
        return Numeric.asarray(MA.getmaskarray(self.__ratio), Numeric.Int)


    def getNumReplicas_nonFiltered(self, mergeLevel, callback):
        """Returns (..., 2) Numeric array where rows represent different probes and columns:
            0: number of all probes
            1: number of non-filtered probes
        """
        na = Numeric.transpose(Numeric.asarray([Numeric.ones(self.getFilter().shape), Numeric.logical_not(self.getFilter())]))
        return self._mergeFunc[mergeLevel](na, Numeric.add.reduce, callback)

    def getNumFilteredControlsMaxCV(self):
        if type(self.__filterMaxCV) == types.NoneType:
            self._setFilterMaxCV()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxCV, self._isControlArr()))

    def getNumFilteredOthersMaxCV(self):
        if type(self.__filterMaxCV) == types.NoneType:
            self._setFilterMaxCV()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxCV, self._isOtherArr()))

    def getNumFilteredControlsMinRatio(self):
        if type(self.__filterMinRatio) == types.NoneType:
            self._setFilterMinRatio()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMinRatio, self._isControlArr()))

    def getNumFilteredOthersMinRatio(self):
        if type(self.__filterMinRatio) == types.NoneType:
            self._setFilterMinRatio()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMinRatio, self._isOtherArr()))

    def getNumFilteredControlsMaxFGInt(self):
        if type(self.__filterMaxFGInt) == types.NoneType:
            self._setFilterMaxFGInt()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxFGInt, self._isControlArr()))
    
    def getNumFilteredOthersMaxFGInt(self):
        if type(self.__filterMaxFGInt) == types.NoneType:
            self._setFilterMaxFGInt()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxFGInt, self._isOtherArr()))

    def getNumFilteredControlsMaxBGInt(self):
        if type(self.__filterMaxBGInt) == types.NoneType:
            self._setFilterMaxBGInt()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxBGInt, self._isControlArr()))
    
    def getNumFilteredOthersMaxBGInt(self):
        if type(self.__filterMaxBGInt) == types.NoneType:
            self._setFilterMaxBGInt()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxBGInt, self._isOtherArr()))


    def getNumProbes(self):
        return self.__ratio.shape[0]

    def getNumProbesControls(self):
        return Numeric.add.reduce(Numeric.greater(self._isControlArr(), 0))

    def getNumProbesOthers(self):
        return Numeric.add.reduce(self._isOtherArr())


    def getNumProbesControls_filtered(self):
        return Numeric.add.reduce(Numeric.logical_and(self._isControlArr(), self.getFilter()))

    def getNumProbesControls_nonFiltered(self):
        return Numeric.add.reduce(Numeric.logical_and(self._isControlArr(), Numeric.logical_not(self.getFilter())))

    def getNumProbesOthers_filtered(self):
        return Numeric.add.reduce(Numeric.logical_and(self._isOtherArr(), self.getFilter()))

    def getNumProbesOthers_nonFiltered(self):
        return Numeric.add.reduce(Numeric.logical_and(self._isOtherArr(), Numeric.logical_not(self.getFilter())))


    def getNumProbes_nonFiltered(self, probe):
        return Numeric.add.reduce(Numeric.take(Numeric.logical_not(self.getFilter()), probe.getDataIndices()))


    def getNumProbes_nonFiltered_plotted(self):
        return Numeric.add.reduce(Numeric.logical_and(self.__plotted, Numeric.logical_not(self.getFilter())))

    def getNumProbesControls_nonFiltered_plotted(self):
        return Numeric.add.reduce(Numeric.logical_and(self._isControlArr(), Numeric.logical_and(self.__plotted, Numeric.logical_not(self.getFilter()))))

    def getNumProbesOthers_nonFiltered_plotted(self):
        return Numeric.add.reduce(Numeric.logical_and(self._isOtherArr(), Numeric.logical_and(self.__plotted, Numeric.logical_not(self.getFilter()))))


    def getNumProbesControls_indexed(self, indices):
        return Numeric.add.reduce(Numeric.logical_and(self._isControlArr(), NumExtn.indices2condition(indices, self.__sigSmpl.shape[0])))

    def getNumProbesOthers_indexed(self, indices):
        return Numeric.add.reduce(Numeric.logical_and(self._isOtherArr(), NumExtn.indices2condition(indices, self.__sigSmpl.shape[0])))

    def getNumProbesControls_nonFiltered_indexed(self, indices):
        if D6: print "getNumProbesControls_nonFiltered_indexed", Numeric.add.reduce(self._isControlArr()), Numeric.add.reduce(self.getFilter())
        return Numeric.add.reduce(Numeric.logical_and(self._isControlArr(), Numeric.logical_and(NumExtn.indices2condition(indices, self.__sigSmpl.shape[0]), Numeric.logical_not(self.getFilter()))))

    def getNumProbesOthers_nonFiltered_indexed(self, indices):
        return Numeric.add.reduce(Numeric.logical_and(self._isOtherArr(), Numeric.logical_and(NumExtn.indices2condition(indices, self.__sigSmpl.shape[0]), Numeric.logical_not(self.getFilter()))))


    ############################################
    # MARKER, RATIO
    ############################################

    def setMarker(self, probe, symbol, color, refresh=True):
        """Set symbol or color to None to remove probe set from the graph;
        """
        if D1 or D2 or D3: print "Probes.setMarker"
        probe.symbol = symbol
        probe.color = color
        self._replotProbeCurve(probe, refresh)


    def setRatio(self, probe, ratioExpr, recalc, refresh=True):
        """Sets self.__ratio and dict[pKey].ratioExpr;
        ratioExpr should be None in order not to use probeSet as a control.
        """
        if D1 or D2 or D3: print "Probes.setRatio"
        try:
            ratio = float(eval(ratioExpr))
            if ratio <= 0:
                ratio = 1
                newRatioExpr = ""
            else:
                newRatioExpr = ratioExpr
        except:
            ratio = 1
            newRatioExpr = ""
        # if self.__ratio different from ratio
        if probe.ratioExpr <> newRatioExpr:
            probe.ratioExpr = newRatioExpr
            MA.put(self.__ratio, probe.getDataIndices(), ratio)
##            # if we succeed in setting the ratio: put probe among controls
##            Numeric.put(self._control, probe.getDataIndices(), probe.ratioExpr <> "")
            self._replotProbeCurve(probe, refresh and not recalc)
            # recalc and replot norm. curves
            if recalc:                    
                self.updateReplotNormCurves(refresh)
        return ratio


    def getRatioStr(self, probe):
        if D1: print "Probes.getRatioStr"
        if probe.ratioExpr and len(probe.getDataIndices())>0:
            return "%.2f" % self.__ratio.filled(1)[probe.getDataIndices()[0]]
        else:
            return ""


    def getRatioSortingKey(self, probe):
        """Returns a string with leading spaces followed by str(val), whose length is at least len."""
        if probe.ratioExpr and len(probe.getDataIndices())>0:
            return "%15.7f" % self.__ratio.filled(1)[probe.getDataIndices()[0]]
        else:
            return ""


    def getActiveCurveKeys(self):
        if D1: print "Probes.getActiveCurveKeys"
        return self._active.keys()


    ############################################
    # PROBE CURVES
    ############################################

    def _removeProbeCurve(self, probe, refresh=True):
        if D1 or D2 or D3: print "Probes._removeProbeCurve"
        if probe and probe.curve:
            Numeric.put(self.__plotted, probe.getDataIndices(), 0)
            self.graph.removeCurve(probe.curve)
            probe.curve = None
        if refresh: self.graph.replot()


    def _removeAllProbeCurves(self, refresh=True):
        if D1 or D2 or D6: print "Probes._removeAllProbeCurves"
        for probe in self.values():
            if probe.curve:
                self.graph.removeCurve(probe.curve)
                probe.curve = None
        self.__plotted *= 0
        if refresh and len(self)>0: self.graph.replot()


    def _replotProbeCurve(self, probe, refresh=True):
        if D1 or D2 or D3: print "Probes._replotProbeCurve"
        change = False
        if probe.curve:
            self._removeProbeCurve(probe, False)
            change = True
        if probe.symbol <> QwtSymbol.None:
            probe.curve = self.graph.insertCurve(probe.valA, probe.pKey)
            Numeric.put(self.__plotted, probe.getDataIndices(), 1)
            M,A = self.getMA(probe.pKey, True)
            self.graph.setCurveData(probe.curve, A, M)
            self.graph.setCurveStyle(probe.curve, QwtCurve.NoCurve)
            self._setProbeCurveSymbol(probe, False)
            change = True
        if change and refresh: self.graph.replot()


    def _setProbeCurveSymbol(self, probe, refresh=True):
        """sets graph marker symbol
        """
        if probe.curve:
            if self._active.has_key(probe.pKey):
                pen = QPen(QColor(0,0,0),ProbeSet.PenWidthActiveProbe)
            else:
                pen = QPen(QColor(0,0,0),ProbeSet.PenWidthInactiveProbe)
            qSymbol = QwtSymbol(probe.symbol, QBrush(probe.color, QBrush.SolidPattern), pen, QSize(self.markerSize,self.markerSize))
            self.graph.setCurveSymbol(probe.curve, qSymbol)
            if refresh: self.graph.replot()


    def replotProbeCurves(self, refresh=True):
        """iterate all probes, remove their curves (if exist) and replot them (if symbol <> None)
        """
        if D1 or D2 or D6: print "Probes.replotProbeCurves"
        for probe in self.values():
            self._replotProbeCurve(probe, False)
        if refresh: self.graph.replot()


    ############################################
    # GRAPH: NORM. CURVES; NORMALIZED LOG2 RATIOS
    ############################################

    def updateReplotNormCurves(self, refresh=True):
        """updates self._Mnorm and replots norm curves
        """
        if D1 or D2 or D6: print "Probes.updateReplotNormCurves"
        self._clearNormCurves(False)
        self._Mnorm = MA.zeros(self.__sigSmpl.shape, MA.Float) * MA.masked  # log2ratio
        if self.__sigSmpl:
            condColors = [QColor(0,0,0), QColor(0,0,255), QColor(0,255,0)]
            # fill self._ncdd (NormCurveDataDict)
            if self._normRange == Probes.NormRangeGlobal:
                self._ncdd.add(Probes.NormCurveNameGlobal, range(len(self.__sigSmpl)), range(len(self.__sigSmpl)), self.values())
            elif self._normRange == Probes.NormRangeLocal:
                for valB, ind in self._valB2ind.items():
                    self._ncdd.add(valB, ind, ind, self._valB2probes[valB])
            elif self._normRange == Probes.NormRangeCombined:
                for valB, ind in self._valB2ind.items():
                    if self.getNumProbesControls_nonFiltered_indexed(ind) < self._minNumControlProbes:
                        self._ncdd.add(Probes.NormCurveNameGlobal, range(len(self.__sigSmpl)), ind, self._valB2probes[valB])
                    else:
                        self._ncdd.add(valB, ind, ind, self._valB2probes[valB])
            else:
                raise ValueError, "Unknown Probes._normRange: %s" % str(self._normRange)
            # cycle all norm curve data, normalize and plot norm. curves
            for ncKey, ncData in self._ncdd.items():
                Ac, Mc, An_masked, Mn_masked = self._getNormData_ActrlMctrl_AM_masked(ncData.computInd)
                # condition for plotting ticks
                condTicks = NumExtn.indices2condition(ncData.normInd, An_masked.shape[0])
                if len(Ac) <= 0:
                    continue
                minAc = min(Ac)
                maxAc = max(Ac)
                # plot norm. curve
                notMask_AnMn = MA.logical_not(MA.getmaskarray(An_masked+Mn_masked))
                # condList: [interpolated part, extrapolated lower part, extrapolated upper part of norm. curve]
                condList =  [MA.logical_and(MA.logical_and(MA.greater_equal(An_masked, minAc), MA.less_equal(An_masked, maxAc)), notMask_AnMn),
                             MA.logical_and(MA.less_equal(An_masked, minAc), notMask_AnMn),
                             MA.logical_and(MA.greater_equal(An_masked, maxAc), notMask_AnMn)]
                # add curveLists to self._ncdd items
                for condIdx, cond in enumerate(condList):
                    if MA.add.reduce(MA.asarray(cond, MA.Float)) > 1:
                        # plot normalization curve
                        normCurve = self.graph.insertCurve("Norm. curve %i: %s" % (condIdx, str(ncKey)), ncKey)
                        ncData.curveList.append(normCurve)
                        Aplot = Numeric.asarray(MA.compress(cond, An_masked))
                        Aargsort = Numeric.argsort(Aplot)
                        Mplot = Numeric.asarray(MA.compress(cond, Mn_masked))
                        self.graph.setCurveData(normCurve, Numeric.take(Aplot, Aargsort), Numeric.take(Mplot, Aargsort))
                        pen = QPen(condColors[condIdx],ProbeSet.PenWidthInactiveCurve)
                        self.graph.setCurvePen(normCurve, pen)
                        self.graph.setCurveStyle(normCurve, self.normCurveStyle)

                        # plot a normalization curve consisting only of ticks corresponding to the "right" probes
                        normCurveTicks = self.graph.insertCurve("Norm. curve ticks %i: %s" % (condIdx, str(ncKey)), ncKey)
                        ncData.curveList.append(normCurveTicks)
                        cond_condTicks = MA.logical_and(cond, condTicks)
                        Aplot = Numeric.asarray(MA.compress(cond_condTicks, An_masked))
                        Aargsort = Numeric.argsort(Aplot)
                        Mplot = Numeric.asarray(MA.compress(cond_condTicks, Mn_masked))
                        self.graph.setCurveData(normCurveTicks, Numeric.take(Aplot, Aargsort), Numeric.take(Mplot, Aargsort))
                        pen = QPen(condColors[condIdx],ProbeSet.PenWidthInactiveCurve)
                        self.graph.setCurvePen(normCurveTicks, pen)
                        self.graph.setCurveStyle(normCurveTicks, QwtCurve.NoCurve)
                        # add markers: 5x5 circles
                        qSymbol = QwtSymbol(1, QBrush(QColor(255,255,255), QBrush.SolidPattern), pen, QSize(5,5))
                        self.graph.setCurveSymbol(normCurveTicks, qSymbol)

                # compute normalized log2 ratio for indices ncData.normInd
                Mdata, Adata = self._getMA_masked_indexed(ncData.normInd, False)
                if self.logAxisY:
                    Mdata -= Mn_masked
                else:
                    Mdata /= Mn_masked
                    Mdata = MA.log(Mdata) / math.log(2)
                # store self._Mnorm
                self._Mnorm = MA.where(NumExtn.indices2condition(ncData.normInd, self._Mnorm.shape[0]), Mdata, self._Mnorm)
        if refresh: self.graph.replot()


    def _clearNormCurves(self, refresh=True):
        if D1 or D2 or D6: print "Probes._clearNormCurves"
        changed = False
        for ncData in self._ncdd.values():
            for curve in ncData.curveList:
                # remove curve from self._active list
                plotCurve = self.graph.curve(curve)
                if self._active.has_key(plotCurve.key):
                    self._active.pop(plotCurve.key)
                # remove curve from graph
                self.graph.removeCurve(curve)
                changed = True
            ncData.curveList = []
        self._ncdd = NormCurveDataDict()
        if refresh and changed:
            self.graph.replot()


    ############################################
    # GRAPH: ACTIVATE PROBE & NORM. CURVES
    ############################################

    def switchCurveActive(self, curveKey, refresh=True):
        if D1 or D3 or D6: print "Probes.switchCurveActive"
        self.setCurveActive(curveKey, not(self._active.has_key(curveKey)), refresh)


    def setCurveActiveList(self, curveKeyList, refresh=True):
        """Deactivetes currently active, activates those from the given list;
        curveKeyList : [curveKey1, curveKey2,...] | None
        """
        if D1 or D3 or D6: print "Probes.setCurveActiveList", curveKeyList
        if curveKeyList:
            for curveKey in self._active.keys():
                if curveKey not in curveKeyList:
                    self.setCurveActive(curveKey, False, False)
            for curveKey in curveKeyList:
                self.setCurveActive(curveKey, True, False)                    
        else:
            for curveKey in self._active.keys():
                self.setCurveActive(curveKey, False, False)
        if refresh: self.graph.replot()


    def setCurveActive(self, curveKey, active, refresh=True):
        """activate either probeSet or normalization curve
        """
        if D1 or D3 or D6: print "Probes.setCurveActive"
        # if curveKey represents a normalization curve
        if self._ncdd.has_key(curveKey):
            # activate probe curves that match with the normalization curve
            for probe in self._ncdd[curveKey].probeList:
                self._setProbeCurveActive(probe, active, False)
            # activate the normalization curve
            self._setNormCurveActive(curveKey, active, refresh)
        # if curveKey represents a ProbeSet
        elif self.has_key(curveKey):
            # activate corresponding norm. curve (if exists)
            if self._ncdd.has_key(curveKey):
                self._setNormCurveActive(curveKey, active, False)
            elif self._ncdd.has_key(Probes.NormCurveNameGlobal):
                self._setNormCurveActive(Probes.NormCurveNameGlobal, active, False)
            # activate probe
            self._setProbeCurveActive(self[curveKey], active, refresh)
        else:
            # do not raise error, this can occur on mouseOver when self._ncdd is being updated
            print "Warning: unknown curveKey: %s, known %s" % (str(curveKey), str(self._ncdd.keys()))


    def _setProbeCurveActive(self, probe, active, refresh=True):
        if D1 or D3: print "Probes._setProbeCurveActive"
        if probe.curve <> None and active <> self._active.has_key(probe.pKey):
            if active:
                self._active[probe.pKey] = probe.pKey
            else:
                self._active.pop(probe.pKey)
            self._setProbeCurveSymbol(probe, refresh)


    def _setNormCurveActive(self, curveKey, active, refresh=True):
        if D1: print "Probes._setNormCurveActive"
        if self._ncdd.has_key(curveKey) and self._ncdd[curveKey].curveList and active <> self._active.has_key(curveKey):
            for curve in self._ncdd[curveKey].curveList:
                pen = self.graph.curve(curve).pen() # curve is actually a long curve key
                pen.setWidth(ProbeSet.PenWidths[active])
                self.graph.setCurvePen(curve, pen)
                symbol = self.graph.curveSymbol(curve)
                symbol.setPen(pen)
                self.graph.setCurveSymbol(curve, symbol)
            if active:
                self._active[curveKey] = curveKey
            else:
                self._active.pop(curveKey)
            if refresh and len(self._ncdd[curveKey].curveList) > 0:
                self.graph.replot()
                

    ############################################
    # FILTER
    ############################################

    def getFilter(self):
        if D1: print "Probes.getFilter"
        if type(self.__filterMaxCV) == types.NoneType:
            self._setFilterMaxCV()
        if type(self.__filterMinRatio) == types.NoneType:
            self._setFilterMinRatio()
        if type(self.__filterMaxFGInt) == types.NoneType:
            self._setFilterMaxFGInt()
        if type(self.__filterMaxBGInt) == types.NoneType:
            self._setFilterMaxBGInt()
        return Numeric.logical_or(Numeric.logical_or(Numeric.logical_or(self.__filterMaxCV, self.__filterMinRatio), self.__filterMaxFGInt), self.__filterMaxBGInt)


    def _setFilterMaxCV(self):
        if D1 or D2 or D4: print "Probes._setFilterMaxCV"
        if self.__sigSmpl:
            # maxCV: bgSD / sig <= self.maxCV
            if self.__bgSmplSD <> None and self.__bgRefSD <> None:
                self.__filterMaxCV = MA.asarray(self.__bgSmplSD / self.__sigSmpl).filled(Probes.midVal) > self.maxCV
                self.__filterMaxCV += MA.asarray(self.__bgRefSD / self.__sigRef).filled(Probes.midVal) > self.maxCV
                # convert to 0/1
                self.__filterMaxCV = self.__filterMaxCV > 0
            else:
                self.__filterMaxCV = Numeric.zeros(self.__sigSmpl.shape, Numeric.Int)

    def _setFilterMinRatio(self):
        if D1 or D2 or D4: print "Probes._setFilterMinRatio"
        if self.__sigSmpl:
            # minIntRatio: sig / bg >= self.max
            self.__filterMinRatio = MA.asarray(self.__sigSmpl < self.minIntensityRatio * self.__bgSmpl).filled(1)
            self.__filterMinRatio += MA.asarray(self.__sigRef < self.minIntensityRatio * self.__bgRef).filled(1)
            # convert to 0/1
            self.__filterMinRatio = self.__filterMinRatio > 0

    def _setFilterMaxFGInt(self):
        if D1 or D2 or D4: print "Probes._setFilterMaxFGInt"
        if self.__sigSmpl:
            # maxFGIntensity: sig <= maxFGIntensity
            self.__filterMaxFGInt = MA.asarray(self.__sigSmpl > self.maxFGIntensity).filled(1)
            self.__filterMaxFGInt += MA.asarray(self.__sigRef > self.maxFGIntensity).filled(1)
            # convert to 0/1
            self.__filterMaxFGInt = self.__filterMaxFGInt > 0

    def _setFilterMaxBGInt(self):
        if D1 or D2 or D4: print "Probes._setFilterMaxBGInt"
        if self.__bgSmpl:
            # maxBGIntensity: bg <= maxBGIntensity
            self.__filterMaxBGInt = MA.asarray(self.__bgSmpl > self.maxBGIntensity).filled(1)
            self.__filterMaxBGInt += MA.asarray(self.__bgRef > self.maxBGIntensity).filled(1)
            # convert to 0/1
            self.__filterMaxBGInt = self.__filterMaxBGInt > 0


    ##############################################
    # DATA (accounts for filters)
    # _get..._masked(): returns MA array
    # _get...(): returns compressed Numeric array
    ##############################################

    def _sigSmpl_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__sigSmpl)

    def _sigRef_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__sigRef)

    def _bgSmpl_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__bgSmpl)

    def _bgRef_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__bgRef)

    def _ratio_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__ratio.filled(1))


    def _netSmpl_masked(self, condition):
        return self.__netSmpl_masked_func[self.subtrBG](condition)
    
    def _netRef_masked(self, condition):
        return self.__netRef_masked_func[self.subtrBG](condition)
    

    def _sigSmpl(self, condition):
        return Numeric.asarray(self._sigSmpl_masked(condition).compressed())

    def _sigRef(self, condition):
        return Numeric.asarray(self._sigRef_masked(condition).compressed())

    def _bgSmpl(self, condition):
        return Numeric.asarray(self._bgSmpl_masked(condition).compressed())

    def _bgRef(self, condition):
        return Numeric.asarray(self._bgRef_masked(condition).compressed())

    def _ratio(self, condition):
        return Numeric.asarray(self._ratio_masked(condition).compressed())


    def sigSmpl(self, pKey):
        cond = NumExtn.indices2condition(dict.__getitem__(self, pKey).getDataIndices(), self.__sigSmpl.shape[0])
        return self._sigSmpl(cond)

    def sigRef(self, pKey):
        cond = NumExtn.indices2condition(dict.__getitem__(self, pKey).getDataIndices(), self.__sigRef.shape[0])
        return self._sigRef(cond)

    def bgSmpl(self, pKey):
        cond = NumExtn.indices2condition(dict.__getitem__(self, pKey).getDataIndices(), self.__bgSmpl.shape[0])
        return self._bgSmpl(cond)

    def bgRef(self, pKey):
        cond = NumExtn.indices2condition(dict.__getitem__(self, pKey).getDataIndices(), self.__bgRef.shape[0])
        return self._bgRef(cond)

    def ratio(self, pKey):
        cond = NumExtn.indices2condition(dict.__getitem__(self, pKey).getDataIndices(), self.__ratio.shape[0])
        return self._ratio(cond)


    def __compM_masked(self, netSmpl, netRef, ratios):
        M = netSmpl / netRef / ratios
        if self.logAxisY:
            M = MA.log(M) / math.log(2)
        return M

##    def _getM_masked(self, condition):
##        netSmpl = self._netSmpl_masked(condition)
##        netRef = self._netRef_masked(condition)
##        ratios = self._ratio_masked(condition)
##        return self.__compM_masked(netSmpl, netRef, ratios)
##
##    def _getM(self, condition):
##        return Numeric.asarray(self._getM_masked(condition).compressed())



    def __compA_masked(self, netSmpl, netRef):
        """TODO: check why condition is not used!!!
        """
        A = MA.sqrt(netSmpl*netRef)
        if self.logAxisX:
            A = MA.log(A) / math.log(2)
        return A

    def _getA_masked(self, condition):
        netSmpl = self._netSmpl_masked(condition)
        netRef = self._netRef_masked(condition)
        return self.__compA_masked(netSmpl, netRef)

##    def _getA(self, condition):
##        return Numeric.asarray(self._getA_masked(condition).compressed())



    def _getMA_masked(self, condition, center):
        """Returns MA arrays: M/ratio, A (masked by filter and condition);
        if center: center M by ratios.
        """
        if D1 or D2: print "Probes._getMA_masked"
        netSmpl = self._netSmpl_masked(condition)
        netRef = self._netRef_masked(condition)
        if center:
            ratios = self._ratio_masked(condition)
        else:
            ratios = Numeric.ones(condition.shape)
        return self.__compM_masked(netSmpl, netRef, ratios), self.__compA_masked(netSmpl, netRef)

    def _getMA_masked_indexed(self, indices, center):
        return self._getMA_masked(NumExtn.indices2condition(indices, self.__sigSmpl.shape[0]), center)

    def _getMA_compressed(self, condition, center):
        """Returns Numeric arrays: M/ratio, A (compressed by filter, mask and condition).
        """
        M,A = self._getMA_masked(condition, center)
        noMask = Numeric.logical_not(Numeric.logical_or(MA.getmaskarray(A), MA.getmaskarray(M)))
        return Numeric.asarray(MA.compress(noMask, M)), Numeric.asarray(MA.compress(noMask, A))

    def getMA(self, pKey, center):
        """Returns Numeric arrays: M/ratio, A (compressed by filter, condition and mask)
        """
        if D1 or D2: print "Probes.getMA"
        cond = NumExtn.indices2condition(dict.__getitem__(self, pKey).getDataIndices(), self.__sigSmpl.shape[0])
        return self._getMA_compressed(cond, center)


    ############################################
    # TOOLTIP
    ############################################

    def showDataTooltip(self, curveKey, x, y):
        if D1: print "Probes.showDataTooltip"
        if D4: print "# self.__plotted and self.getFilter(): ", Numeric.add.reduce(self.__plotted == self.getFilter())
        out = ""
        # curve represents a probeSet, curveKey corresponds to probe.pKey
        if self.has_key(curveKey):
            probe = self.get(curveKey)
            out = "%s (%s)" % (str(probe.valA), str(probe.valAAlias))
            if self.varNameB <> "<none>":
                out += ", " + str(probe.valB)
            out += "\nSmpl (signal - bg = net) / Ref (signal - bg = net)\n"
            for ss,sr,bs,br in zip(self.sigSmpl(curveKey), self.sigRef(curveKey), self.bgSmpl(curveKey), self.bgRef(curveKey)):
                out += "%5.f - %5.f = %5.f  /  %5.f - %5.f = %5.f\n" % (ss,bs,ss-bs, sr,br,sr-br)
        # curve represents one of the normalization curves, curveKey corresponds to curve.key
        elif self._ncdd.has_key(curveKey):
            indices = self._ncdd[curveKey].normInd
            out = "%s: %s\nProbes in total:\t%d control,\t%d other.\nAccepted:\t%d control,\t%d other.\nNormalized:\t%d control,\t%d other." % \
                  (self.varNameB, curveKey,
                   self.getNumProbesControls_indexed(self._ncdd[curveKey].computInd), self.getNumProbesOthers_indexed(self._ncdd[curveKey].computInd),
                   self.getNumProbesControls_nonFiltered_indexed(self._ncdd[curveKey].computInd), self.getNumProbesOthers_nonFiltered_indexed(self._ncdd[curveKey].computInd),
                   self.getNumProbesControls_nonFiltered_indexed(self._ncdd[curveKey].normInd), self.getNumProbesOthers_nonFiltered_indexed(self._ncdd[curveKey].normInd))
        else:
            raise ValueError, "Unknown curveKey: %s" % str(curveKey)
        if out:
            out = out[:-1]
            xPoints = self.graph.transform(QwtPlot.xBottom, x)
            yPoints = self.graph.transform(QwtPlot.yLeft, y)
            rect = QRect(xPoints+self.graph.canvas().frameGeometry().x()-self.markerSize/2, yPoints+self.graph.canvas().frameGeometry().y()-self.markerSize/2, self.markerSize, self.markerSize)
            MyQToolTip.setRect(self.graph.tooltip, rect, out)


    ############################################
    # NORMALIZED DATA for plotting curves
    ############################################

    def _getNormData_ActrlMctrl_AM_masked(self, nameInd):
        """returns A and M of controls (Ac, Mc), A of all named probes and normalized M in points A;
        Ac, Mc: compressed; values of controls
        A, M masked; normalization curve computed in all "named" probes (or numPoints from min(Aprobes), max(Aprobes)
        """
        if  D6: print "Probes._getNormData_ActrlMctrl_AM_masked"
        condName = NumExtn.indices2condition(nameInd, self.__sigSmpl.shape[0])
        condControlName = Numeric.logical_and(self._isControlArr(), condName)
        Mc,Ac = self._getMA_compressed(condControlName, True)
        # if A not given, then it should equal to all probes which correspond to name
        A = self._getA_masked(condName)
        M = MA.zeros(A.shape, MA.Float) * MA.masked
        # proceed if we have at least one control probe (we account for _minNumControlProbes in _getNormCurveName2Ind())
        if Mc.shape[0] >= 1:
            Mcompressed = self._approxFunction(A.compressed(), Mc, Ac)
            if Mcompressed:
                MA.put(M, NumExtn.condition2indices(Numeric.logical_not(MA.getmaskarray(A))), Mcompressed)
        return  Ac, Mc, A, M


    def _getNormCurveLoess(self, A, Mc, Ac):
        return NumExtn.lowess2(Ac, Mc, A, f=self.loessWindow/100.)


    def _getNormCurveLinReg(self, A, Mc, Ac):
        X = Numeric.reshape(Numeric.asarray(Ac), (len(Ac),1))
        y = Numeric.asarray(Mc)
        X = Numeric.concatenate((Numeric.ones((X.shape[0],1), Numeric.Float), X), 1)
        XT = Numeric.transpose(X)
        try:
            XTXinv = LinearAlgebra.inverse(Numeric.dot(XT,X))
        except LinearAlgebra.LinAlgError:
            print "Warning: singular matrix, using generalized_inverse"
            XTX = Numeric.dot(XT,X)   # store the singuar matrix
            XTXinv = LinearAlgebra.generalized_inverse(XTX)
        b = Numeric.dot(Numeric.dot(XTXinv, XT), y)
        return b[0] + b[1]*A


    def _getNormCurveMedian(self, A, Mc, Ac):
        if D4 or D5: print "Probes._getNormCurveMedian, value:", NumExtn.median(Mc)
        return Numeric.resize(NumExtn.median(Mc), A.shape)
        

    #################################################
    # MERGE REPLICAS, CONCATENATE LISTS OF STRINGS
    #################################################

    def __mergeReplicasNone(self, ma, mergeFunction, callback):
        callback(100)
        return ma

    def __concatReplicasNone(self, lst, removeDupl):
        return lst


    def __mergeReplicasPerVarsAB(self, ma, mergeFunction, callback):
        """merge by pKey (varA and varB)
        """
        if D1: print "Probes.__mergeReplicasPerVarsAB"
        shp = list(ma.shape)
        shp[0] = len(self.values())
        maMerged = MA.zeros(shp, ma.typecode())
        pbStep = 100./len(self)
        for idx, probe in enumerate(self.values()):
            maMerged[idx] = mergeFunction(MA.take(ma, probe.getDataIndices()))
            callback(pbStep)
        return maMerged

    def __concatReplicasPerVarsAB(self, lst, removeDupl):
        """concatenate list of strings by pKey (var A and var B)
        """
        if D1: print "Probes.__concatReplicasPerVarsAB"
        lstMerged = [""]*len(self.values())
        for idx, probe in enumerate(self.values()):
            # remove duplicates, sort, merge as CSV
            subLst = Numeric.take(Numeric.asarray(lst, Numeric.PyObject), probe.getDataIndices()).tolist()
            if removeDupl:
                subLst = dict(zip(subLst,subLst)).keys()
            subLst.sort()
            lstMerged[idx] = reduce(lambda a,b: "%s, %s" % (a,b), subLst, "")[2:]
        return lstMerged


    def __mergeReplicasPerVarA(self, ma, mergeFunction, callback):
        """merge by var A
        """
        if D1: print "Probes.__mergeReplicasPerVarA"
        shp = list(ma.shape)
        shp[0] = len(self._valA2ind)
        maMerged = MA.zeros(shp, ma.typecode())
        pbStep = 100./len(self._valA2ind)
        for idx, dataInd in enumerate(self._valA2ind.values()):
            maMerged[idx] = mergeFunction(MA.take(ma, dataInd))
            callback(pbStep)
        return maMerged

    def __concatReplicasPerVarA(self, lst, removeDupl):
        """concatenate list of strings by var A
        """
        if D1: print "Probes.__concatReplicasPerVarA"
        lstMerged = [""]*len(self._valA2ind)
        for idx, dataInd in enumerate(self._valA2ind.values()):
            # remove duplicates, sort, merge as CSV
            subLst = Numeric.take(Numeric.asarray(lst, Numeric.PyObject), dataInd).tolist()
            if removeDupl:
                subLst = dict(zip(subLst,subLst)).keys()
            subLst.sort()
            lstMerged[idx] = reduce(lambda a,b: "%s, %s" % (a,b), subLst, "")[2:]
        return lstMerged


    ############################################
    # NORMALIZED DATA for output
    ############################################

##    def getLog2Ratio_normCentered_masked(self, mergeLevel, mergeType, callback):
##        """Returns masked array of normalized and centered log2ratio of individual probes;
##        accounts for filters and for ratios;
##        mergeType: 0:None, 1:mean, 2:median
##        """
##        if D6: print "Probes.getLog2Ratio_normCentered_masked"
##        return self._mergeFunc[mergeLevel](self._MnormCentered, Probes.mergeTypes[mergeType], callback)
    

    def getLog2Ratio_norm_masked(self, mergeLevel, mergeType, callback, center):
        """Returns masked array of normalized log2ratio of individual probes;
        accounts for filters, but NOT for ratios;
        mergeType: 0:None, 1:mean, 2:median
        """
        if D6: print "Probes.getLog2Ratio_norm_masked"
        if center:
            Mnorm = self._Mnorm - MA.log(self.__ratio.filled(1))/math.log(2)
        else:
            Mnorm = self._Mnorm
        return self._mergeFunc[mergeLevel](Mnorm, Probes.mergeTypes[mergeType], callback)


    def getNetIntensity_smpl_ref(self, mergeLevel, mergeType, callback):
        """For output, return MA array (#probes, 2) where columns correspond to sample & reference signals.
        """
        if D6: print "Probes.getNetIntensity_smpl_ref"
        netS = self._netSmpl_masked(1)
        netR = self._netRef_masked(1)
        iSR = MA.concatenate([MA.reshape(netS, (netS.shape[0], 1)), MA.reshape(netR, (netR.shape[0], 1))], 1)
        # merge and return
        return self._mergeFunc[mergeLevel](iSR, Probes.mergeTypes[mergeType], callback)
    

    def getLog2Ratio_raw_masked(self, mergeLevel, mergeType, callback, center):
        """returns non-normalized log2 ratio, accounts for filters
        """
        if D6: print "Probes.getLog2Ratio_raw_masked"
        if center:
            l2r = MA.log(self._netSmpl_masked(1) / self._netRef_masked(1) / self._ratio_masked(1)) / math.log(2)
        else:
            l2r = MA.log(self._netSmpl_masked(1) / self._netRef_masked(1)) / math.log(2)
        # merge and return
        return self._mergeFunc[mergeLevel](l2r, Probes.mergeTypes[mergeType], callback)


##    def getLog2Ratio_rawCentered_masked(self, mergeLevel, mergeType, callback):
##        """returns non-normalized centered log2 ratio, accounts for filters
##        """
##        if D6: print "Probes.getLog2Ratio_rawCentered_masked"
##        # merge and return
##        return self._mergeFunc[mergeLevel](l2r, Probes.mergeTypes[mergeType], callback)
##

    def getA_masked(self, mergeLevel, mergeType, callback):
        """returns log2 average (net) intensity, accounts for filters
        """
        if D6: print "Probes.getA_masked"
        A = MA.log(MA.sqrt(self._netSmpl_masked(1)*self._netRef_masked(1))) / math.log(2)
        # merge and return
        return self._mergeFunc[mergeLevel](A, Probes.mergeTypes[mergeType], callback)


    def getControlRatios(self, mergeLevel, mergeType, callback):
        """returns ratios of controls, DK for others; does not account for filters
        """
        return self._mergeFunc[mergeLevel](self.__ratio, Probes.mergeTypes[mergeType], callback)
        


    ############################################
    # valsA, valsB
    ############################################

    def getValsA(self, mergeLevel):
        if D1: print "Probes.getValsA"
        if mergeLevel == OWNormalize.MergeLevelNone:
            return self._valAList
        elif mergeLevel == OWNormalize.MergeLevelPerVarsAB:
            return map(lambda x: x.valA, self.values())
        elif mergeLevel == OWNormalize.MergeLevelPerVarA:
            return self._valA2ind.keys()
        else:
            raise AttributeError, "unknown merge level: %s" % str(mergeLevel)


    def getValsAAlias(self, mergeLevel):
        if D1: print "Probes.getValsAAlias"
        aliasList = list(self._valAList)
        for pr in self.values():
            for idx in pr.getDataIndices():
                aliasList[idx] = pr.valAAlias
        return self._concatFunc[mergeLevel](aliasList, removeDupl=True)


    def getValsB(self, mergeLevel):
        if mergeLevel == OWNormalize.MergeLevelNone:
            return self._valBList
        elif mergeLevel == OWNormalize.MergeLevelPerVarsAB:
            return map(lambda x: x.valB, self.values())
        elif mergeLevel ==OWNormalize.MergeLevelPerVarA:
            raise ValueError, "cannot return probe vals B if mergeLevel == %i" % OWNormalize.MergeLevelPerVarA
        else:
            raise AttributeError, "unknown merge level: %s" % str(mergeLevel)




if __name__=="__main__":

    def test_minNumControlProbes(numC):
        import Numeric, LinearAlgebra, statc
        A = [-1,0,0.5,1,1.5, 3]
        Mc = Numeric.arange(0,numC,1.)
        Ac = Numeric.arange(0,numC,1.)
        X = Numeric.reshape(Numeric.asarray(Ac), (len(Ac),1))
        y = Numeric.asarray(Mc)
        X = Numeric.concatenate((Numeric.ones((X.shape[0],1), Numeric.Float), X), 1)
        XT = Numeric.transpose(X)
        try:
            XTXinv = LinearAlgebra.inverse(Numeric.dot(XT,X))
        except LinearAlgebra.LinAlgError:
            print "Warning: singular matrix, using generalized_inverse"
            XTX = Numeric.dot(XT,X)   # store the singuar matrix
            XTXinv = LinearAlgebra.generalized_inverse(XTX)
        b = Numeric.dot(Numeric.dot(XTXinv, XT), y)
        print "linReg: ", b[0] + b[1]*Numeric.asarray(A)
        print "loess:  ", Numeric.asarray(statc.loess(zip(Ac, Mc), Numeric.asarray(A).tolist(), 0.5))[:,1]


    def test_widget():        
        import OWDataTable, orngSignalManager
        signalManager = orngSignalManager.SignalManager(0)
        a=QApplication(sys.argv)
        ow=OWNormalize(signalManager = signalManager)
        a.setMainWidget(ow)
        ow.show()

        # settings    
##        ow.outNonNormLogRatio = True
##        ow.normRange = Probes.NormRangeLocal
##        ow.approxFunction = 0

        # DATA 1: horizontal line in the middle of the slide
##        ow.onDataInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\PB, cholesterol\Tadeja 2nd image analysis\10vs10mg original data\0449yPos.txt", DC="<NO DATA>"))
####        ow.varNameB = "yPos"
####        ow.varABChange()
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\Sterolgene v0 mouse probeRatios (ID v0).tab"))

##        # DATA 2: extremely low signal (only few genes pass the filters)
##        ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\PB, cholesterol\Tadeja drago\05vs10mg\chol.diet\2537.txt"))
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\sterolgene v.0 mouse controlGeneRatios 2.tab"))
##
##        # DATA 3: high noise
##        ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\PB, cholesterol\Tadeja drago\05vs10mg\control\6033.txt"))
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\sterolgene v.0 mouse controlGeneRatios 2.tab"))
##
##        # DATA 4: krizstina
##        ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 human\Krisztina\13217291-A01.txt", DC="<NO DATA>", noClass=1, noCodedDiscrete=1))
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 human\Sterolgene v1 ControlGeneRatios REVERSED 2.tab"))
##
##        # DATA 5: patologija, predzadnje spotiranje (redcene kontrole)
##        ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\MF other\Patologija\2005-11-21\136926 results.txt", DC="<NO DATA>"))
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\MF other\Patologija\2005-11-21\gasper controlRatios genes 2.tab"))
##
##        # DATA 6: patologija, zadnje spotiranje
##        ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\MF other\Patologija\2006-01-18\2006-01-18-rezultati.txt", DC="<NO DATA>"))
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\MF other\Patologija\2005-11-21\gasper controlRatios genes 2.tab"))

##        # DATA 7: Agilent
##        ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\TNF, starved\Agilent\6135_A01.txt"))

        # DATA 8: Krisztina, June 2006
##        ow.onDataInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 human\Krisztina\2006-06-09\Analysis\13217311-top.txt"))
##        ow.onDataInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 human\Krisztina\2006-06-09\Analysis\13217311-bottom.txt"))
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 human\Sterolgene v1 ControlGeneRatios 2.tab"))

##        # DATA 9: Jana, 28.6.2006
##        ow.onDataInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 mouse\Jana\2006-07-11 data\PMEA\13221205bottom.txt"))
##        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 mouse\Sterolgene v1 mouse probeRatios (ID).tab"))

        # DATA 10: Viola, 2006-07-27
        ow.defNameA = "SteroltalkID"
        ow.defNameB = ""
        ow.defNameSmpl1 = "635"
        ow.defNameSmpl2 = ""
        ow.defNameRef1 = "532"
        ow.defNameRef2 = ""
        ow.defNameForeground = "F"
        ow.defNameBackground = "B"
        ow.defNameMean = "Median"
        ow.defNameSD = "SD"
        ow.onDataInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 mouse\2006-07-26 Viola Tamasi\2007-03 data from Viola\PB CAR- 01.tab", DK="Error"))
        ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 mouse\ST1m probeRatios Lucidea Luty1,2.tab"))

        # OWDataTable
        dt = OWDataTable.OWDataTable(signalManager = signalManager)
        signalManager.addWidget(ow)
        signalManager.addWidget(dt)
        signalManager.setFreeze(1)
##        signalManager.addLink(ow, dt, "Probe Data", "Examples", 1)
        signalManager.addLink(ow, dt, "Expression Data", "Examples", 1)
        signalManager.setFreeze(0)
        dt.show()

##        # save
##        orange.saveTabDelimited(r"c:\Documents and Settings\peterjuv\My Documents\Orange\OWNormalize\test comp 2\output.tab", dt.data.values()[0])

        # exec, save settings 
        a.exec_loop()
        ow.saveSettings()


##############################################################
##    test_minNumControlProbes(2)
    test_widget()
