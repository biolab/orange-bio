"""
<name>Normalize Microarray Data</name>
<description>Normalization of custom cDNA microarray data.</description>
<icon>icons/Unknown.png</icon>
<priority>1150</priority>
<author>Peter Juvan (peter.juvan@fri.uni-lj.si)</author>
"""

"""
TODO: settingsList
"""

import string, math
import Numeric, MA
import LinearAlgebra
import NumExtn
import orange
from qttable import *
from OWWidget import *
import OWGUI, OWToolbars
from OWGraph import *
from OWGraphTools import *      # color palletes, user defined curves, ...
import ColorPalette             # ColorButton
import statc
import Types

import chipstat

##from Meda.Anova import MultLinReg


# global debugging variables
D1 = False
D2 = False
D3 = False
D4 = False

class OWNormalize(OWWidget):

##    settingsList = [""]
    # constants
    stylesIndices = zip(["Circle","Rect","Diamond","Triangle","DTriangle","UTriangle","LTriangle","RTriangle","Cross","XCross"], range(1,12))
    sizeButtonColor = 18
    tcPKey = 0      # hidden column for probe keys
    tcMarker = 1
    tcRatio = 2
    tcID = 3
    tcName = 4      # hidden if probe name is <none>


    def __init__(self, parent = None, signalManager = None, name = "Normalize Microarray Data"):
        if D1: print "OWNormalize.__init__"
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True, wantStatusBar=False)  #initialize base class

        # set channels
        self.inputs = [("Examples", ExampleTable, self.onDataInput), ("Probes", ExampleTable, self.onProbesInput)]
        self.outputs = [("Examples", ExampleTable), ("Probes", ExampleTable)]

        # general settings
        self.controlName = ""
        self.ratioStr = ""
        self._probeSymbolComboIdx = 0    # index for cmbProbeSymbol, equals to index+1 for QwtSymbol.Style
        self.probeColor = QColor(0,0,255)
        self.subtrBG = False
        self.useCV = True
        self.maxCV = 0.5
        self.useMinIntensity = True
        self.minIntensityRatio = 1.5
        self.useMaxIntensity = True
        self.maxIntensity = 65536
        self.normRange = 0  # 0: global, 1: local, per probe names
        self.normType = 2   # 0: median, 1: LR, 2: LOESS
        self.loessWindow = 60
        self.loessWeight = 0.2
        # graph
        self.logAxisX = True
        self.logAxisY = True
        self.markerSize = 9
        self.mergeReplGraph = False
        self.showLegend = False
        self.tracking = True
        # output
##        self.mergeRepl = True
        self.mergeReplType = 1      #0: none, 1: mean, 2: median
        self.mergeOtherType = 0     #0: use first, 1: concatenate
        self.outNumProbes = True
        self.outNetSignal = True
        self.outNonNormLogRatio = False
        self.autoSendSelection = 1
        # load previous settings
        self.loadSettings()
        
        # context-specific settings
        self.data = None
        self.dataProbes = None
        self.varsAll = {}
        # variable indices selected within combos
        self.varNameID = None
        self.varNameName = "<none>"
        self.varNameSignalSmpl = None
        self.varNameSignalRef = None
        self.varNameBGSmpl = None
        self.varNameBGRef = None
        self.varNameBGSmplSD = "<none>"
        self.varNameBGRefSD = "<none>"
        # names of selected vars from listBox
        self.varsOtherSelected = {}

        # GUI
        self.controlArea.setFixedWidth(260)
        self.resize(1000, 800)

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
        self.tabs.insertTab(boxVars, "Variables")
        self.cmbVarID = OWGUI.comboBox(boxVars, self, "varNameID", "Probe ID", callback=self.varIDChange, sendSelectedValue=1, valueType=str)
        self.cmbVarName = OWGUI.comboBox(boxVars, self, "varNameName", "Probe Name", callback=self.varIDChange, sendSelectedValue=1, valueType=str)
        self.cmbVarSignalSmpl = OWGUI.comboBox(boxVars, self, "varNameSignalSmpl", "Sample foreground intensity", callback=lambda x="varSignalSmpl": self.varDataChange(x), sendSelectedValue=1, valueType=str)
        self.cmbVarSignalRef = OWGUI.comboBox(boxVars, self, "varNameSignalRef", "Reference foreground intensity", callback=lambda x="varSignalRef": self.varDataChange(x), sendSelectedValue=1, valueType=str)
        self.cmbVarBGSmpl = OWGUI.comboBox(boxVars, self, "varNameBGSmpl", "Sample background intensity", callback=lambda x="varBGSmpl": self.varDataChange(x), sendSelectedValue=1, valueType=str)
        self.cmbVarBGRef = OWGUI.comboBox(boxVars, self, "varNameBGRef", "Reference background intensity", callback=lambda x="varBGRef": self.varDataChange(x), sendSelectedValue=1, valueType=str)
        self.cmbVarBGSmplSD = OWGUI.comboBox(boxVars, self, "varNameBGSmplSD", "Sample background std. dev.", callback=lambda x="varBGSmplSD": self.varSDChange(x), sendSelectedValue=1, valueType=str)
        self.cmbVarBGRefSD = OWGUI.comboBox(boxVars, self, "varNameBGRefSD", "Reference background std. dev.", callback=lambda x="varBGRefSD": self.varSDChange(x), sendSelectedValue=1, valueType=str)
        OWGUI.button(boxVars, self, "Search for Default Variables", callback=self.defaultVarAssignmentClick)

        # tab 2: table probe/ratio/marker
        boxProbes = QVGroupBox(boxVars)
        self.tabs.insertTab(boxProbes, "Probes")
        self.tblControls = QTable(boxProbes)
        self.tblControls.setNumCols(5)
        self.tblControls.setColumnWidth(OWNormalize.tcID, 15)
        self.tblControls.setColumnWidth(OWNormalize.tcName, 15)
        self.tblControls.setColumnWidth(OWNormalize.tcMarker, 20)
        self.tblControls.setColumnWidth(OWNormalize.tcRatio, 15)
        self.connect(self.tblControls, SIGNAL("valueChanged(int,int)"), self.tblControlsValueChange)
        self.connect(self.tblControls, SIGNAL("currentChanged(int, int)"), self.tblControlsCurrentChanged)
        self.connect(self.tblControls , SIGNAL('selectionChanged()'), self.tblControlsSelectionChanged)
        self.connect(self.tblControls, SIGNAL("doubleClicked(int, int, int, const QPoint &)"), self.tblControlsDoubleClicked)
        hheader=self.tblControls.horizontalHeader()
        self.connect(hheader,SIGNAL("clicked(int)"), self.tblControlsHHeaderClicked)
        self.sortby = 0
        hheader.setLabel(OWNormalize.tcID, "Probe ID")
        hheader.setLabel(OWNormalize.tcName, "Probe Name")
        hheader.setLabel(OWNormalize.tcMarker, "")
        hheader.setLabel(OWNormalize.tcRatio, "Ratio")
        hheader.setMovingEnabled(False)
        # hide vertical header and columns pKey, name
        self.tblControls.setLeftMargin(0)
        self.tblControls.verticalHeader().hide()
        self.tblControls.hideColumn(OWNormalize.tcPKey)
        self.tblControls.hideColumn(OWNormalize.tcName)

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
        self.cmbProbeSymbol = OWGUI.comboBox(boxBtns11, self, "_probeSymbolComboIdx", callback=self.cmbProbeSymbolActivated)
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

        # tab 3: settings
        boxSettings = QVGroupBox(self)
        self.tabs.insertTab(boxSettings, "Settings")

        # tab 3: settings: filters
        boxFilters = QVGroupBox('Filters', boxSettings)
        self.cbSubtrBG = OWGUI.checkBox(boxFilters, self, "subtrBG", "Subtract background", callback=self.settingsSubstrBGChange)
        # tab 3: settings: filters: CV
        self.boxMaxCV = QVGroupBox('Max. coeff. of variation (CV)', boxFilters)
        OWGUI.checkBox(self.boxMaxCV, self, "useCV", "Enabled", callback=self.settingsFilterMaxCVChange)
        sldMaxCV = OWGUI.qwtHSlider(self.boxMaxCV, self, "maxCV", minValue=0, maxValue=2, step=0.01, precision=2, callback=None, logarithmic=0, ticks=0, maxWidth=140)
        self.connect(sldMaxCV, SIGNAL("sliderReleased()"), self.settingsFilterMaxCVChange)
        self.lblInfoFilterMaxCV = QLabel("\n", self.boxMaxCV)
        # tab 3: settings: filters: minIntensityRatio
        boxMinIntRatio = QVGroupBox('Min. signal to background ratio', boxFilters)
        OWGUI.checkBox(boxMinIntRatio, self, "useMinIntensity", "Enabled", callback=self.settingsFilterMinIntRatioChange)
        sldMinInt = OWGUI.qwtHSlider(boxMinIntRatio, self, "minIntensityRatio", minValue=0, maxValue=5, step=0.01, precision=2, callback=None, logarithmic=0, ticks=0, maxWidth=140)
        self.connect(sldMinInt, SIGNAL("sliderReleased()"), self.settingsFilterMinIntRatioChange)
        self.lblInfoFilterMinIntRatio = QLabel("\n", boxMinIntRatio)
        # tab 3: settings: filters: maxIntensity
        boxMaxIntensity = QVGroupBox('Max. foreground intensity', boxFilters)
        OWGUI.checkBox(boxMaxIntensity, self, "useMaxIntensity", "Enabled", callback=self.settingsFilterMaxIntChange)
        sldMaxInt = OWGUI.qwtHSlider(boxMaxIntensity, self, "maxIntensity", minValue=0, maxValue=65536, step=1, precision=0, callback=None, logarithmic=0, ticks=0, maxWidth=140)
        self.connect(sldMaxInt, SIGNAL("sliderReleased()"), self.settingsFilterMaxIntChange)
        self.lblInfoFilterMaxInt = QLabel("\n", boxMaxIntensity)
        # tab 3: settings: filters: info about all filters
##        self.lblInfoFiltersAll = QLabel("\n\n", boxFilters)
        

        # tab 3: settings: normalization
        boxNorm = QVGroupBox('Normalization', boxSettings)
        # tab 3: settings: normalization: range, type
        self.boxNormRange = OWGUI.radioButtonsInBox(boxNorm, self, value="normRange", box='Range', btnLabels=["Global", "Local, per probe name"], callback=self.settingsNormalizationChange)
        self.boxNormRange.setEnabled(False)
        boxNormType = OWGUI.radioButtonsInBox(boxNorm, self, value="normType", box='Approx. function', btnLabels=["Median (intensity independent)", "Linear regression", "Loess"], callback=self.settingsNormalizationChange)
        # tab 3: settings: normalization type: loess settings
        sldLoessWindow = OWGUI.qwtHSlider(boxNormType, self, "loessWindow", box="Window size (% of points)", minValue=1, maxValue=100, step=1, precision=0, logarithmic=0, ticks=0, maxWidth=140)
        self.connect(sldLoessWindow, SIGNAL("sliderReleased()"), self.settingsNormalizationChange)
        boxSldLoessWeight = QHBox(boxNormType)
        sldLoessWeight = OWGUI.qwtHSlider(boxSldLoessWeight, self, "loessWeight", box="Weight of non-control probes [0,1]", minValue=0, maxValue=1, step=0.01, precision=2, logarithmic=0, ticks=0, maxWidth=140)
        self.connect(sldLoessWeight, SIGNAL("sliderReleased()"), self.settingsNormalizationChange)
        boxSldLoessWeight.setEnabled(False)

        # tab 4: output
        boxOutput = QVGroupBox(self)
        self.tabs.insertTab(boxOutput, "Output")
        # tab 3: settings: graph
        boxGraph = QVGroupBox('Graph', boxOutput)
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
        self.zoomSelectToolbar = OWToolbars.ZoomSelectToolbar(self, boxGraph, self.graph, self.autoSendSelection)
        # tab 4: output: merge replicas
        self.rbgMergeReplType = OWGUI.radioButtonsInBox(boxOutput, self, value="mergeReplType", btnLabels=["None", "Mean", "Median"], box="Merge probe intensities", callback=self.settingsReplicasChange)
        # tab 4: output: other variables
        boxOtherVars = QVGroupBox('Other variables', boxOutput)
        self.lbVarOthers = QListBox(boxOtherVars)
        self.lbVarOthers.setSelectionMode(QListBox.Multi)
        self.connect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
        self.rbgMergeOtherType = OWGUI.radioButtonsInBox(boxOtherVars, self, value="mergeOtherType", btnLabels=["Use first value", "Concatenate values"], box="Merge other variables", callback=self.settingsReplicasChange)
        # tab 4: output: other additional variables
        self.cbOutNumProbes = OWGUI.checkBox(boxOutput, self, "outNumProbes", "Number of probes", callback=self.settingsOutputChange)
        OWGUI.checkBox(boxOutput, self, "outNetSignal", "Net intensities", callback=self.settingsOutputChange)
        OWGUI.checkBox(boxOutput, self, "outNonNormLogRatio", "Non-normalized log2 ratio", callback=self.settingsOutputChange)

        # control area: info
        boxProbeInfo = QVGroupBox("Info", self.controlArea)
        self.lblProbeInfo = QLabel("\n\n", boxProbeInfo)

        # INITIALIZATION: controls/ratios, probe info, filters, filter info
        self.probes = Probes(self.graph, self.subtrBG, self.logAxisX, self.logAxisY, self.markerSize)
        self.setInfoProbes()
        self.setProbeFilters()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
        self.setNormalization()


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
        self.varsAll = {}
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
                varName = var.name
                if self.varsAll.has_key(varName):
                    print "Warning: domain contains two variables with the same name: %s; the first one will be used"  % varName
                else:
                    self.varsAll[varName] = var
        # fill combos and listBox with variables
        self.fillCmbVars()
        self.fillLbVarOthers()
        self.searchDefaultVarAssignment()
        self.initProbes()
        self.setNormalization()
        # send data
        self.sendData()
        qApp.restoreOverrideCursor()


    def defaultVarAssignmentClick(self):
        """Select default variables in combos based on their names.
        """
        if D1: print "OWNormalize.defaultVarAssignmentClick"
        if self.data:
            qApp.restoreOverrideCursor()
            qApp.setOverrideCursor(QWidget.waitCursor)
            self.searchDefaultVarAssignment()
            self.initProbes()
            self.setNormalization()
            # send data
            self.sendData()
            qApp.restoreOverrideCursor()
            
        
    def fillCmbVars(self):
        """ Fills combos with variables; add "<none>":None for 2nd Probe combo.
        """
        if D1: print "OWNormalize.fillCmbVars"
        cmbNamesAllVars = ["cmbVarID", "cmbVarName"]
        cmbNamesFloatVars = ["cmbVarSignalSmpl", "cmbVarSignalRef", "cmbVarBGSmpl", "cmbVarBGRef", "cmbVarBGSmplSD", "cmbVarBGRefSD"]
        if self.data:
            varNamesAllVars = zip(map(lambda x: x.lower(), self.varsAll.keys()), self.varsAll.keys())
            varNamesAllVars.sort()
            varNamesFloatVars = []
            for varName, var in self.varsAll.items():
                if var.varType == orange.VarTypes.Continuous:
                    varNamesFloatVars.append((varName.lower(),varName))
            varNamesFloatVars.sort()
            # fill Probe ID & Probe Name combos
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
        self.cmbVarName.insertItem("<none>", 0)
        self.cmbVarBGSmplSD.insertItem("<none>", 0)
        self.cmbVarBGRefSD.insertItem("<none>", 0)


    def searchDefaultVarAssignment(self):
        """Select default variables in combos based on their names.
        """
        if D1 or D2: print "OWNormalize.searchDefaultVarAssignment"
        if self.data and len(self.varsAll) > 0:
            # separate continuous and float variables
            varsFloat = {}
            varsEnum = {}
            for varName, var in self.varsAll.items():
                if var.varType == orange.VarTypes.Continuous:
                    varsFloat[varName] = var
                elif var.varType == orange.VarTypes.Discrete:
                    varsEnum[varName] = var
            # smart variable assignment: ID
            self.varNameID = None
            for name in varsEnum:
                if "id" in name.lower():
                    self.varNameID = name
                    break
            if self.varNameID == None:
                if len(self.varsEnum) > 0:
                    self.varNameID = self.varsEnum.keys()[0]
                else:
                    self.varNameID = self.varsAll.keys()[0]
            # smart variable assignment: Name
            self.varNameName = "<none>"
            for name in varsEnum:
                if "name" in name.lower():
                    self.varNameName = name
                    break
            # smart variable assignment: signal, background, background s.d. (smpl & ref)
            colNames = ["raw intensity (med) {smpl}", "raw intensity (med) {ref}", "background (med) {smpl}", "background (med) {ref}", "background (st.dev.) {smpl}", "background (st.dev.) {ref}"]
            vars = [None, None, None, None, None, None]
            varsFloatNames = varsFloat.keys()
            for cIdx, cName in enumerate(colNames):
                for name in varsFloatNames:
                    if cName in name.lower():
                        vars[cIdx] = name
                        varsFloatNames.remove(name)
                        break
                if vars[cIdx] == None:
                    if len(varsFloatNames) > 0:
                        vars[cIdx] = varsFloatNames[0]
                        varsFloatNames.pop(0)
                    else:
                        vars[cIdx] = self.varsAll.keys()[0]
            # select vars in combos, update listbox with other variables
            self.varNameSignalSmpl = vars[0]
            self.varNameSignalRef = vars[1]
            self.varNameBGSmpl = vars[2]
            self.varNameBGRef = vars[3]
            self.varNameBGSmplSD = vars[4]
            self.varNameBGRefSD = vars[5]
            # enable/disable subtract BG checkbox & Max. CV slider
##            self.cbSubtrBG.setEnabled(self.varNameBGSmplSD != "<none>" and self.varNameBGRefSD != "<none>")
            self.boxMaxCV.setEnabled(self.varNameBGSmplSD != "<none>" and self.varNameBGRefSD != "<none>")
            # select other variables where name contains "name"
            self.disconnect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
            for idx in range(self.lbVarOthers.count()):
                if "name" in self.lbVarOthers.item(idx).text().lower():
                    self.lbVarOthers.setSelected(idx, True)
            self.fillVarsOtherSelected()
            self.connect(self.lbVarOthers , SIGNAL('selectionChanged()'), self.varOthersChange)
            


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
        if D1 or D2: print "OWNormalize.initProbes"
        if self.data:
            self.probes.initProbes(self.data, self.varNameID, self.varNameName, self.varNameSignalSmpl, self.varNameSignalRef, self.varNameBGSmpl, self.varNameBGRef, self.varNameBGSmplSD, self.varNameBGRefSD)
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
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()


    def sendData(self):
        """Compute norm. factors, plot them, normalize data and send out normalized data.
        """
        if D1 or D2: print "OWNormalize.sendData"
        if self.probes:
            # normalized log2 ratios -> new example table
            varList = [orange.FloatVariable("Log2 ratio")]            
            l2r = self.probes.getNormalizedLog2Ratio_masked(self.mergeReplType)
            maData = MA.reshape(l2r, (l2r.shape[0], 1))
            if self.outNumProbes:
                varList += [orange.FloatVariable("Num. probes"), orange.FloatVariable("Num. accepted probes")]
                maData = MA.concatenate([maData, self.probes.getNumReplicas_nonFiltered(self.mergeReplType)], 1)
            if self.outNetSignal:
                varList += [orange.FloatVariable("Net intensity (Smpl)"), orange.FloatVariable("Net intensity (Ref)")]
                maData = MA.concatenate([maData, self.probes.getNetIntensity_smpl_ref(self.mergeReplType)], 1)
            if self.outNonNormLogRatio:
                varList.append(orange.FloatVariable("Non-normalized log2 ratio"))
                maData = MA.concatenate([maData, MA.reshape(self.probes.getRawLog2Ratio_masked(self.mergeReplType), (maData.shape[0], 1))], 1)
            etNew = chipstat.ma2orng(maData, orange.Domain(varList, None))

            # data table with pKey, ID and name
            varList = [orange.StringVariable("pKey"), self.data.domain[self.varNameID]]
            IDList = self.probes.getIDs(self.mergeReplType)
            if self.varNameName <> "<none>":
                varList.append(self.data.domain[self.varNameName])
                nameList = self.probes.getNames(self.mergeReplType)
                valListList = map(lambda i,n: [str(i)+str(n), i, n], IDList, nameList)
            else:
                valListList = map(lambda i: [str(i), i], IDList)
            etPKeyIDName = orange.ExampleTable(orange.Domain(varList, None), valListList)            

            # etSubset: subset of examples/attributes from self.data
            domainOtherVars = orange.Domain(self.varsOtherSelected.values(), None)
            if self.mergeReplType:
                if self.mergeOtherType == 0:
                    # use first value: subset of examples
                    etSubset = orange.ExampleTable(self.data.domain)
                    for e in etPKeyIDName:
                        dataIdx0 = self.probes[str(e["pKey"])].getDataIndices()[0]
                        etSubset.append(self.data[dataIdx0])
                    # subset of attributes
                    etOtherVars = orange.ExampleTable(domainOtherVars, etSubset)
                else:
                    # create new string attributes, concatenate values
                    strDomain = orange.Domain(map(lambda varName: orange.StringVariable(varName+" List"), self.varsOtherSelected.keys()), None)
                    etOtherVars = orange.ExampleTable(strDomain)
                    for probe in self.probes.values():
                        vals = [[]]*len(self.varsOtherSelected)
                        for eIdx in probe.getDataIndices():
                            for vIdx, vName in enumerate(self.varsOtherSelected.keys()):
                                vals[vIdx].append(str(self.data[eIdx][vName].native()))
                        vals = map(lambda l: string.join(l, ", "), vals)
                        etOtherVars.append(orange.Example(strDomain, vals))
            else:
                etOtherVars = orange.ExampleTable(domainOtherVars, self.data)

            # final example table (leave out pKey)
            etOut = orange.ExampleTable([orange.ExampleTable(orange.Domain(varList[1:], None), etPKeyIDName), etNew, etOtherVars])
            etOut.name = (self.data.name + " normalized").strip()
        else:
            etOut = None
        self.send("Examples", etOut)
                

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
        # process external probe data
        self.processDataProbes()
        # fill / update probe table & info
        self.fillProbeTable()
        self.setInfoProbes()
        # filters info
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
        # send data and probes data
        self.sendData()
        self.sendProbes()
        qApp.restoreOverrideCursor()


##    def processDataProbes(self):
##        """Copy data from orange.ExampleTable to self.probes.
##        """
##        if D1 or D2: print "OWNormalize.processDataProbes"
##        if self.dataProbes and self.varNameID:
##            varNames = map(lambda var: var.name, self.dataProbes.domain.variables)
##            if self.varNameID in varNames:
##                for e in self.dataProbes:
##                    if self.varNameName <> "<none>" and self.varNameName in varNames:
##                        pKey = str(e[self.varNameID].native()) + str(e[self.varNameName].native())
##                    else:
##                        pKey = str(e[self.varNameID].native())
##                    c = str(e["ColorRGB"])
##                    color = QColor(int(c[0:2],16), int(c[2:4],16), int(c[4:6],16), QColor.Rgb)
##                    try:
##                        symbol = int(e["Symbol"].native())
##                    except TypeError:
##                        symbol = ProbeSet.NoSymbol
##                    ratio = str(e["Ratio"])
##                    if self.probes.has_key(pKey):
##                        print "self.probes.has_key(%s)" %pKey
##                        self.probes.setMarker(pKey, symbol, color, refresh=False)
##                        self.probes.setRatio(pKey, ratio, replot=False)
##                    else:
##                        print "getKeysStartsWith: %s" % pKey
##                        for pk in self.probes.getKeysStartsWith(pKey):
##                            self.probes.setMarker(pk, symbol, color, refresh=False)
##                            self.probes.setRatio(pk, ratio, replot=False)
##                if len(self.dataProbes) > 0:
##                    self.probes.replotAllCurves(refresh=False)
##                    self.probes.replotNormCurve(refresh=True)
##            else:
##                print "Warning: probe data is missing attribute named " + str(self.varNameID)

    def processDataProbes(self):
        """Copy data from orange.ExampleTable to self.probes.
        """
        if D1 or D2: print "OWNormalize.processDataProbes"
        if self.dataProbes:
            varNames = map(lambda var: var.name, self.dataProbes.domain.variables)
            if "Probe ID" in varNames and "ColorRGB" in varNames and "Ratio" in varNames and "Symbol" in varNames:
                # list of probe keys
                if self.varNameName <> "<none>" and "Probe Name" in varNames:
                    pKeysIDsNames = map(lambda e: (str(e["Probe ID"].native()) + str(e["Probe Name"].native()), str(e["Probe ID"].native()), str(e["Probe Name"].native())), self.dataProbes)
                else:
                    pKeysIDsNames = map(lambda e: (str(e["Probe ID"].native()), str(e["Probe ID"].native()), ""), self.dataProbes)
                # get color, symbol and ratio of probes
                for e, (pKey, pID, pName) in zip(self.dataProbes, pKeysIDsNames):
                    c = str(e["ColorRGB"])
                    color = QColor(int(c[0:2],16), int(c[2:4],16), int(c[4:6],16), QColor.Rgb)
                    try:
                        symbol = int(e["Symbol"].native())
                    except TypeError:
                        symbol = ProbeSet.NoSymbol
                    ratio = str(e["Ratio"])
                    # if pKey exists
                    if self.probes.has_key(pKey):
##                        print "self.probes.has_key(%s)" %pKey
                        self.probes.setMarker(pKey, symbol, color, refresh=False)
                        self.probes.setRatio(pKey, ratio, replot=False)
                    # if pKey does not exist
                    else:
##                        print "getKeysFromIDsNames(%s, %s): %s" % (pID, pName, self.probes.getKeysFromIDsNames(pID, pName))
                        for pk in self.probes.getKeysFromIDsNames(pID, pName):
                            self.probes.setMarker(pk, symbol, color, refresh=False)
                            self.probes.setRatio(pk, ratio, replot=False)
                if len(self.dataProbes) > 0:
                    self.probes.replotAllCurves(refresh=False)
                    self.probes.replotNormCurve(refresh=True)
            else:
                print "Warning: probe data must consist of attributes 'Probe ID', 'Ratio', 'ColorRGB' and 'Symbol'; 'Probe Name' is optional"


    def sendProbes(self):
        """Sends out example table with currently used probes.
        """
        if D1 or D2: print "OWNormalize.sendProbes"
        if self.data:
            vals = []
            pKeys = self.probes.keys()
            pKeys.sort()
            for pKey in pKeys:
                probe = self.probes[pKey]
                vals.append([probe.ID, probe.name, probe.ratioExpr, probe.getColorRgbStr(), probe.symbol])
            # create domain, clone and rename variables to "Probe ID" and "Probe Name"
            if self.varNameName <> "<none>":
                varName = self.varsAll[self.varNameName].clone()
                varName.name = "Probe Name"
            else:
                varName = orange.StringVariable("Probe Name")
            domain = orange.Domain([self.varsAll[self.varNameID].clone(), varName, orange.StringVariable("Ratio"), orange.StringVariable("ColorRGB"), orange.FloatVariable("Symbol")], None)
            domain[0].name = "Probe ID"
            et = orange.ExampleTable(domain, vals)
            # remove 'Probe Name' if <none>
            if self.varNameName == "<none>":
                et = orange.ExampleTable(orange.Domain([et.domain[0]] + et.domain[2:], None), et)
            self.send("Probes", et)
        else:
            self.send("Probes", None)


    ###################################################################################
    ## SELECTION OF VARIABLES
    ###################################################################################

    def varIDChange(self):
        """Refresh listbox containing other variables and refill self.tblControls.
        """
        if D1: print "OWNormalize.varIDChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.initProbes()
        self.setNormalization()
        # send data
        self.sendData()
        self.sendProbes()
        qApp.restoreOverrideCursor()


    def varDataChange(self, memberVarName):
        """Refresh listbox containing other variables.
        """
        if D1: print "OWNormalize.varDataChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.updateProbeData(self.data, self.varNameSignalSmpl, self.varNameSignalRef, self.varNameBGSmpl, self.varNameBGRef, self.varNameBGSmplSD, self.varNameBGRefSD, replot=True)
        # update info
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
        self.sendData()
        qApp.restoreOverrideCursor()


    def varSDChange(self, memberVarName):
        """Refresh listbox containing other variables; enables/disables Subtract BG checkbox and Max. CV slider.
        """
        if D1: print "OWNormalize.varSDChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        # enable/disable subtract BG checkbox & Max. CV slider
##        self.cbSubtrBG.setEnabled(self.varNameBGSmplSD != "<none>" and self.varNameBGRefSD != "<none>")
        self.boxMaxCV.setEnabled(self.varNameBGSmplSD != "<none>" and self.varNameBGRefSD != "<none>")
        # update
        self.probes.updateProbeData(self.data, self.varNameSignalSmpl, self.varNameSignalRef, self.varNameBGSmpl, self.varNameBGRef, self.varNameBGSmplSD, self.varNameBGRefSD, replot=True)
        # update info
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
        self.sendData()
        qApp.restoreOverrideCursor()


    def varOthersChange(self):
        """Updates list of selected other vars (lbVarOthers -> self.varsOtherSelected).
        """
        self.fillVarsOtherSelected()
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
        if D1 or D2: print "OWNormalize.fillProbeTable"
        if self.probes:
            # table row items and their number, header labels, len of sorting keys where spaces are added in front
            allProbes = self.probes.values()
            numProbes = len(allProbes)
            self.tblControls.setNumRows(numProbes)
            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcID, self.varNameID)
            sortKeyLen = int(math.log(numProbes, 10))+1
            if self.varNameName <> "<none>":
                # generate sorting keys: (ID, name)
                idIdx = zip(map(lambda pr: (pr.ID, pr.name), allProbes), range(numProbes))
                idIdx.sort()
                idRank = dict(zip([x[1] for x in idIdx], range(numProbes)))
                nameIdx = zip(map(lambda pr: (pr.name, pr.ID), allProbes), range(numProbes))
                nameIdx.sort()
                nameRank = dict(zip([x[1] for x in nameIdx], range(numProbes)))
                # fill rows
                for row, (key,pr) in enumerate(self.probes.items()):
                    # store the table row index of the probe
                    pr.tblRowIdx = row
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcID, str(pr.ID), editType=QTableItem.Never, sortingKey=self.sortingKey(idRank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcName, str(pr.name), editType=QTableItem.Never, sortingKey=self.sortingKey(nameRank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcPKey, key)
                    self.fillProbeTableItemMarkerRatio(row, pr)
                # show, label and adjust column tcName
                self.tblControls.showColumn(OWNormalize.tcName)
                self.tblControls.horizontalHeader().setLabel(OWNormalize.tcName, self.varNameName)
                # adjust columns' width
                self.tblControls.adjustColumn(OWNormalize.tcRatio)
                self.tblControls.adjustColumn(OWNormalize.tcID)
                self.tblControls.adjustColumn(OWNormalize.tcName)
                self.tblControls.setColumnWidth(OWNormalize.tcName, max(self.tblControls.columnWidth(OWNormalize.tcName), self.tblControls.visibleWidth() - self.tblControls.columnWidth(OWNormalize.tcRatio) - self.tblControls.columnWidth(OWNormalize.tcMarker) - self.tblControls.columnWidth(OWNormalize.tcID)))
            else:
                # generate sorting keys: ID
                idIdx = zip(map(lambda pr: pr.ID, allProbes), range(numProbes))
                idIdx.sort()
                idRank = dict(zip([x[1] for x in idIdx], range(numProbes)))
                nameRank = range(numProbes)
                # fill rows
                for row, (key,pr) in enumerate(self.probes.items()):
                    # store the table row index of the probe
                    pr.tblRowIdx = row
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcID, str(pr.ID), editType=QTableItem.Never, sortingKey=self.sortingKey(idRank[row], 20))#, background=QColor(160,160,160))
                    OWGUI.tableItem(self.tblControls, row, OWNormalize.tcPKey, key)
                    self.fillProbeTableItemMarkerRatio(row, pr)
                # hide column tcName
                self.tblControls.hideColumn(OWNormalize.tcName)
                # adjust columns' width
                self.tblControls.adjustColumn(OWNormalize.tcRatio)
                self.tblControls.adjustColumn(OWNormalize.tcID)
                self.tblControls.setColumnWidth(OWNormalize.tcID, max(self.tblControls.columnWidth(OWNormalize.tcID), self.tblControls.visibleWidth() - self.tblControls.columnWidth(OWNormalize.tcRatio) - self.tblControls.columnWidth(OWNormalize.tcMarker)))
        else:
            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcID, "Probe ID")
            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcName, "Probe Name")
            self.tblControls.setNumRows(0)

##    def fillProbeTable(self):
##        # init self.tblControls
##        if D1 or D2: print "OWNormalize.fillProbeTable"
##        if self.probes:
##            # table row items and their number, header labels, len of sorting keys where spaces are added in front
##            allProbes = self.probes.values()
##            numProbes = len(allProbes)
##            self.tblControls.setNumRows(numProbes)
##            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcID, self.varNameID)
##            sortKeyLen = int(math.log(numProbes, 10))+1
##            # generate sorting keys
##            idIdx = zip(map(lambda pr: pr.ID, allProbes), range(numProbes))
##            idIdx.sort()
##            idRank = dict(zip([x[1] for x in idIdx], range(numProbes)))
##            if self.varNameName <> "<none>":
##                nameIdx = zip(map(lambda pr: pr.name, allProbes), range(numProbes))
##                nameIdx.sort()
##                nameRank = dict(zip([x[1] for x in nameIdx], range(numProbes)))
##            else:
##                nameRank = range(numProbes)
##            # fill rows
##            for row, (key,pr) in enumerate(self.probes.items()):
##                # store the table row index of the probe
##                pr.tblRowIdx = row
##                OWGUI.tableItem(self.tblControls, row, OWNormalize.tcID, str(pr.ID), editType=QTableItem.Never, sortingKey=self.sortingKey(idRank[row], 20))#, background=QColor(160,160,160))
##                OWGUI.tableItem(self.tblControls, row, OWNormalize.tcName, str(pr.name), editType=QTableItem.Never, sortingKey=self.sortingKey(nameRank[row], 20))#, background=QColor(160,160,160))
##                OWGUI.tableItem(self.tblControls, row, OWNormalize.tcPKey, key)
##                self.fillProbeTableItemMarkerRatio(row, pr)
##            # adjust columns' width
##            if self.varNameName <> "<none>":
##                # show, label and adjust column tcName
##                self.tblControls.showColumn(OWNormalize.tcName)
##                self.tblControls.horizontalHeader().setLabel(OWNormalize.tcName, self.varNameName)
##                # adjust columns' width
##                self.tblControls.adjustColumn(OWNormalize.tcRatio)
##                self.tblControls.adjustColumn(OWNormalize.tcID)
##                self.tblControls.adjustColumn(OWNormalize.tcName)
##                self.tblControls.setColumnWidth(OWNormalize.tcName, max(self.tblControls.columnWidth(OWNormalize.tcName), self.tblControls.visibleWidth() - self.tblControls.columnWidth(OWNormalize.tcRatio) - self.tblControls.columnWidth(OWNormalize.tcMarker) - self.tblControls.columnWidth(OWNormalize.tcID)))
##            else:
##                # hide column tcName
##                self.tblControls.hideColumn(OWNormalize.tcName)
##                # adjust columns' width
##                self.tblControls.adjustColumn(OWNormalize.tcRatio)
##                self.tblControls.adjustColumn(OWNormalize.tcID)
##                self.tblControls.setColumnWidth(OWNormalize.tcID, max(self.tblControls.columnWidth(OWNormalize.tcID), self.tblControls.visibleWidth() - self.tblControls.columnWidth(OWNormalize.tcRatio) - self.tblControls.columnWidth(OWNormalize.tcMarker)))
##        else:
##            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcID, "Probe ID")
##            self.tblControls.horizontalHeader().setLabel(OWNormalize.tcName, "Probe Name")
##            self.tblControls.setNumRows(0)


    def fillProbeTableItemMarkerRatio(self, row, probe):
        """Updates a given row of self.tblControls: marker & ratio.
        """
        pxm = probe.getSymbolPixmap(self.tblControls.cellGeometry(row, OWNormalize.tcMarker))
        ratioSortKey = self.probes.getRatioSortingKey(probe.pKey)
##        if probe.color <> ProbeSet.NoColor:
        ch = probe.color.hsv()
        markerSortKey = "%02d%03d%03d%03d" % (probe.symbol, ch[0], ch[1], ch[2])
##        else:
##            markerSortKey = "%02d" % probe.symbol
        OWGUI.tableItem(self.tblControls, row, OWNormalize.tcRatio, self.probes.getRatioStr(probe.pKey), editType=QTableItem.OnTyping, sortingKey=ratioSortKey+markerSortKey) #, background=QColor(160,160,160))
        OWGUI.tableItem(self.tblControls, row, OWNormalize.tcMarker, "", editType=QTableItem.Never, sortingKey=markerSortKey+ratioSortKey, pixmap=pxm)#, background=QColor(160,160,160))


    def sortingKey(self, val, len):
        """Returns a string with leading spaces followed by str(val), whose length is at least len.
        """
        s = "%" + str(int(len)) + "." + str(int(len/2)) + "f"
        return s % val


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
        updates self.internContRatios and sends out data and control/ratios.
        """
        if D1 or D3: print "OWNormalize.tblControlsValueChange"
        pKey = str(self.tblControls.item(row, OWNormalize.tcPKey).text())
        ratio = str(self.tblControls.item(row, OWNormalize.tcRatio).text())
        if ratio <> self.probes[pKey].ratioExpr:
            # set ratio
            newRatio = self.probes.setRatio(pKey, ratio, replot=True)
            probe = self.probes[pKey]
########            # if ratio and no marker: set marker
########            if newRatio:
########                if probe.symbol == ProbeSet.NoSymbol:
########                    symbol = self.getProbeSymbolIdx()
########                else:
########                    symbol = probe.symbol
########                if probe.color == ProbeSet.NoColor:
########                    color = self.probeColor
########                else:
########                    color = probe.color
########                self.probes.setMarker(pKey, symbol, color, refresh=True)
            # update table & info
            self.fillProbeTableItemMarkerRatio(row, probe)
            self.setInfoProbes()
            self.setInfoFilterMaxCV()
            self.setInfoFilterMinRatio()
            self.setInfoFilterMaxInt()
##            self.setInfoFiltersAll()
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
##            self.probes.activate(pKey, True)


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
            if probe.symbol <> ProbeSet.NoSymbol:
                self.setProbeSymbolIdx(probe.symbol)
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
            if string.lower(self.controlName) in string.lower(self.tblControls.item(idx, OWNormalize.tcID).text()):
                sel = QTableSelection()
                sel.init(idx,OWNormalize.tcID)
                sel.expandTo(idx,OWNormalize.tcID)
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
        sel.expandTo(self.tblControls.numRows()-1,OWNormalize.tcName)
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
        self.updateSelectedProbes(self.ratioStr, self.probeColor, self.getProbeSymbolIdx(), True, False, False)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
        self.sendData()
        self.sendProbes()
        qApp.restoreOverrideCursor()


    def probeColorClick(self):
        """Show color-selection window, update button color and color of the symbols of the selected probes.
        """
        if D1: print "OWNormalize.probeColorCick"
        probeColor = QColorDialog.getColor(self.probeColor, self)
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        if probeColor.isValid():
            self.probeColor = probeColor
        self.btnProbeColor.pixmap().fill(self.probeColor)
        self.btnProbeColor.repaint()
        self.updateSelectedProbes(self.ratioStr, self.probeColor, self.getProbeSymbolIdx(), False, True, False)
        self.sendData()
        self.sendProbes()
        qApp.restoreOverrideCursor()


    def cmbProbeSymbolActivated(self):
        """Update symbol for the selected probes.
        """
        if D1: print "OWNormalize.cmbProbeSymbolActivated"
        self.updateSelectedProbes(self.ratioStr, self.probeColor, self.getProbeSymbolIdx(), False, False, True)
        # update infos (if probe gets plotted)
        self.setInfoProbes()
        self.sendData()
        self.sendProbes()

        
    def btnSetProbesClick(self):
        """Sets ratio, color and symbol for the selected controls.
        """
        if D1: print "OWNormalize.btnSetProbesClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.updateSelectedProbes(self.ratioStr, self.probeColor, self.getProbeSymbolIdx(), True, True, True)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
        self.sendData()
        self.sendProbes()
        qApp.restoreOverrideCursor()


    def btnClearProbesClick(self):
        """Clears ratios for the selected controls
        """
        if D1: print "OWNormalize.btnClearProbesClick"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.updateSelectedProbes("", ProbeSet.NoColor, ProbeSet.NoSymbol, True, True, True)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
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
                    self.probes.setRatio(pKey, ratioStr, replot=True)
                if updateColor or updateSymbol:
                    if updateColor:
                        newColor = color
                    else:
                        newColor = probe.color
                    if updateSymbol:
                        newSymbol = symbol
                    else:
                        newSymbol = probe.symbol
                    self.probes.setMarker(pKey, newSymbol, newColor, refresh=True)
                # update table
                self.fillProbeTableItemMarkerRatio(row, probe)


    def activateSelectedProbes(self):
        """Activate currently selected probes;
        """
        if D1 or D2: print "OWNormalize.activateSelectedProbes"
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
##            self.lblProbeInfo.setText("%d probes, %d plotted, %d controls." % (self.probes.getNumProbes(), self.probes.getNumProbesPlotted_nonFiltered(), self.probes.getNumProbesControls_nonFiltered()))
##            self.lblProbeInfo.setText("In total:\t%d control and %d other probes.\nAccepted:\t%d control and %d other probes.\nPlotted:\t%d control and %d other probes." %
##                                      (self.probes.getNumProbesControls(), self.probes.getNumProbesOthers(),
##                                       self.probes.getNumProbesControls_nonFiltered(), self.probes.getNumProbesOthers_nonFiltered(),
##                                       self.probes.getNumProbesControls_nonFiltered_plotted(), self.probes.getNumProbesOthers_nonFiltered_plotted()))
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
        titles = {False: ["Ratio"]*2 + ["Average intensity"]*2, True: ["M: Log2 ratio"]*2 + ["A: Log2 average intensity"]*2}
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
        if D1: print "OWNormalize.onMouseMoved"
        (curveKey, distPoints, x, y, pointKey) = self.graph.closestCurve(e.x(), e.y())
        curve = self.graph.curve(curveKey) 
        if curve and curve.__dict__.has_key("key"):
            pKey = curve.key
            probe = self.probes.get(pKey)
            if probe and distPoints<=self.markerSize/2:
                # activata probe (if not already active)
                self.probes.setCurveActiveList([pKey], refresh=True)
                # show tooltip
                xPoints = self.graph.transform(QwtPlot.xBottom, x)
                yPoints = self.graph.transform(QwtPlot.yLeft, y)
                rect = QRect(xPoints+self.graph.canvas().frameGeometry().x()-self.markerSize/2, yPoints+self.graph.canvas().frameGeometry().y()-self.markerSize/2, self.markerSize, self.markerSize)
                MyQToolTip.setRect(self.graph.tooltip, rect, str(probe.ID) + ", " + str(probe.name) + "\nSmpl (signal - bg = net) / Ref (signal - bg = net)\n" + self.probes.getDataTooltipStr(pKey))
            else:
                self.probes.setCurveActiveList(None, refresh=True)


    def onMousePressed(self, e):
        """If left button is pressed, the active control is made current in self.tblControls;
        first, OWGraph.onMousePressed(e) is executed (by default), followed by this code.
        """
        if D1: print "OWNormalize.onMousePressed"
        aProbes = self.probes.getActiveProbes()
        if e.button() == Qt.LeftButton and len(aProbes) > 0:
            # clear selections from tblControls (without fireing events)
            self.disconnect(self.tblControls , SIGNAL('selectionChanged()'), self.tblControlsSelectionChanged)
            self.tblControls.clearSelection(True)
            self.connect(self.tblControls , SIGNAL('selectionChanged()'), self.tblControlsSelectionChanged)
            # set current cell to the active probe (hopefully there is only one active)
            self.tblControls.setCurrentCell(aProbes[aProbes.keys()[0]].tblRowIdx, OWNormalize.tcID)


    ###################################################################################
    ## SETTINGS
    ###################################################################################

    def settingsSubstrBGChange(self):
        if D1: print "OWNormalize.settingsSubstrBGChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.setSubtrBG(self.subtrBG, replot=True)
        self.setInfoProbes()
        self.setInfoFilterMaxCV()
        self.setInfoFilterMinRatio()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
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
##        self.setInfoFiltersAll()
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
##        self.setInfoFiltersAll()
        self.sendData()
        qApp.restoreOverrideCursor()


    def settingsFilterMaxIntChange(self):
        """Handles changes of filter settings, which affects graph curves and ouput data.
        """
        if D1: print "OWNormalize.settingsFiltersChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.setProbeFilters()
        self.setInfoProbes()
        self.setInfoFilterMaxInt()
##        self.setInfoFiltersAll()
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
        if self.useMaxIntensity:
            maxIntensity = self.maxIntensity
        else:
            maxIntensity = Probes.bigVal
        self.probes.setFilterParameters(maxCV, minIntensityRatio, maxIntensity, replot=True)


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


    def setInfoFilterMaxInt(self):
        if D4: print "OWNormalize.setInfoFilterMaxInt"
        if self.probes:
            if self.probes.getNumProbesControls() > 0:
                ratioControls = 100. * self.probes.getNumFilteredControlsMaxInt() / self.probes.getNumProbesControls()
            else:
                ratioControls = 0
            if self.probes.getNumProbesOthers() > 0:
                ratioOthers = 100. * self.probes.getNumFilteredOthersMaxInt() / self.probes.getNumProbesOthers()
            else:
                ratioOthers = 0
            self.lblInfoFilterMaxInt.setText("%d (%.2f%s) control probes removed.\n%d (%.2f%s) other probes removed." % (self.probes.getNumFilteredControlsMaxInt(), ratioControls, "%", self.probes.getNumFilteredOthersMaxInt(), ratioOthers, "%"))
        else:
            self.lblInfoFilterMaxInt.setText("No data on input.\n")


##    def setInfoFiltersAll(self):
##        if D4: print "OWNormalize.setInfoFiltersAll"
##        if self.probes:
##            numC = self.probes.getNumProbesControls()
##            numC_filtered = self.probes.getNumProbesControls_filtered()
##            if numC > 0:
##                ratioControls = 100. * numC_filtered / numC
##            else:
##                ratioControls = 0
##            numO = self.probes.getNumProbesOthers()
##            numO_filtered = self.probes.getNumProbesOthers_filtered()
##            if self.probes.getNumProbesOthers() > 0:
##                ratioOthers = 100. * numO_filtered / numO
##            else:
##                ratioOthers = 0
##            self.lblInfoFiltersAll.setText("Total of %d control and %d other probes.\n%d (%.2f%s) control probes removed.\n%d (%.2f%s) other probes removed." % (numC, numO, numC_filtered, ratioControls, "%", numO_filtered, ratioOthers, "%"))
##        else:
##            self.lblInfoFiltersAll.setText("No data on input.\n\n")


    def settingsNormalizationChange(self):
        """Handles changes of normalization type, which affects normalization curve and output data.
        """
        if D1: print "OWNormalize.settingsNormalizationChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.setNormalization()
        self.sendData()
        qApp.restoreOverrideCursor()


    def setNormalization(self):        
        if D1: print "OWNormalize.setNormalization"
        if self.varNameName == "<none>":
            self.normRange = 0
        self.boxNormRange.setEnabled(self.varNameName <> "<none>")
        if self.probes:
            self.probes.setNormalizationParameters(self.normRange, self.normType, self.loessWindow, self.loessWeight, replot=True)


    def settingsReplicasChange(self):
        """Handles changes of replicas settings, which affects output data.
        """
        if D1: print "OWNormalize.settingsReplicasChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.rbgMergeOtherType.setEnabled(self.mergeReplType)
        self.sendData()
        qApp.restoreOverrideCursor()


    def settingsGraphAxisChange(self, axis):
        """Handles changes of graph axis settings; replot axes and curves.
        """
        if D1: print "OWNormalize.settingsGraphAxisChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.setPlotParameters(self.logAxisX, self.logAxisY, self.markerSize, replot=True)
        self.setGraphAxes([axis])
        qApp.restoreOverrideCursor()


    def settingsGraphChange(self):
        """Handles changes of graph settings; replot curves.
        """
        if D1: print "OWNormalize.settingsGraphChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.probes.setPlotParameters(self.logAxisX, self.logAxisY, self.markerSize, replot=True)
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


    def settingsOutputChange(self):    
        """Handles changes of output settings; send out data.
        """
        if D1: print "OWNormalize.settingsOutputChange"
        qApp.restoreOverrideCursor()
        qApp.setOverrideCursor(QWidget.waitCursor)
        self.sendData()
        # D1
        for probe in self.probes.values():
            pass
        qApp.restoreOverrideCursor()


    ###################################################################################
    ## UTILITY FUNCTIONS
    ###################################################################################

    def getProbeSymbolIdx(self):
        return self._probeSymbolComboIdx + 1


    def setProbeSymbolIdx(self, idx):
        self._probeSymbolComboIdx = idx - 1




class QwtPlotCurveKey(QwtPlotCurve):
    """QwtPlotCurve with additional member variable: key.
    """

    def __init__(self, parent, name, key):
        self.key = key
        QwtPlotCurve.__init__(self, parent, name)



class OWGraphMA(OWGraph):
    """OWGraph for MA plot with curves of type QwtPlotCurveKey;
    additionally handles zooming and mouse clicks.
    """

    def insertCurve(self, title, key, xAxis=QwtPlot.xBottom, yAxis=QwtPlot.yLeft):
        curve = QwtPlotCurveKey(self, title, key)
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





class ProbeSet:
    """A set of equivalent probes, their intensities, and their expected ratio, display symbol & color;
    curve represents normalization factors where color and symbol is used to plot the curve on the graph.
    """
    _activePenWidth = 3
    _inactivePenWidth = 1
    NoSymbol = 0
    NoColor = QColor(0,0,0)
    
    def __init__(self, ID, name, pKey):
        if D1: print "ProbeSet.__init__"
        self.ID = ID                # string
        self.name = name            # string
        self.pKey = pKey            # ID / ID+name
        # ratio
        self.ratioExpr = ""         # string from which ratio is evaluated
        # marker
        self.color = ProbeSet.NoColor           # QColor
        self.symbol = ProbeSet.NoSymbol         # int (QwtSymbol.Style)
##        self._penColor = QColor(255,255,255)    # color of symbols from table
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


    def addProbeData(self, dataIndices):
        if D1: print "ProbeSet.addProbeData"
        self._dataIndices.extend(dataIndices)


    def addProbeDatum(self, dataIdx):
        if D1: print "ProbeSet.addProbeDatum"
        self._dataIndices.append(dataIdx)


    def clearProbeData(self):
        if D1: print "ProbeSet.clearProbeData"
        self._dataIndices = []


    ############################################
    # TABLE
    ############################################

    def getSymbolPixmap(self, rect):
        """Input: QRect instance; output: QPixmap instance.
        """
        if D1: print "ProbeSet.getSymbolPixmap"
        # init pixmap for table item
        symbolPixmap = QPixmap(rect.width(),rect.height())
        symbolPixmap.fill(QColor(255,255,255))
        if self.symbol <> ProbeSet.NoSymbol:
            painter = QPainter(symbolPixmap)
            symbol = QwtSymbol(self.symbol, QBrush(self.color, QBrush.SolidPattern), QPen(QColor(0,0,0),1), QSize(8,8))
            symbol.draw(painter, QPoint(rect.width()/2,rect.height()/2))
            painter.end()
        return symbolPixmap


    def getColorRgbStr(self):
        """returns hexadecimal color string, e.g. "00ff10
        """
        return string.replace("%2s%2s%2s" % (hex(self.color.red())[2:], hex(self.color.green())[2:], hex(self.color.blue())[2:]), " ", "0")




class Probes(dict):
    """Dictionary of ProbeSet items.
    replot / refresh: recalculate MA data / 
    """

    """TODO:
    - setFilterParameters: divide into three functions, one for each parameter
    """

    mergeFunctions = {1:MA.average, 2:NumExtn.medianMA}
    bigVal = 1e20
    midVal = 1e10

    def __init__(self, graph, subtrBG, logAxisX, logAxisY, markerSize):
        if D1 or D2: print "Probes.__init__"
##        self._plotted = {}  # if probe.curve <> None
        self._active = {}   # currently selected probe sets
        self.graph = graph
        # plot parameters
        self.subtrBG = subtrBG
        self.logAxisX = logAxisX
        self.logAxisY = logAxisY
        self.markerSize = markerSize
        # var names
        self.varNameID = None
        self.varNameName = "<none>"
        # IDs, names
        self._IDList = []
        self._nameList = []
        # data: Numeric arrays
        self.__sigSmpl = None
        self.__sigRef = None
        self.__bgSmpl = None
        self.__bgRef = None
        self.__bgSmplSD = None
        self.__bgRefSD = None
        # ratios
        self.__ratio = None
        # Numeric array: 0: OK, 1: filtered out
        self.__filterMaxCV = None
        self.__filterMinRatio = None
        self.__filterMaxInt = None
        self._control = Numeric.zeros((0,), Numeric.Int)
        self.__plotted = Numeric.zeros((0,), Numeric.Int)
        # normalization functions
        self._normFuncDict = {0:self._getNormCurveMedian, 1:self._getNormCurveLinReg, 2:self._getNormCurveLoess}
        self._normFunction = self._normFuncDict[2]
        # normalization curves
        self._normCurve = None
        self._normCurveLow = None
        self._normCurveHigh = None
        # default parameters
        self.normRange = 0  # 0: global, 1: local, per probe names
        self.maxCV = Probes.bigVal
        self.minIntensityRatio = 0
        self.maxIntensity = Probes.bigVal
        self.loessWindow = 60
        self.loessWeight = 0
        # number of probes (controls, others) filtered by different filters


##    def getKeysStartsWith(self, pKeySub):
##        """Return a list of probe keys which contain a given subkey.
##        """
##        if D1: print "Probes.getKeyStartsWith"
##        pKeysStartsWith = []
##        for pKey in dict.keys(self):
##            if pKey.startswith(pKeySub):
##                pKeysStartsWith.append(pKey)
##        return pKeysStartsWith

    def getKeysFromIDsNames(self, IDSub, nameSub):
        """Return a list of probe keys which contain a given part of ID and name (optionally);
        """
        if D1: print "Probes.getKeysFromIDsNames"
        probesIdName = {}   # probes where both ID and name matches
        probesId = {}       # probes where only ID matches
##        probesName = {}     # probes where only name matches
        IDSub = str(IDSub)
        nameSub = str(nameSub)
        for pKey, probe in self.items():
            if IDSub and (IDSub in str(probe.ID)):
                probesId[pKey] = probe
                if nameSub and (nameSub in str(probe.name)):
                    probesIdName[pKey] = probe
##            else:
##                if nameSub and nameSub in probe.name:
##                    probesName[pKey] = probe
##        if len(probesIdName) > 0:
##            return probesIdName.keys()
##        elif len(probesId) > 0:
##            return probesId.keys()
##        else:
##            return probesName.keys()
        if len(probesIdName) > 0:
            return probesIdName.keys()
        else:
            return probesId.keys()


    def setSubtrBG(self, subtrBG, replot=True):
        if D1 or D2: print "Probes.setSubtrBG"
        if self.subtrBG <> subtrBG:
            self.subtrBG = subtrBG
            if replot:
                self.replotAllCurves(False)
                self.replotNormCurve(True)
        

    def setPlotParameters(self, logAxisX, logAxisY, markerSize, replot=True):
        if D1 or D2: print "Probes.setPlotParameters"
        if self.logAxisX <> logAxisX or  self.logAxisY <> logAxisY or self.markerSize <> markerSize:
            self.logAxisX = logAxisX
            self.logAxisY = logAxisY
            self.markerSize = markerSize
            if replot:
                self.replotAllCurves(False)
                self.replotNormCurve(True)


    def setFilterParameters(self, maxCV, minIntensityRatio, maxIntensity, replot=True):
        if D1 or D2: print "Probes.setFilterParameters"
        if maxCV <> self.maxCV or minIntensityRatio <> self.minIntensityRatio or maxIntensity <> self.maxIntensity:
            self.maxCV = maxCV                        
            self.minIntensityRatio = minIntensityRatio
            self.maxIntensity = maxIntensity
            self.__filterMaxCV = None
            self.__filterMinRatio = None
            self.__filterMaxInt = None
            if replot:
                self.replotAllCurves(False)
                self.replotNormCurve(True)


    def setNormalizationParameters(self, normRange, normType, loessWindow, loessWeight, replot=True):
        """ normRange: 0: global, 1: per probe name
            normType: 0: median, 1: LR, 2: LOESS
        """
        if D1 or D2: print "Probes.setNormalizationParameters"
        if self._normFuncDict[normType] <> self._normFunction or loessWindow <> self.loessWindow or loessWeight <> self.loessWeight or normRange <> self._normRange:
            self._normRange = normRange
            self._normFunction = self._normFuncDict[normType]
            self.loessWindow = loessWindow
            self.loessWeight = loessWeight
            if replot:
                self.replotNormCurve(True)


    def initProbes(self, data, varNameID, varNameName, varNameSignalSmpl, varNameSignalRef, varNameBGSmpl, varNameBGRef, varNameBGSmplSD, varNameBGRefSD):
        """Input: orange.ExampleTable, orange.Variable.names;
        stores self.varNameID and self.varNameName.
        """
        if D1 or D2: print "Probes.initProbes"
        self.clear(refresh=True)
        self.varNameID = varNameID
        self.varNameName = varNameName
        self._control = Numeric.zeros((len(data),), Numeric.Int)
        self.__plotted = Numeric.zeros((len(data),), Numeric.Int)
        # update data and probe data indices
        self.updateProbeData(data, varNameSignalSmpl, varNameSignalRef, varNameBGSmpl, varNameBGRef, varNameBGSmplSD, varNameBGRefSD, replot=False)
        if self.varNameName <> "<none>":
            for eIdx, e in enumerate(data):
                self._addProbeData(str(e[varNameID].native())+str(e[varNameName].native()), e[varNameID].native(), e[varNameName].native(), eIdx)
        else:
            for eIdx, e in enumerate(data):
##                self._addProbeData(str(e[varNameID].native()), e[varNameID].native(), None, eIdx)
                self._addProbeData(str(e[varNameID].native()), e[varNameID].native(), "", eIdx)


    def clear(self, refresh=True):
        if D1 or D2: print "Probes.clear"
        self.removeNormCurves(False)
        self.removeAllCurves(refresh)
        dict.clear(self)
##        self._plotted = {}
        self._active = {}
        self._IDList = []
        self._nameList = []
        self.__sigSmpl = None
        self.__sigRef = None
        self.__bgSmpl = None
        self.__bgRef = None
        self.__bgSmplSD = None
        self.__bgRefSD = None
        # ratio
        self.__ratio = None
        # Numeric array: 0: OK, 1: filtered out
        self.__filterMaxCV = None
        self.__filterMinRatio = None
        self.__filterMaxInt = None
        self._control = Numeric.zeros((0,), Numeric.Int)
        self.__plotted = Numeric.zeros((0,), Numeric.Int)


    def updateProbeData(self, data, varNameSignalSmpl, varNameSignalRef, varNameBGSmpl, varNameBGRef, varNameBGSmplSD, varNameBGRefSD, replot=True):
        """Update signal and background of the selected probes, construct new filter.
        """
        if D1 or D2: print "Probes.updateProbeData"
        if data:
            dataMA = data.toMA("a")[0]
            self.__sigSmpl = dataMA[:,data.domain.index(varNameSignalSmpl)]
            self.__sigRef = dataMA[:,data.domain.index(varNameSignalRef)]
            self.__bgSmpl = dataMA[:,data.domain.index(varNameBGSmpl)]
            self.__bgRef = dataMA[:,data.domain.index(varNameBGRef)]
            self.__ratio = Numeric.ones((len(data),), Numeric.Float)
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
            self.__filterMaxInt = None
            if replot:
                self.replotAllCurves(False)
                self.replotNormCurve(True)


    def _addProbeData(self, pKey, ID, name, dataIdx):
        if D1: print "Probes._addProbeData"
        self._IDList.append(ID)
        self._nameList.append(name)
        if dict.has_key(self, pKey):
            ps = dict.get(self, pKey)
        else:
            ps = ProbeSet(ID, name, pKey)
            dict.__setitem__(self, pKey, ps)
        ps.addProbeDatum(dataIdx)

        
    ############################################
    # NUMBER OF PROBES
    ############################################

    def getNumReplicas_nonFiltered(self, merge):
        """Returns (..., 2) Numeric array where rows represent different probes and columns:
            0: number of all probes
            1: number of non-filtered probes
        """
        if merge:
            return self.__mergeReplicas(Numeric.transpose(Numeric.asarray([Numeric.ones(self.getFilter().shape), Numeric.logical_not(self.getFilter())])), Numeric.add.reduce)
        else:
            return Numeric.transpose(Numeric.asarray([Numeric.ones(self.getFilter().shape), Numeric.logical_not(self.getFilter())]))


    def getNumFilteredControlsMaxCV(self):
        if type(self.__filterMaxCV) == Types.NoneType:
            self._setFilterMaxCV()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxCV, self._control))

    def getNumFilteredOthersMaxCV(self):
        if type(self.__filterMaxCV) == Types.NoneType:
            self._setFilterMaxCV()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxCV, Numeric.logical_not(self._control)))

    def getNumFilteredControlsMinRatio(self):
        if type(self.__filterMinRatio) == Types.NoneType:
            self._setFilterMinRatio()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMinRatio, self._control))

    def getNumFilteredOthersMinRatio(self):
        if type(self.__filterMinRatio) == Types.NoneType:
            self._setFilterMinRatio()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMinRatio, Numeric.logical_not(self._control)))

    def getNumFilteredControlsMaxInt(self):
        if type(self.__filterMaxInt) == Types.NoneType:
            self._setFilterMaxInt()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxInt, self._control))
    
    def getNumFilteredOthersMaxInt(self):
        if type(self.__filterMaxInt) == Types.NoneType:
            self._setFilterMaxInt()
        return Numeric.add.reduce(Numeric.logical_and(self.__filterMaxInt, Numeric.logical_not(self._control)))


    def getNumProbes(self):
        return self._control.shape[0]

    def getNumProbesControls(self):
        return Numeric.add.reduce(Numeric.greater(self._control, 0))

    def getNumProbesOthers(self):
        return Numeric.add.reduce(Numeric.logical_not(self._control))


    def getNumProbesControls_filtered(self):
        return Numeric.add.reduce(Numeric.logical_and(self._control, self.getFilter()))

    def getNumProbesControls_nonFiltered(self):
        return Numeric.add.reduce(Numeric.logical_and(self._control, Numeric.logical_not(self.getFilter())))

    def getNumProbesOthers_filtered(self):
        return Numeric.add.reduce(Numeric.logical_and(Numeric.logical_not(self._control), self.getFilter()))

    def getNumProbesOthers_nonFiltered(self):
        return Numeric.add.reduce(Numeric.logical_and(Numeric.logical_not(self._control), Numeric.logical_not(self.getFilter())))

    def getNumProbes_nonFiltered_plotted(self):
        return Numeric.add.reduce(Numeric.logical_and(self.__plotted, Numeric.logical_not(self.getFilter())))

    def getNumProbesControls_nonFiltered_plotted(self):
        return Numeric.add.reduce(Numeric.logical_and(self._control, Numeric.logical_and(self.__plotted, Numeric.logical_not(self.getFilter()))))

    def getNumProbesOthers_nonFiltered_plotted(self):
        return Numeric.add.reduce(Numeric.logical_and(Numeric.logical_not(self._control), Numeric.logical_and(self.__plotted, Numeric.logical_not(self.getFilter()))))

    ############################################
    # MARKER, RATIO
    ############################################

    def setMarker(self, pKey, symbol, color, refresh=True):
        """Set symbol or color to None to remove probe set from the graph.
        """
        if D1 or D2 or D3: print "Probes.setMarker"
        if dict.has_key(self, pKey):
            ps = dict.get(self, pKey)
            ps.symbol = symbol
            ps.color = color
            self._replotCurve(ps, refresh)


    def setRatio(self, pKey, ratioExpr, replot=True):
        """Sets self.__ratio, self._control and dict[pKey].ratioExpr;
        ratioExpr should be None in order not to use probe set as control.
        """
        if D1 or D2 or D3: print "Probes.setRatio"
        if dict.has_key(self, pKey):
            ps = dict.get(self, pKey)
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
            if ps and ps.ratioExpr <> newRatioExpr:
                ps.ratioExpr = newRatioExpr
                Numeric.put(self.__ratio, ps.getDataIndices(), ratio)
                # if we succeed in setting the ratio: put probe among controls
                Numeric.put(self._control, ps.getDataIndices(), ps.ratioExpr <> "")
                if replot:                    
                    self._replotCurve(ps, False)
                    self.replotNormCurve(True)
        return ratio


    def getRatioStr(self, pKey):
        if D1: print "Probes.getRatioStr"
        probe = self.__getitem__(pKey)
        if probe and probe.ratioExpr and len(probe.getDataIndices())>0:
            return "%.2f" % self.__ratio[probe.getDataIndices()[0]]
        else:
            return ""


    def getRatioSortingKey(self, pKey):
        """Returns a string with leading spaces followed by str(val), whose length is at least len."""
        probe = self.__getitem__(pKey)
        if probe and probe.ratioExpr and len(probe.getDataIndices())>0:
            return "%15.7f" % self.__ratio[probe.getDataIndices()[0]]
        else:
            return ""


    def getActiveProbes(self):
        if D1: print "Probes.getActiveProbes"
        return self._active


##    def getPlotedProbes(self):
##        if D1: print "Probes.getPlotedProbes"
##        return self._plotted


    ############################################
    # GRAPH
    ############################################

    def _removeCurve(self, probe, refresh=True):
        if D1 or D2 or D3: print "Probes._removeCurve"
        if probe and probe.curve:
##            self._plotted.pop(pKey)
            Numeric.put(self.__plotted, probe.getDataIndices(), 0)
            self.graph.removeCurve(probe.curve)
            probe.curve = None
        if refresh: self.graph.replot()


    def removeAllCurves(self, refresh=True):
        if D1 or D2: print "Probes.removeAllCurves"
##        for probe in self._plotted.values():
        for probe in self.values():
            if probe.curve:
                self.graph.removeCurve(probe.curve)
                probe.curve = None
##        self._plotted = {}
        self.__plotted *= 0
        if refresh and len(self)>0: self.graph.replot()


    def _replotCurve(self, probe, refresh=True):
        if D1 or D2 or D3: print "Probes._replotCurve"
        change = False
        if probe.curve:
            self._removeCurve(probe, False)
            change = True
        if probe.symbol <> ProbeSet.NoSymbol:
            probe.curve = self.graph.insertCurve(probe.ID, probe.pKey)
            Numeric.put(self.__plotted, probe.getDataIndices(), 1)
            M,A = self.getMA(probe.pKey)
            self.graph.setCurveData(probe.curve, A, M)
            self.graph.setCurveStyle(probe.curve, QwtCurve.NoCurve)
            self._setCurveSymbol(probe, False)
##            self._plotted[pKey] = probe
            change = True
        if change and refresh: self.graph.replot()


    def _setCurveSymbol(self, probe, refresh=True):
        """sets graph marker symbol
        """
        if probe.curve:
            if self._active.has_key(probe.pKey):
                pen = QPen(QColor(0,0,0),ProbeSet._activePenWidth)
            else:
                pen = QPen(QColor(0,0,0),ProbeSet._inactivePenWidth)
            qSymbol = QwtSymbol(probe.symbol, QBrush(probe.color, QBrush.SolidPattern), pen, QSize(self.markerSize,self.markerSize))
            self.graph.setCurveSymbol(probe.curve, qSymbol)
            if refresh: self.graph.replot()


    def replotAllCurves(self, refresh=True):
        if D1 or D2: print "Probes.replotAllCurves"
##        for pKey in self._plotted.keys():
        for probe in self.values():
            self._replotCurve(probe, False)
        if refresh: self.graph.replot()


    def setCurveActive(self, pKey, active, refresh=True):
        if D1 or D3: print "Probes.setCurveActive"
        probe = self.__getitem__(pKey)
        if probe <> None and probe.curve <> None:
            if active <> self._active.has_key(pKey):
                if active:
                    self._active[pKey] = probe
                else:
                    self._active.pop(pKey)
                self._setCurveSymbol(probe, refresh)


    def switchCurveActive(self, pKey, refresh=True):
        if D1 or D3: print "Probes.switchCurveActive"
        self.setCurveActive(pKey, not(self._active.has_key(pKey)), refresh)


    def setCurveActiveList(self, pKeyList, refresh=True):
        """Deactivetes currently active, activates those from the given list;
        pKeyList : [ [pKey1, pKey2,...] | None ]
        """
        if D1 or D3: print "Probes.setCurveActiveList"
        if pKeyList:
            for pKey in pKeyList:
                if not self._active.has_key(pKey):
                    self.setCurveActive(pKey, True, False)
            for pKey in self._active.keys():
                if pKey not in pKeyList:
                    self.setCurveActive(pKey, False, False)
        else:
            for pKey in self._active.keys():
                self.setCurveActive(pKey, False, False)
        if refresh: self.graph.replot()


    ############################################
    # FILTER
    ############################################

    def getFilter(self):
        if D1: print "Probes.getFilter"
        if type(self.__filterMaxCV) == Types.NoneType:
            self._setFilterMaxCV()
        if type(self.__filterMinRatio) == Types.NoneType:
            self._setFilterMinRatio()
        if type(self.__filterMaxInt) == Types.NoneType:
            self._setFilterMaxInt()
        return Numeric.logical_or(Numeric.logical_or(self.__filterMaxCV, self.__filterMinRatio), self.__filterMaxInt)


    def _setFilterMaxCV(self):
        if D1 or D2 or D4: print "Probes._setFilterMaxCV"
        if self.__sigSmpl:
            # maxCV: bgSD / sig <= self.maxCV
            if self.__bgSmplSD <> None and self.__bgRefSD <> None:
                self.__filterMaxCV = MA.asarray(self.__bgSmplSD / self.__sigSmpl).filled(Probes.midVal) > self.maxCV
##                print Numeric.add.reduce(Numeric.take(Numeric.logical_not(self.__filter), self.get("ST_sbspotting_solution").getDataIndices())),
##                print MA.add.reduce(MA.take(MA.logical_not(self.__bgSmplSD / self.__sigSmpl), self.get("ST_sbspotting_solution").getDataIndices()))
                self.__filterMaxCV += MA.asarray(self.__bgRefSD / self.__sigRef).filled(Probes.midVal) > self.maxCV
##                print Numeric.add.reduce(Numeric.take(Numeric.logical_not(self.__filter), self.get("ST_sbspotting_solution").getDataIndices())),
##                print MA.add.reduce(MA.take(MA.logical_not(self.__bgRefSD / self.__sigRef), self.get("ST_sbspotting_solution").getDataIndices()))
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
##            print Numeric.add.reduce(Numeric.take(Numeric.logical_not(self.__filter), self.get("ST_sbspotting_solution").getDataIndices()))
##            print MA.add.reduce(MA.take(MA.logical_not(self.__bgRefSD / self.__sigRef), self.get("ST_sbspotting_solution").getDataIndices()))
            # convert to 0/1
            self.__filterMinRatio = self.__filterMinRatio > 0

    def _setFilterMaxInt(self):
        if D1 or D2 or D4: print "Probes._setFilterMaxInt"
        if self.__sigSmpl:
            # maxIntensity: sig <= maxIntensity
            self.__filterMaxInt = MA.asarray(self.__sigSmpl > self.maxIntensity).filled(1)
##            print Numeric.add.reduce(Numeric.take(Numeric.logical_not(self.__filter), self.get("ST_sbspotting_solution").getDataIndices()))
            self.__filterMaxInt += MA.asarray(self.__sigRef > self.maxIntensity).filled(1)
##            print Numeric.add.reduce(Numeric.take(Numeric.logical_not(self.__filter), self.get("ST_sbspotting_solution").getDataIndices()))
            # convert to 0/1
            self.__filterMaxInt = self.__filterMaxInt > 0
##            print Numeric.add.reduce(Numeric.take(Numeric.logical_not(self.__filter), self.get("ST_sbspotting_solution").getDataIndices()))


    ############################################
    # DATA (accounts for filters)
    # _get..._masked: returns MA array
    # _get...: returns compressed Numeric array
    ############################################

    def _sigSmpl_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__sigSmpl)

    def _sigRef_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__sigRef)

    def _bgSmpl_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__bgSmpl)

    def _bgRef_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__bgRef)

    def _ratio_masked(self, condition):
        return MA.masked_where(Numeric.logical_or(Numeric.logical_not(condition), self.getFilter()), self.__ratio)

    def _getNet_masked(self, condition):
        netSmpl = self._sigSmpl_masked(condition)
        netRef =  self._sigRef_masked(condition)
        if self.subtrBG:
            netSmpl -= self._bgSmpl_masked(condition)
            netRef -=  self._bgRef_masked(condition)
        return netSmpl, netRef



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

    def _getNet(self, condition):
        netSmpl, netRef = self._getNet_masked(condition)
        return netSmpl.compressed(), netRef.compressed()



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



    def __getM_masked(self, condition, netSmpl, netRef):
        M = netSmpl / netRef / self._ratio_masked(condition)
        if self.logAxisY:
            M = MA.log(M) / math.log(2)
        return M

    def _getM_masked(self, condition):
        netSmpl, netRef = self._getNet_masked(condition)
        return self.__getM_masked(condition, netSmpl, netRef)

    def _getM(self, condition):
        return Numeric.asarray(self._getM_masked(condition).compressed())



    def __getA_masked(self, condition, netSmpl, netRef):
        A = MA.sqrt(netSmpl*netRef)
        if self.logAxisX:
            A = MA.log(A) / math.log(2)
        return A

    def _getA_masked(self, condition):
        netSmpl, netRef = self._getNet_masked(condition)
        return self.__getA_masked(condition, netSmpl, netRef)

    def _getA(self, condition):
        return Numeric.asarray(self._getA_masked(condition).compressed())



    def _getMA_masked(self, condition):
        """Returns MA arrays: M/ratio, A (compressed by filter and condition).
        """
        if D1 or D2: print "Probes._getMA_masked"
        netSmpl, netRef = self._getNet_masked(condition)
        return self.__getM_masked(condition, netSmpl, netRef), self.__getA_masked(condition, netSmpl, netRef)

    def _getMA(self, condition):
        """Returns Numeric arrays: M/ratio, A (compressed by filter, mask and condition).
        """
        M,A = self._getMA_masked(condition)
        noMask = Numeric.logical_not(Numeric.logical_or(MA.getmaskarray(A), MA.getmaskarray(M)))
        return Numeric.asarray(MA.compress(noMask, M)), Numeric.asarray(MA.compress(noMask, A))

    def getMA(self, pKey):
        """Returns Numeric arrays: M/ratio, A (compressed by filter, condition and mask)
        """
        if D1 or D2: print "Probes.getMA"
        cond = NumExtn.indices2condition(dict.__getitem__(self, pKey).getDataIndices(), self.__sigSmpl.shape[0])
        return self._getMA(cond)


    ############################################
    # TOOLTIP
    ############################################

    def getDataTooltipStr(self, pKey):
        if D1: print "Probes.getDataTooltipStr"
        if D4: print "# self.__plotted and self.getFilter(): ", Numeric.add.reduce(self.__plotted == self.getFilter())
        out = ""
        for ss,sr,bs,br in zip(self.sigSmpl(pKey), self.sigRef(pKey), self.bgSmpl(pKey), self.bgRef(pKey)):
            out += "%5.f - %5.f = %5.f  /  %5.f - %5.f = %5.f\n" % (ss,bs,ss-bs, sr,br,sr-br)
        return out[:-1]


    ############################################
    # NORMALIZATION CURVE
    ############################################

    def _getNormCurve(self, A):
        """Returns normalized M in points A;
        """
        Mc,Ac = self._getMA(self._control)
        if Mc.shape[0] > 1:
            return self._normFunction(A, Mc, Ac)
        else:
            return Numeric.array([])


    def _getNormCurveLoess(self, A, Mc, Ac):
        return Numeric.asarray(statc.loess(zip(Ac, Mc), Numeric.asarray(A).tolist(), self.loessWindow/100.))[:,1]


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
        if D4: print "Probes._getNormCurveMedian, value:", NumExtn.median(Mc)
        return Numeric.resize(NumExtn.median(Mc), A.shape)
        

    ############################################
    # GRAPH: NORM. CURVE
    ############################################

    def replotNormCurve(self, refresh=True):
        if D1 or D2: print "Probes.replotNormCurve"
        self.removeNormCurves(False)
        if self.__sigSmpl:
            A1 = self._getA(1)
            if A1.shape[0] > 1:
                minA1 = min(A1)
                maxA1 = max(A1)
                loessA = Numeric.arange(minA1, maxA1, (maxA1-minA1)/(10000./self.loessWindow))
                loessM = self._getNormCurve(loessA)
                if loessM.shape == loessA.shape:
                    # M,A from controls (to compute minAc & maxAc)
                    Ac = self._getA(self._control)
                    minAc = min(Ac)
                    maxAc = max(Ac)
                    # plot norm. curve
                    self._normCurve = self.graph.insertCurve("Norm. curve", None)
                    cond = Numeric.logical_and(Numeric.greater_equal(loessA, minAc), Numeric.less_equal(loessA, maxAc))
                    self.graph.setCurveData(self._normCurve, Numeric.compress(cond, loessA), Numeric.compress(cond, loessM))
                    pen = QPen(QColor(0,0,0),ProbeSet._activePenWidth)
                    self.graph.setCurvePen(self._normCurve, pen)
                    # plot extrapolated lower part of norm. curve
                    cond = Numeric.less_equal(loessA, minAc)
                    if Numeric.add.reduce(cond) > 0:
                        cond[Numeric.add.reduce(cond)] = 1  # add another point to connect with norm. curve
                        self._normCurveLow = self.graph.insertCurve("Norm. curve (low-int. extrap.)", None)
                        self.graph.setCurveData(self._normCurveLow, Numeric.compress(cond, loessA), Numeric.compress(cond, loessM))
                        pen = QPen(QColor(0,0,255),ProbeSet._activePenWidth)
                        self.graph.setCurvePen(self._normCurveLow, pen)
                    # plot extrapolated upper part of norm. curve
                    cond = Numeric.greater_equal(loessA, maxAc)
                    if Numeric.add.reduce(cond) > 0:
                        cond[cond.shape[0]-Numeric.add.reduce(cond)-1] = 1  # add another point to connect with norm. curve
                        self._normCurveHigh = self.graph.insertCurve("Norm. curve (high-int. extrap.)", None)
                        self.graph.setCurveData(self._normCurveHigh, Numeric.compress(cond, loessA), Numeric.compress(cond, loessM))
                        pen = QPen(QColor(0,255,0),ProbeSet._activePenWidth)
                        self.graph.setCurvePen(self._normCurveHigh, pen)
        if refresh: self.graph.replot()


    def removeNormCurves(self, refresh=True):
        if D1 or D2: print "Probes.removeNormCurves"
        change = False
        if self._normCurve <> None:
            self.graph.removeCurve(self._normCurve)
            self._normCurve = None
            change = True
        if self._normCurveLow <> None:
            self.graph.removeCurve(self._normCurveLow)
            self._normCurveLow = None
            change = True
        if self._normCurveHigh <> None:
            self.graph.removeCurve(self._normCurveHigh)
            self._normCurveHigh = None
            change = True
        if change and refresh:
            self.graph.replot()


    ############################################
    # MERGE REPLICAS
    ############################################

    def __mergeReplicas(self, ma, mergeFunction):
        shp = list(ma.shape)
        shp[0] = len(self.values())
        maMerged = MA.zeros(shp, ma.typecode())
        for idx, probe in enumerate(self.values()):
            maMerged[idx] = mergeFunction(MA.take(ma, probe.getDataIndices()))
        return maMerged


    ############################################
    # NORMALIZED DATA
    ############################################

    def getNormalizedLog2Ratio_masked(self, mergeType):
        """Returns masked array of log2ratio of individual probes;
        accounts for filters, but NOT for ratios;
        mergeType: 0:None, 1:mean, 2:median
        """
        M,A = self._getMA_masked(1)
        nonMaskedIndices = NumExtn.condition2indices(Numeric.logical_not(MA.getmaskarray(A)))
        normCurve = self._getNormCurve(MA.compress(Numeric.logical_not(MA.getmaskarray(M)), A))
        if normCurve.shape == M.compressed().shape:
            # put normalized M values back to masked array M
            if self.logAxisY:
                MA.put(M, nonMaskedIndices, M.compressed() - normCurve)
            else:
                MA.put(M, nonMaskedIndices, M.compressed() / normCurve)
                M = MA.log(M) / math.log(2)
            M += MA.log(self._ratio_masked(1)) / math.log(2)
            if mergeType:
                return self.__mergeReplicas(M, Probes.mergeFunctions[mergeType])
            else:
                return M
        else:
            if mergeType:
                return MA.zeros(len(self.values()), MA.Float) * MA.masked
            else:
                return MA.zeros(M.shape, MA.Float) * MA.masked


    def getNetIntensity_smpl_ref(self, mergeType):
        netS, netR = self._getNet_masked(1)
        iSR = MA.concatenate([MA.reshape(netS, (netS.shape[0], 1)), MA.reshape(netR, (netR.shape[0], 1))], 1)
        if mergeType:
            return self.__mergeReplicas(iSR, Probes.mergeFunctions[mergeType])
        else:
            return iSR
    

    def getRawLog2Ratio_masked(self, mergeType):
        l2r = MA.log(self.__sigSmpl / self.__sigRef) / math.log(2)
        if mergeType:
            return self.__mergeReplicas(l2r, Probes.mergeFunctions[mergeType])
        else:
            return l2r

    ############################################
    # IDs, NAMEs
    ############################################

    def getIDs(self, mergeType):
        if mergeType:
            return map(lambda x: x.ID, self.values())
        else:
            return self._IDList


    def getNames(self, mergeType):
        if mergeType:
            return map(lambda x: x.name, self.values())
        else:
            return self._nameList





if __name__=="__main__":
    import OWDataTable, orngSignalManager
    signalManager = orngSignalManager.SignalManager(0)
    a=QApplication(sys.argv)
    ow=OWNormalize(signalManager = signalManager)
    a.setMainWidget(ow)
    ow.show()
##    ow.onDataInput(orange.ExampleTable(r"C:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\Tadeja 2nd image analysis\10vs10mg original data\0449.txt", DC="<NO DATA>"))
##    ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\Tadeja drago\05vs10mg\chol.diet\2537.txt"))
    ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\Tadeja drago\05vs10mg\control\6033.txt"))
    ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.0 mouse\sterolgene v.0 mouse controlGeneRatios 2.tab"))

##    ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 human\Krisztina\13217291-A01.txt", DC="<NO DATA>", noClass=1, noCodedDiscrete=1))
##    ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\STEROLTALK\Sterolgene v.1 human\Sterolgene v1 ControlGeneRatios REVERSED 2.tab"))

##    ow.onDataInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\MF other\Patologija\2005-11-21\136926 results.txt", DC="<NO DATA>"))
##    ow.onProbesInput(orange.ExampleTable(r"c:\Documents and Settings\peterjuv\My Documents\MF other\Patologija\2005-11-21\gasper controlRatios genes 2.tab"))

    ow.outNonNormLogRatio = True
    ow.settingsOutputChange()
    ow.normType = 0
    ow.settingsNormalizationChange()
    # OWDataTable
    dt = OWDataTable.OWDataTable(signalManager = signalManager)
    signalManager.addWidget(ow)
    signalManager.addWidget(dt)
    signalManager.setFreeze(1)
    signalManager.addLink(ow, dt, 'Probes', 'Examples', 1)
    signalManager.setFreeze(0)
    orange.saveTabDelimited(r"c:\Documents and Settings\peterjuv\My Documents\Orange\OWNormalize\test comp 2\output.tab", dt.data.values()[0])
    dt.show()
    # exec
    a.exec_loop()
    #save settings 
    ow.saveSettings()
