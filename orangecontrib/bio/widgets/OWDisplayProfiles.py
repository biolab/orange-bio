"""
<name>Data Profiles</name>
<description>Visualization of data profiles (e.g., time series).</description>
<contact>Tomaz Curk</contact>
<icon>icons/ExpressionProfiles.svg</icon>
<priority>1030</priority>
"""

from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWGraph import *
from Orange.OrangeWidgets.OWToolbars import ZoomSelectToolbar
from Orange.OrangeWidgets.OWTools import *
from Orange.OrangeWidgets.OWWidget import *
import statc

NAME = "Data Profiles"
DESCRIPTION = "Visualization of data profiles (e.g., time series)."
ICON = "icons/ExpressionProfiles.svg"
PRIORITY = 1030

INPUTS = [("Examples", Orange.data.Table, "data", Multiple + Default)]
OUTPUTS = [("Examples", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWDisplayProfiles.OWDisplayProfiles"]

## format of data:
##     largestAdjacentValue
##     yq3
##     yavg
##     yqM
##     yq1
##     smallestAdjacentValue

class boxPlotQwtPlotCurve(QwtPlotCurve):
    def __init__(self, text = None, connectPoints = 1, tickXw = 1.0/5.0):
        QwtPlotCurve.__init__(self, QwtText(text))
        self.connectPoints = connectPoints
        self.tickXw = tickXw
        self.boxPlotPenWidth = 2

    def draw(self, p, xMap, yMap, f=0, t=-1):
        # save ex settings
##        pen = p.pen()
##        brush = p.brush()
        p.save()
        if type(f)==QRect:
            f = 0 
        
        p.setBackgroundMode(Qt.OpaqueMode)
        p.setPen(self.pen())
        if self.style() == QwtPlotCurve.UserCurve:
            back = p.backgroundMode()
            p.setBackgroundMode(Qt.OpaqueMode)
            if t < 0: t = self.dataSize() - 1

            if f % 6 <> 0: f -= f % 6
            if t % 6 <> 0:  t += 6 - (t % 6)
            ## first connect averages
            first = 1
            if self.connectPoints: 
                for i in range(f, t, 6):
#                    py = yqM = yMap.transform(self.y(i + 3))
                    py = yavg = yMap.transform(self.y(i + 2))
                    px = xMap.transform(self.x(i))
                    if first:
                        first = 0
                    else:
                        p.drawLine(ppx, ppy, px, py)
                    ppx = px
                    ppy = py

            ## then draw boxes
            if self.boxPlotPenWidth > 0:
                np = QPen(self.pen())
                np.setWidth(self.boxPlotPenWidth)
                p.setPen(np)
                for i in range(f, t, 6):
                    largestAdjVal = yMap.transform(self.y(i))
                    yq3 = yMap.transform(self.y(i + 1))
                    yavg = yMap.transform(self.y(i + 2))
                    yqM = yMap.transform(self.y(i + 3))
                    yq1 = yMap.transform(self.y(i + 4))
                    smallestAdjVal = yMap.transform(self.y(i + 5))

                    px = xMap.transform(self.x(i))
                    wxl = xMap.transform(self.x(i) - self.tickXw/2.0)
                    wxr = xMap.transform(self.x(i) + self.tickXw/2.0)

                    p.drawLine(wxl, largestAdjVal, wxr,   largestAdjVal) ## - upper whisker
                    p.drawLine(px, largestAdjVal, px, yq3)               ## | connection between upper whisker and q3
                    p.drawRect(wxl, yq3, wxr - wxl, yq1 - yq3)           ## box from q3 to q1
                    p.drawLine(wxl, yqM, wxr, yqM)                       ## median line
                    p.drawLine(px, yq1, px, smallestAdjVal)              ## | connection between q1 and lower whisker
                    p.drawLine(wxl, smallestAdjVal, wxr, smallestAdjVal) ## _ lower whisker

                    ## average line (circle)
                    p.drawEllipse(px - 3, yavg - 3, 6, 6)

            p.setBackgroundMode(back)
        else:
            QwtPlotCurve.draw(self, p, xMap, yMap, f, t)

        # restore ex settings
        p.restore()
##        p.setPen(pen)
##        p.setBrush(brush)

    def setBoxPlotPenWidth(self, v):
        self.boxPlotPenWidth = v
        self.itemChanged()

class profilesGraph(OWGraph):
    def __init__(self, master, parent = None, name = None, title = ""):
        OWGraph.__init__(self, parent, name)
        self.master = master
        self.setYRlabels(None)
        self.enableGridXB(0)
        self.enableGridYL(0)
        self.setAxisMaxMajor(QwtPlot.xBottom, 10)
        self.setAxisMaxMinor(QwtPlot.xBottom, 0)
        self.setAxisMaxMajor(QwtPlot.yLeft, 10)
        self.setAxisMaxMinor(QwtPlot.yLeft, 5)
        self.setShowMainTitle(1)
        self.setMainTitle(title)
        self.setAxisAutoScale(QwtPlot.xBottom)
        self.setAxisAutoScale(QwtPlot.xTop)
        self.setAxisAutoScale(QwtPlot.yLeft)
        self.setAxisAutoScale(QwtPlot.yRight)

        self.showAverageProfile = 1
        self.showSingleProfiles = 0
        self.curveWidth = 1
        self.renderAntialiased = True
##        self.groups = [('grp1', ['0', '2', '4']), ('grp2', ['4', '6', '8', '10', '12', '14']), ('grp3', ['16', '18'])]

        self.selectedCurves = []
        self.highlightedCurve = None
        self.removeCurves()
##        self.connect(self, SIGNAL("plotMouseMoved(const QMouseEvent &)"), self.onMouseMoved)

    def removeCurves(self):
##        OWGraph.removeCurves(self)
        self.clear()
        self.classColor = None
        self.classBrighterColor = None
        self.profileCurveKeys = []
        self.averageProfileCurveKeys = []
        self.showClasses = []

    def setData(self, data, classColor, classBrighterColor, ShowAverageProfile, ShowSingleProfiles, progressBar = None):
        self.removeCurves()
        self.classColor = classColor
        self.classBrighterColor = classBrighterColor
        self.showAverageProfile = ShowAverageProfile
        self.showSingleProfiles = ShowSingleProfiles

        self.groups = [('grp', data.domain.attributes)]
        ## remove any non continuous attributes from list
        ## at the same time convert any attr. string name into orange var type
        filteredGroups = []
        for (grpname, grpattrs) in self.groups:
            filteredGrpAttrs = []
            for a in grpattrs:
                var = data.domain[a]
                if (var.varType == orange.VarTypes.Continuous):
                    filteredGrpAttrs.append(var)
                else:
                    print "warning, skipping attribute:", a
            if len(filteredGrpAttrs) > 0:
                filteredGroups.append( (grpname, filteredGrpAttrs) )
        self.groups = filteredGroups
        
        ## go group by group
        avgCurveData = []
        boxPlotCurveData = []
        ccn = 0
        if data.domain.classVar and data.domain.classVar.varType <> orange.VarTypes.Discrete:
            print "error, class variable not discrete:", data.domain.classVar
            return
        classes = data.domain.classVar.values if data.domain.classVar else ["(No class)"]
        allc = len(data.domain.classVar.values) if data.domain.classVar else 1
##        print data, len(data)
        for i, c in enumerate(classes):
            if progressBar <> None: progressBar(i*100.0/len(classes))
            classSymb = QwtSymbol(QwtSymbol.Ellipse, QBrush(self.classBrighterColor[i]), QPen(self.classBrighterColor[i]), QSize(7,7)) ##self.black
            self.showClasses.append(0)

            self.profileCurveKeys.append([])
            self.averageProfileCurveKeys.append([])
            allg = len(self.groups)
            gcn = 0
            grpcnx = 0
            for (grpname, grpattrs) in self.groups:
                oneClassData = data.select({data.domain.classVar.name:c}) if data.domain.classVar else data
                oneGrpData = oneClassData.select(orange.Domain(grpattrs, oneClassData.domain))

                ## single profiles
                nativeData = oneGrpData.native(2)
                yVals = [[] for cn in range(len(grpattrs))]
                alle = len(nativeData)
                ecn = 0
                for e, example in zip(nativeData, oneClassData):
                    ecn += 1
                    progress = 100.0*(ccn + (float(gcn)/allg) * (float(ecn)/alle))
                    progress = int(round(progress / allc))
                    if progressBar <> None: progressBar(progress)
                    y = []
                    x = []
                    xcn = grpcnx
                    vcn = 0
                    en = e.native(1)
                    for v in en:
                        if not v.isSpecial():
                            yVal = v.native()
                            yVals[vcn].append( yVal )
                            y.append( yVal )
                            x.append( xcn )
                        xcn += 1
                        vcn += 1
                    curve = self.addCurve('', style = QwtPlotCurve.Lines)
                    curve.example = example
                    if self.showAverageProfile:
                        curve.setPen(QPen(self.classColor[ccn], 1))
                    else:
                        curve.setPen(QPen(self.classBrighterColor[ccn], 1))
                    curve.setData(x, y)
                    curve.setSymbol(classSymb)
                    self.profileCurveKeys[-1].append(curve)

                ## average profile and box plot
                BPx = []
                BPy = []
                xcn = grpcnx
                vcn = 0
                dist = orange.DomainDistributions(oneGrpData)
                for a in dist:
                    if a and len(a) > 0:
                        ## box plot data
                        yavg = a.average()
                        yq1 = a.percentile(25)
                        yqM = a.percentile(50)
                        yq3 = a.percentile(75)

                        iqr = yq3 - yq1
                        yLowerCutOff = yq1 - 1.5 * iqr
                        yUpperCutOff = yq3 + 1.5 * iqr
                        
                        yVals[vcn].sort() 
                        ## find the smallest value above the lower inner fence
                        smallestAdjacentValue = None
                        for v in yVals[vcn]:
                            if v >= yLowerCutOff:
                                smallestAdjacentValue = v
                                break

                        yVals[vcn].reverse()
                        ## find the largest value below the upper inner fence
                        largestAdjacentValue = None
                        for v in yVals[vcn]:
                            if v <= yUpperCutOff:
                                largestAdjacentValue = v
                                break
                        BPy.append( largestAdjacentValue )
                        BPy.append( yq3 )
                        BPy.append( yavg )
                        BPy.append( yqM )
                        BPy.append( yq1 )
                        BPy.append( smallestAdjacentValue )
                        BPx.append( xcn )
                        BPx.append( xcn )
                        BPx.append( xcn )
                        BPx.append( xcn )
                        BPx.append( xcn )
                        BPx.append( xcn )

                    xcn += 1
                    vcn += 1

                boxPlotCurveData.append( (BPx, BPy, ccn) )
                grpcnx += len(grpattrs)
                gcn +=1
            ccn += 1

        for (x, y, tmpCcn) in boxPlotCurveData:
            classSymb = QwtSymbol(QwtSymbol.Cross, QBrush(self.classBrighterColor[tmpCcn]), QPen(self.classBrighterColor[tmpCcn]), QSize(8,8))
            curve = boxPlotQwtPlotCurve('', connectPoints = 1, tickXw = 1.0/5.0)
            curve.attach(self)
            curve.setPen(QPen(self.classBrighterColor[tmpCcn], 3))
            curve.setStyle(QwtPlotCurve.UserCurve)
            curve.setSymbol(classSymb)
            curve.setData(x, y)
            self.averageProfileCurveKeys[tmpCcn].append(curve)

        ## generate labels for attributes
        labels = []
        for (grpname, grpattrs) in self.groups:
            for a in grpattrs:
                labels.append( a.name)

        if None: 
            self.setXlabels(labels)
        self.updateCurveDisplay()

    def updateCurveDisplay(self):
        profilesCount = 0
        for cNum in range(len(self.showClasses or [1])):
            showCNum = (self.showClasses[cNum] <> 0) if len(self.showClasses) > 1 else True
##            print cNum, showCNum
            ## single profiles
            bSingle = showCNum and self.showSingleProfiles
            bAve = showCNum and self.showAverageProfile ## 1 = show average profiles for now
            for ckey in self.profileCurveKeys[cNum]:
                curve =  ckey #self.curve(ckey)
                if curve <> None:
                    curve.setVisible(bSingle)
                    profilesCount += 1 if bSingle else 0
##                    qp = self.curvePen(ckey)
                    qp = curve.pen()
                    if not(bAve):
                        qp.setColor(self.classBrighterColor[cNum])
                    else:
                        qp.setColor(self.classColor[cNum])
##                    curve.setPen(qp)

            ## average profiles
            for ckey in self.averageProfileCurveKeys[cNum]:
                curve =  ckey #self.curve(ckey)
                if curve <> None: curve.setVisible(bAve)

##        self.updateLayout()
        self.master.infoLabel.setText("Showing %i profiles" % profilesCount)
        self.setCurveRenderHints()
        self.replot()

    def setShowClasses(self, list):
        self.showClasses = list
        self.updateCurveDisplay()

    def setShowAverageProfile(self, v):
        self.showAverageProfile = v
        self.updateCurveDisplay()

    def setShowSingleProfiles(self, v):
        self.showSingleProfiles = v
        self.updateCurveDisplay()

    def setPointWidth(self, v):
        for cNum in range(len(self.showClasses)):
            for curve in self.profileCurveKeys[cNum]:
                symb = curve.symbol()
                symb.setSize(v, v)
##                curve.setSymbol(symb)
        self.replot()

    def setCurveWidth(self, v):
        fix = dict.fromkeys(self.selectedCurves, 2)
        fix[self.highlightedCurve] = 3
        for cNum in range(len(self.showClasses)):
            for curve in self.profileCurveKeys[cNum]:
                qp = curve.pen()
                qp.setWidth(v + fix.get(curve, 0))
##                curve.setPen(qp)
        self.curveWidth = v
        self.replot()

    def setAverageCurveWidth(self, v):
        for cNum in range(len(self.showClasses)):
            for curve in self.averageProfileCurveKeys[cNum]:
                qp = curve.pen()
                qp.setWidth(v)
##                curve.setPen(qp)
        self.replot()

    def setBoxPlotWidth(self, v):
        for cNum in range(len(self.showClasses)):
            for curve in self.averageProfileCurveKeys[cNum]:
##                c = self.curve(ckey)
                curve.setBoxPlotPenWidth(v)
        self.replot()

    def sizeHint(self):
        return QSize(170, 170)

    def closestCurve(self, point):
        pointDistances = [(curve,) + curve.closestPoint(point) for curve in self.itemList() if isinstance(curve, QwtPlotCurve) and curve.isVisible()]
        return min(pointDistances, key=lambda t:t[-1]) if pointDistances else (None, -1, sys.maxint)

    def curveUnderMouse(self, pos):
        pos = self.canvas().mapFrom(self, pos)
        curve, i, d = self.closestCurve(pos)
        if d <= 5 and hasattr(curve, "example"):
            return curve
        else:
            return None
    
    def mouseMoveEvent(self, event):
        curve = self.curveUnderMouse(event.pos())
        if self.highlightedCurve != curve:
            if self.highlightedCurve:
                self.highlightedCurve.setPen(QPen(self.highlightedCurve.pen().color(), self.curveWidth + (2 if self.highlightedCurve in self.selectedCurves else 0)))
            if curve:
                curve.setPen(QPen(curve.pen().color(), self.curveWidth + 3))
            self.highlightedCurve = curve
            self.replot()
            
            if curve:
                QToolTip.showText(event.globalPos(), "")
                QToolTip.showText(event.globalPos(), str(curve.example[self.master.profileLabel]))
                                         
        return OWGraph.mouseMoveEvent(self, event)
    
    def mousePressEvent(self, event):
        curve = self.curveUnderMouse(event.pos())
        if curve and event.button() == Qt.LeftButton and self.state == SELECT:
            if self.master.ctrlPressed:
                if curve in self.selectedCurves:
                    self.setSelectedCurves([c for c in self.selectedCurves if c is not curve])
                else:
                    self.setSelectedCurves(self.selectedCurves + [curve])
            else:
                self.setSelectedCurves([curve])
            self.replot()
                
        return OWGraph.mousePressEvent(self, event)
    
    def setCurveRenderHints(self):
        for item in self.itemList():
            if isinstance(item, QwtPlotCurve):
                item.setRenderHint(QwtPlotCurve.RenderAntialiased, self.renderAntialiased)
        self.replot()
        
    def setSelectedCurves(self, curves=[]):
        for c in self.selectedCurves:
            c.setPen(QPen(c.pen().color(), self.curveWidth))
        for c in curves:
            c.setPen(QPen(c.pen().color(), self.curveWidth + 2))
        self.selectedCurves = curves
        self.master.commitIf()
        
    def removeLastSelection(self):
        self.setSelectedCurves(self.selectedCurves[:-1])
        self.replot()
        
    def removeAllSelections(self):
        self.setSelectedCurves([])
        self.replot()
        
    def sendData(self):
        self.master.commitIf()


class OWDisplayProfiles(OWWidget):
    settingsList = ["SelectedClasses", "PointWidth", "CurveWidth", "AverageCurveWidth", "BoxPlotWidth", "ShowAverageProfile", "ShowSingleProfiles", "CutEnabled", "CutLow", "CutHigh", "autoSendSelected"]
    contextHandlers = {"": DomainContextHandler("", [ContextField("SelectedClasses")], loadImperfect=False)}
    def __init__(self, parent=None, signalManager = None):
#        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, signalManager, 'Expression Profiles', 1)

        #set default settings
        self.ShowAverageProfile = 1
        self.ShowSingleProfiles = 0
        self.PointWidth = 2
        self.CurveWidth = 1
        self.AverageCurveWidth = 4
        self.BoxPlotWidth = 2
        self.SelectedClasses = []
        self.Classes = []
        self.autoSendSelected = 0
        self.selectionChangedFlag = False

        self.CutLow = 0; self.CutHigh = 0; self.CutEnabled = 0
        
        self.profileLabel = None

        #load settings
        self.loadSettings()

        # GUI
        self.graph = profilesGraph(self, self.mainArea, "")
        self.mainArea.layout().addWidget(self.graph)
        self.graph.hide()
        self.connect(self.graphButton, SIGNAL("clicked()"), self.graph.saveToFile)

        # GUI definition
        self.tabs = OWGUI.tabWidget(self.space)

        # GRAPH TAB
        GraphTab = OWGUI.createTabPage(self.tabs, "Graph")

        ## display options
        self.infoLabel = OWGUI.widgetLabel(OWGUI.widgetBox(GraphTab, "Info"), "No data on input.")
        displayOptBox = OWGUI.widgetBox(GraphTab, "Display") #QVButtonGroup("Display", GraphTab)
        displayOptButtons = ['Majority Class', 'Majority Class Probability', 'Target Class Probability', 'Number of Instances']
        OWGUI.checkBox(displayOptBox, self, 'ShowSingleProfiles', 'Expression Profiles', tooltip='', callback=self.updateShowSingleProfiles)
        OWGUI.checkBox(displayOptBox, self, 'ShowAverageProfile', 'Box Plot', tooltip='', callback=self.updateShowAverageProfile)

        ## class selection (classQLB)
        self.classQVGB = OWGUI.widgetBox(GraphTab, "Classes") #QVGroupBox(GraphTab)
##        self.classQVGB.setTitle("Classes")
        self.classQLB = OWGUI.listBox(self.classQVGB, self, "SelectedClasses", "Classes", selectionMode=QListWidget.MultiSelection, callback=self.classSelectionChange) #QListBox(self.classQVGB)
##        self.classQLB.setSelectionMode(QListBox.Multi)
        self.unselectAllClassedQLB = OWGUI.button(self.classQVGB, self, "Unselect all", callback=self.SUAclassQLB) #QPushButton("(Un)Select All", self.classQVGB)
##        self.connect(self.unselectAllClassedQLB, SIGNAL("clicked()"), self.SUAclassQLB)
##        self.connect(self.classQLB, SIGNAL("selectionChanged()"), self.classSelectionChange)

        ## show single/average profile
##        self.showAverageQLB = QPushButton("Box Plot", self.classQVGB)
##        self.showAverageQLB.setToggleButton(1)
##        self.showAverageQLB.setOn(self.ShowAverageProfile)
##        self.showSingleQLB = QPushButton("Single Profiles", self.classQVGB)
##        self.showSingleQLB.setToggleButton(1)
##        self.showSingleQLB.setOn(self.ShowSingleProfiles)
##        self.connect(self.showAverageQLB, SIGNAL("toggled(bool)"), self.updateShowAverageProfile)
##        self.connect(self.showSingleQLB, SIGNAL("toggled(bool)"), self.updateShowSingleProfiles)

##        self.tabs.insertTab(GraphTab, "Graph")

        # SETTINGS TAB
        SettingsTab = OWGUI.createTabPage(self.tabs, "Settings") #QVGroupBox(self)

        self.profileLabelComboBox = OWGUI.comboBox(SettingsTab, self, 'profileLabel', 'Profile Labels', sendSelectedValue=True, valueType=str)
        OWGUI.hSlider(SettingsTab, self, 'PointWidth', box='Point Width', minValue=0, maxValue=9, step=1, callback=self.updatePointWidth, ticks=1)
        OWGUI.hSlider(SettingsTab, self, 'CurveWidth', box='Profile Width', minValue=1, maxValue=9, step=1, callback=self.updateCurveWidth, ticks=1)
        OWGUI.hSlider(SettingsTab, self, 'AverageCurveWidth', box='Average Profile Width', minValue=1, maxValue=9, step=1, callback=self.updateAverageCurveWidth, ticks=1)
        OWGUI.hSlider(SettingsTab, self, 'BoxPlotWidth', box='Box Plot Width', minValue=0, maxValue=9, step=1, callback=self.updateBoxPlotWidth, ticks=1)
        OWGUI.checkBox(SettingsTab, self, 'graph.renderAntialiased', "Render antialiased", callback=self.graph.setCurveRenderHints)


    ## graph y scale
        box = OWGUI.widgetBox(SettingsTab, "Threshold/ Values") #QVButtonGroup("Threshol/ Values", SettingsTab)
        OWGUI.checkBox(box, self, 'CutEnabled', "Enabled", callback=self.setCutEnabled)
        self.sliderCutLow = OWGUI.qwtHSlider(box, self, 'CutLow', label='Low:', labelWidth=33, minValue=-20, maxValue=0, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.updateYaxis)
        self.sliderCutHigh = OWGUI.qwtHSlider(box, self, 'CutHigh', label='High:', labelWidth=33, minValue=0, maxValue=20, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.updateYaxis)
        if not self.CutEnabled:
            self.sliderCutLow.box.setDisabled(1)
            self.sliderCutHigh.box.setDisabled(1)

##        self.tabs.insertTab(SettingsTab, "Settings")

        self.toolbarSelection = ZoomSelectToolbar(self, self.controlArea, self.graph, self.autoSendSelected, buttons=(ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan, ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconSpace,
                                                                       ZoomSelectToolbar.IconRemoveLast, ZoomSelectToolbar.IconRemoveAll, ZoomSelectToolbar.IconSendSelection))
        cb = OWGUI.checkBox(self.controlArea, self, "autoSendSelected", "Auto send on selection change")
        OWGUI.setStopper(self, self.toolbarSelection.buttonSendSelections, cb, "selectionChangedFlag", self.commit)
        
        # inputs
        # data and graph temp variables
        
        self.inputs = [("Examples", ExampleTable, self.data, Default + Multiple)]
        self.outputs = [("Examples", ExampleTable, Default)]

        # temp variables
        self.MAdata = []
        self.classColor = None
        self.classBrighterColor = None
        self.numberOfClasses  = 0
        self.classValues = []
        self.ctrlPressed = False
        self.attrIcons = self.createAttributeIconDict()

        self.graph.canvas().setMouseTracking(1)

#        self.zoomStack = []
##        self.connect(self.graph,
##                     SIGNAL('plotMousePressed(const QMouseEvent&)'),
##                     self.onMousePressed)
##        self.connect(self.graph,
##                     SIGNAL('plotMouseReleased(const QMouseEvent&)'),
##                     self.onMouseReleased)
        self.resize(800,600)

    def updateYaxis(self):
        if (not self.CutEnabled):
            self.graph.setAxisAutoScale(QwtPlot.yLeft)
            self.graph.replot()
            return
                                                                                                                                     
        self.graph.setAxisScale(QwtPlot.yLeft, self.CutLow, self.CutHigh)
        self.graph.replot()

    def setCutEnabled(self):
        self.sliderCutLow.box.setDisabled(not self.CutEnabled)
        self.sliderCutHigh.box.setDisabled(not self.CutEnabled)
        self.updateYaxis()

#    def onMousePressed(self, e):
#        if Qt.LeftButton == e.button():
#            # Python semantics: self.pos = e.pos() does not work; force a copy
#            self.xpos = e.pos().x()
#            self.ypos = e.pos().y()
#            self.graph.enableOutline(1)
#            self.graph.setOutlinePen(QPen(Qt.black))
#            self.graph.setOutlineStyle(Qwt.Rect)
#            self.zooming = 1
#            if self.zoomStack == []:
#                self.zoomState = (
#                    self.graph.axisScale(QwtPlot.xBottom).lBound(),
#                    self.graph.axisScale(QwtPlot.xBottom).hBound(),
#                    self.graph.axisScale(QwtPlot.yLeft).lBound(),
#                    self.graph.axisScale(QwtPlot.yLeft).hBound(),
#                    )
#        elif Qt.RightButton == e.button():
#            self.zooming = 0
#        # fake a mouse move to show the cursor position
#
#    # onMousePressed()
#
#    def onMouseReleased(self, e):
#        if Qt.LeftButton == e.button():
#            xmin = min(self.xpos, e.pos().x())
#            xmax = max(self.xpos, e.pos().x())
#            ymin = min(self.ypos, e.pos().y())
#            ymax = max(self.ypos, e.pos().y())
#            self.graph.setOutlineStyle(Qwt.Cross)
#            xmin = self.graph.invTransform(QwtPlot.xBottom, xmin)
#            xmax = self.graph.invTransform(QwtPlot.xBottom, xmax)
#            ymin = self.graph.invTransform(QwtPlot.yLeft, ymin)
#            ymax = self.graph.invTransform(QwtPlot.yLeft, ymax)
#            if xmin == xmax or ymin == ymax:
#                return
#            self.zoomStack.append(self.zoomState)
#            self.zoomState = (xmin, xmax, ymin, ymax)
#            self.graph.enableOutline(0)
#        elif Qt.RightButton == e.button():
#            if len(self.zoomStack):
#                xmin, xmax, ymin, ymax = self.zoomStack.pop()
#            else:
#                self.graph.setAxisAutoScale(QwtPlot.xBottom)
#                self.graph.setAxisAutoScale(QwtPlot.yLeft)
#                self.graph.replot()
#                return
#
#        self.graph.setAxisScale(QwtPlot.xBottom, xmin, xmax)
#        self.graph.setAxisScale(QwtPlot.yLeft, ymin, ymax)
#        self.graph.replot()
#
#    def saveToFile(self):
#        qfileName = QFileDialog.getSaveFileName("graph.png","Portable Network Graphics (.PNG)\nWindows Bitmap (.BMP)\nGraphics Interchange Format (.GIF)", None, "Save to..")
#        fileName = str(qfileName)
#        if fileName == "": return
#        (fil,ext) = os.path.splitext(fileName)
#        ext = ext.replace(".","")
#        ext = ext.upper()
#        cl = 0
#        for g in self.graphs:
#            if g.isVisible():
#                clfname = fil + "_" + str(cl) + "." + ext
#                g.saveToFileDirect(clfname, ext)
#            cl += 1

    def updateShowAverageProfile(self):
        self.graph.setShowAverageProfile(self.ShowAverageProfile)

    def updateShowSingleProfiles(self):
        self.graph.setShowSingleProfiles(self.ShowSingleProfiles)

    def updatePointWidth(self):
        self.graph.setPointWidth(self.PointWidth)

    def updateCurveWidth(self):
        self.graph.setCurveWidth(self.CurveWidth)

    def updateAverageCurveWidth(self):
        self.graph.setAverageCurveWidth(self.AverageCurveWidth)

    def updateBoxPlotWidth(self):
        self.graph.setBoxPlotWidth(self.BoxPlotWidth)
        
    ##
    def selectUnselectAll(self, qlb):
        self.SelectedClasses = range(len(self.Classes)) if len(self.SelectedClasses) != len(self.Classes)  else []
        
        self.classSelectionChange()

    def SUAclassQLB(self):
        self.selectUnselectAll(self.classQLB)
    ##

    ## class selection (classQLB)
    def classSelectionChange(self):
##        list = []
##        selCls = []
##        for i in range(self.classQLB.count()):
##            if self.classQLB.isSelected(i):
##                list.append( 1 )
##                selCls.append(self.classValues[i])
##            else:
##                list.append( 0 )
        list = [1 if i in self.SelectedClasses else 0 for i in range(self.classQLB.count())]
        self.unselectAllClassedQLB.setText("Select all" if len(self.SelectedClasses) != len(self.Classes) else "Unselect all")
        selCls = [self.classValues[i] for i in self.SelectedClasses]
        self.graph.setShowClasses(list)
#        if selCls == []:
#            self.send("Examples", None)
#        else:
#            if len(self.MAdata) > 1:
#                newet = orange.ExampleTable(self.MAdata[0].domain)
#                for idx, i in enumerate(list):
#                    if i: newet.extend(self.MAdata[idx])
#            elif self.MAdata[0].domain.classVar:
#                newet = self.MAdata[0].filter({self.MAdata[0].domain.classVar.name:selCls})
#            else:
#                newet = self.MAdata[0]
#            self.send("Examples", newet)
    ##

    def calcGraph(self):
        ## compute from self.MAdata
        if len(self.MAdata) == 1:
            d = self.MAdata[0]
        elif len(self.MAdata) > 0:
            combCvals = [str(nd.name) for nd in self.MAdata]
            combClass = orange.EnumVariable('file', values=combCvals)
            combDomain = orange.Domain(self.MAdata[0].domain.attributes + [combClass])
            d = orange.ExampleTable(combDomain)
            
            for tmpd in self.MAdata:
                newtmpd = tmpd.select(combDomain)
                for te in newtmpd:
                    te.setclass(tmpd.name)
                d.extend(newtmpd)
        else:
            return

        self.progressBarInit()
        self.graph.setData(d, self.classColor, self.classBrighterColor, self.ShowAverageProfile, self.ShowSingleProfiles, self.progressBarSet)
        self.graph.setPointWidth(self.PointWidth)
        self.graph.setCurveWidth(self.CurveWidth)
        self.graph.setAverageCurveWidth(self.AverageCurveWidth)
        self.graph.setBoxPlotWidth(self.BoxPlotWidth)

        self.graph.setAxisAutoScale(QwtPlot.xBottom)
##        self.graph.setAxisAutoScale(QwtPlot.yLeft)
        self.updateYaxis()
        self.progressBarFinished()

    def newdata(self):
        self.classQLB.clear()
        self.profileLabelComboBox.clear()
##        if len(self.MAdata) > 1 or (len(self.MAdata) == 1 and self.MAdata[0].domain.classVar.varType == orange.VarTypes.Discrete):
        if len(self.MAdata) >= 1:
            ## classQLB
            if len(self.MAdata) == 1 and self.MAdata[0].domain.classVar and self.MAdata[0].domain.classVar.varType == orange.VarTypes.Discrete:
                self.numberOfClasses = len(self.MAdata[0].domain.classVar.values)
            else:
                self.numberOfClasses = len(self.MAdata)
            self.classColor = ColorPaletteGenerator(self.numberOfClasses)
            self.classBrighterColor = [self.classColor[c, 160] for c in range(self.numberOfClasses)]

            self.calcGraph()
            ## update graphics
            ## classQLB
            self.classQVGB.show()
            if len(self.MAdata) == 1 and self.MAdata[0].domain.classVar and self.MAdata[0].domain.classVar.varType == orange.VarTypes.Discrete:
                self.classValues = self.MAdata[0].domain.classVar.values.native()
            else:
                self.classValues = [str(nd.name) for nd in self.MAdata]

##            selection = QItemSelection()
##            for cn in range(len(self.classValues)):
##                item = QListWidgetItem(QIcon(ColorPixmap(self.classBrighterColor[cn])), self.classValues[cn])
##                self.classQLB.addItem(item)
##                selection.select(self.classQLB.indexFromItem(item), self.classQLB.indexFromItem(item))
##            self.classQLB.selectionModel().select(selection, QItemSelectionModel.Select)
            self.Classes = [(QIcon(ColorPixmap(self.classBrighterColor[cn])), self.classValues[cn])  for cn in range(len(self.classValues))]
            self.SelectedClasses = range(len(self.Classes))
            
##            self.classQLB.selectAll()  ##or: if numberOfClasses > 0: self.classQLB.setSelected(0, 1)
##            self.update()

##            if len(self.MAdata) == 1 and self.MAdata[0].noclass:
##                pass
            self.classQVGB.setDisabled(len(self.MAdata) == 1 and self.MAdata[0].noclass)
            attrs = self.MAdata[0].domain.variables + self.MAdata[0].domain.getmetas().values()
            for attr in attrs:
                self.profileLabelComboBox.addItem(self.attrIcons[attr.varType], attr.name)
            stringAttrs = [attr for attr in attrs if attr.varType == orange.VarTypes.String]
            discAttrs = [attr for attr in attrs if attr.varType == orange.VarTypes.Discrete]
            attrs = stringAttrs[-1:] + discAttrs[-1:] + list(attrs[:1])
            self.profileLabel = attrs[0].name
        else:
            self.classColor = None
            self.classBrighterColor = None
            self.classValues = []
##        self.classQVGB.update()
##        self.classQVGB.layout().activate()
        self.graph.show()
##        self.layout.activate() # this is needed to scale the widget correctly

    def data(self, MAdata, id=None):
        ## if there is no class attribute, create a dummy one
        self.closeContext()
        if MAdata and MAdata.domain.classVar == None:
##            noClass = orange.EnumVariable('file', values=['n', 'y'])
##            newDomain = orange.Domain(MAdata.domain.attributes + [noClass])
##            mname = MAdata.name ## remember name 
##            MAdata = MAdata.select(newDomain) ## because select forgets it
##            MAdata.name = mname
##            for e in MAdata: e.setclass('n')
            MAdata.setattr("noclass", 1) ## remember that there is no class to display
        elif MAdata and MAdata.domain.classVar.varType <> orange.VarTypes.Discrete:
            print "error, ignoring table, because its class variable not discrete:", MAdata.domain.classVar
##            MAdata = None
            MAdata.setattr("noclass", 1) ## remember that there is no class to display
        elif MAdata:
            MAdata.setattr("noclass", 0) ## there are classes by default

        ## handling of more than one data set
        # check if the same domain
        if MAdata and len(self.MAdata) and str(MAdata.domain.attributes) <> str(self.MAdata[0].domain.attributes):
            print "domains:", MAdata.domain.attributes
            print "and,   :", self.MAdata[0].domain.attributes
            print "are not same"
##            MAdata = None

        ids = [d.id for d in self.MAdata]
        if not MAdata:
            if id in ids:
                del self.MAdata[ids.index(id)]
        else:
            MAdata.setattr("id", id)
            if id in ids:
                MAdata.setattr("id", id)
                indx = ids.index(id)
                self.MAdata[indx] = MAdata
            else:
                self.MAdata.append(MAdata)

        if len(self.MAdata) == 0:
            print "hiding graph"
            self.graph.hide()
            self.classQVGB.hide()
            self.infoLabel.setText("No data on input")
            return

        self.newdata()
        self.openContext("", MAdata)
        self.classSelectionChange()

    def sendReport(self):
        self.startReport("%s" % self.windowTitle())
        self.reportSettings("Settings", ([("Selected classes" , ",".join(self.Classes[i][1] for i in self.SelectedClasses))] if len(self.Classes) > 1 else []) +\
                                        [("Show box plot", self.ShowAverageProfile),
                                         ("Show profiles", self.ShowSingleProfiles)])
        self.reportRaw("<p>%s</p>" % self.infoLabel.text())
        self.reportImage(lambda *x: OWChooseImageSizeDlg(self.graph).saveImage(*x))
        
    def keyPressEvent(self, event):
        self.ctrlPressed = event.key() == Qt.Key_Control
        return OWWidget.keyPressEvent(self, event)
    
    def keyReleaseEvent(self, event):
        self.ctrlPressed = self.ctrlPressed and not event.key() == Qt.Key_Control
        return OWWidget.keyReleaseEvent(self, event)
    
    def commitIf(self):
        if self.autoSendSelected:
            self.commit()
        else:
            self.selectionChangedFlag = True
            
    def commit(self):
        data = [c.example for c in self.graph.selectedCurves if hasattr(c, "example")]
        if data:
            data = orange.ExampleTable(data)
        else:
            data = None
            
        self.send("Examples", data)
        self.selectionChangedFlag = False
    
# following is not needed, data handles these cases
##    def cdata(self, MAcdata):
##        if not MAcdata:
##            self.graph.hide()
##            return
##        self.MAdata = MAcdata
##        self.MAnoclass = 0
##        self.newdata()

if __name__ == "__main__":
    a = QApplication(sys.argv)
    owdm = OWDisplayProfiles()
##    a.setMainWidget(owdm)
##    d = orange.ExampleTable('e:\\profiles')
##    d = orange.ExampleTable('e:\\profiles-classes')
    d = orange.ExampleTable('../../../doc/datasets/brown-selected')
    print len(d)
    owdm.data(d)
    owdm.show()
    a.exec_()
    owdm.saveSettings()
