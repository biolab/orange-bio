"""
<name>Expression Profiles</name>
<description>Displays gene expression profiles.</description>
<category>Genomics</category>
<icon>icons/ExpressionProfiles.png</icon>
<priority>10</priority>
"""

from OWTools import *
from OWWidget import *
from OWGraph import *
import OWGUI

from OWGraph import ColorPaletteHSV

import statc

## format of data:
##     largestAdjacentValue
##     yq3
##     yavg
##     yqM
##     yq1
##     smallestAdjacentValue

class boxPlotQwtPlotCurve(QwtPlotCurve):
    def __init__(self, parent = None, text = None, connectPoints = 1, tickXw = 1.0/5.0):
        QwtPlotCurve.__init__(self, parent, text)
        self.connectPoints = connectPoints
        self.tickXw = tickXw
        self.boxPlotPenWidth = 2

    def draw(self, p, xMap, yMap, f, t):
        # save ex settings
        pen = p.pen()
        brush = p.brush()
        
        p.setBackgroundMode(Qt.OpaqueMode)
        p.setPen(self.pen())
        if self.style() == QwtCurve.UserCurve:
            back = p.backgroundMode()
            p.setBackgroundMode(Qt.OpaqueMode)
            if t < 0: t = self.dataSize() - 1

            if f % 6 <> 0: f -= f % 6
            if t % 6 <> 0:  t += 6 - (t % 6)

            ## first connect medians            
            first = 1
            if self.connectPoints: 
                for i in range(f, t, 6):
                    py = yqM = yMap.transform(self.y(i + 3))
                    px = xMap.transform(self.x(i))
                    if first:
                        first = 0
                    else:
                        p.drawLine(ppx, ppy, px, py)
                    ppx = px
                    ppy = py

            ## then draw boxes
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
        p.setPen(pen)
        p.setBrush(brush)

    def setBoxPlotPenWidth(self, v):
        self.boxPlotPenWidth = v

class profilesGraph(OWGraph):
    def __init__(self, parent = None, name = None, title = ""):
        OWGraph.__init__(self, parent, name)
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
##        self.groups = [('grp1', ['0', '2', '4']), ('grp2', ['4', '6', '8', '10', '12', '14']), ('grp3', ['16', '18'])]

        self.removeCurves()
##        self.connect(self, SIGNAL("plotMouseMoved(const QMouseEvent &)"), self.onMouseMoved)

    def removeCurves(self):
        OWGraph.removeCurves(self)
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
        if data.domain.classVar.varType <> orange.VarTypes.Discrete:
            print "error, class variable not discrete:", data.domain.classVar
            return
        allc = len(data.domain.classVar.values)
        for c in data.domain.classVar.values:
            if progressBar <> None: progressBar(int(ccn*100.0/allc))
            classSymb = QwtSymbol(QwtSymbol.Ellipse, QBrush(self.classBrighterColor[ccn]), QPen(self.classBrighterColor[ccn]), QSize(7,7)) ##self.black
            self.showClasses.append(0)

            self.profileCurveKeys.append([])
            self.averageProfileCurveKeys.append([])
            allg = len(self.groups)
            gcn = 0
            grpcnx = 0
            for (grpname, grpattrs) in self.groups:
                oneClassData = data.select({data.domain.classVar.name:c})
                oneGrpData = oneClassData.select(orange.Domain(grpattrs, oneClassData.domain))

                ## single profiles
                nativeData = oneGrpData.native(2)
                yVals = [[] for cn in range(len(grpattrs))]
                alle = len(nativeData)
                ecn = 0
                for e in nativeData:
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
                    ckey = self.insertCurve('')
                    if self.showAverageProfile:
                        self.setCurvePen(ckey, QPen(self.classColor[ccn], 1))
                    else:
                        self.setCurvePen(ckey, QPen(self.classBrighterColor[ccn], 1))
                    self.setCurveData(ckey, x, y)
                    self.setCurveSymbol(ckey, classSymb)
                    self.profileCurveKeys[-1].append(ckey)

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
            curve = boxPlotQwtPlotCurve(self, '', connectPoints = 1, tickXw = 1.0/5.0)
            ckey = self.insertCurve(curve)
            self.setCurvePen(ckey, QPen(self.classBrighterColor[tmpCcn], 3))
            self.setCurveStyle(ckey, QwtCurve.UserCurve)
            self.setCurveSymbol(ckey, classSymb)
            self.setCurveData(ckey, x, y)
            self.averageProfileCurveKeys[tmpCcn].append(ckey)

        ## generate labels for attributes
        labels = []
        for (grpname, grpattrs) in self.groups:
            for a in grpattrs:
                labels.append( a.name)

        self.setXlabels(labels)
        self.updateCurveDisplay()

    def updateCurveDisplay(self):
        for cNum in range(len(self.showClasses)):
            showCNum = (self.showClasses[cNum] <> 0)

            ## single profiles
            bSingle = showCNum and self.showSingleProfiles
            bAve = showCNum and self.showAverageProfile ## 1 = show average profiles for now
            for ckey in self.profileCurveKeys[cNum]:
                curve =  self.curve(ckey)
                if curve <> None:
                    curve.setEnabled(bSingle)
                    qp = self.curvePen(ckey)
                    if not(bAve):
                        qp.setColor(self.classBrighterColor[cNum])
                    else:
                        qp.setColor(self.classColor[cNum])
                    self.setCurvePen(ckey, qp)

            ## average profiles
            for ckey in self.averageProfileCurveKeys[cNum]:
                curve =  self.curve(ckey)
                if curve <> None: curve.setEnabled(bAve)

        self.updateLayout()
        self.update()

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
            for ckey in self.profileCurveKeys[cNum]:
                symb = self.curveSymbol(ckey)
                symb.setSize(v, v)
                self.setCurveSymbol(ckey, symb)
        self.update()

    def setCurveWidth(self, v):
        for cNum in range(len(self.showClasses)):
            for ckey in self.profileCurveKeys[cNum]:
                qp = self.curvePen(ckey)
                qp.setWidth(v)
                self.setCurvePen(ckey, qp)
        self.update()

    def setAverageCurveWidth(self, v):
        for cNum in range(len(self.showClasses)):
            for ckey in self.averageProfileCurveKeys[cNum]:
                qp = self.curvePen(ckey)
                qp.setWidth(v)
                self.setCurvePen(ckey, qp)
        self.update()

    def setBoxPlotWidth(self, v):
        for cNum in range(len(self.showClasses)):
            for ckey in self.averageProfileCurveKeys[cNum]:
                c = self.curve(ckey)
                c.setBoxPlotPenWidth(v)
                c.curveChanged()
        self.update()

    def sizeHint(self):
        return QSize(170, 170)

    def onMouseMoved(self, e):
        (key, foo1, x, y, foo2) = self.closestCurve(e.pos().x(), e.pos().y())
##        print e.pos().x(), e.pos().y(), key, foo1, x, y, foo2
##        print self.invTransform(QwtPlot.xBottom, e.pos().x()), self.invTransform(QwtPlot.yLeft, e.pos().y())


class OWDisplayProfiles(OWWidget):
    settingsList = ["PointWidth", "CurveWidth", "AverageCurveWidth", "BoxPlotWidth", "ShowAverageProfile", "ShowSingleProfiles"]
    def __init__(self,parent=None, name='Expression Profiles'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, name, """None.""", FALSE, TRUE)

        #set default settings
        self.ShowAverageProfile = 1
        self.ShowSingleProfiles = 0
        self.PointWidth = 2
        self.CurveWidth = 1
        self.AverageCurveWidth = 4
        self.BoxPlotWidth = 2

        #load settings
        self.loadSettings()

        # GUI
        self.layout=QVBoxLayout(self.mainArea)
        self.graph = profilesGraph(self.mainArea, "")
        self.layout.add(self.graph)
        self.graph.hide()
        self.connect(self.graphButton, SIGNAL("clicked()"), self.graph.saveToFile)

        # GUI definition
        self.tabs = QTabWidget(self.space, 'tabWidget')

        # GRAPH TAB
        GraphTab = QVGroupBox(self)

        ## display options
        displayOptBox = QVButtonGroup("Display", GraphTab)
        displayOptButtons = ['Majority Class', 'Majority Class Probability', 'Target Class Probability', 'Number of Instances']
        OWGUI.checkBox(displayOptBox, self, 'ShowSingleProfiles', 'Expression Profiles', tooltip='', callback=self.updateShowSingleProfiles)
        OWGUI.checkBox(displayOptBox, self, 'ShowAverageProfile', 'Box Plot', tooltip='', callback=self.updateShowAverageProfile)

        ## class selection (classQLB)
        self.classQVGB = QVGroupBox(GraphTab)
        self.classQVGB.setTitle("Classes")
        self.classQLB = QListBox(self.classQVGB)
        self.classQLB.setSelectionMode(QListBox.Multi)
        self.unselectAllClassedQLB = QPushButton("(Un)Select All", self.classQVGB)
        self.connect(self.unselectAllClassedQLB, SIGNAL("clicked()"), self.SUAclassQLB)
        self.connect(self.classQLB, SIGNAL("selectionChanged()"), self.classSelectionChange)

        ## show single/average profile
##        self.showAverageQLB = QPushButton("Box Plot", self.classQVGB)
##        self.showAverageQLB.setToggleButton(1)
##        self.showAverageQLB.setOn(self.ShowAverageProfile)
##        self.showSingleQLB = QPushButton("Single Profiles", self.classQVGB)
##        self.showSingleQLB.setToggleButton(1)
##        self.showSingleQLB.setOn(self.ShowSingleProfiles)
##        self.connect(self.showAverageQLB, SIGNAL("toggled(bool)"), self.updateShowAverageProfile)
##        self.connect(self.showSingleQLB, SIGNAL("toggled(bool)"), self.updateShowSingleProfiles)


        self.tabs.insertTab(GraphTab, "Graph")

        # SETTINGS TAB
        SettingsTab = QVGroupBox(self)

        OWGUI.hSlider(SettingsTab, self, 'PointWidth', box='Point Width', minValue=0, maxValue=9, step=1, callback=self.updatePointWidth, ticks=1)
        OWGUI.hSlider(SettingsTab, self, 'CurveWidth', box='Profile Width', minValue=1, maxValue=9, step=1, callback=self.updateCurveWidth, ticks=1)
        OWGUI.hSlider(SettingsTab, self, 'AverageCurveWidth', box='Average Profile Width', minValue=1, maxValue=9, step=1, callback=self.updateAverageCurveWidth, ticks=1)
        OWGUI.hSlider(SettingsTab, self, 'BoxPlotWidth', box='Box Plot Width', minValue=1, maxValue=9, step=1, callback=self.updateBoxPlotWidth, ticks=1)

        self.tabs.insertTab(SettingsTab, "Settings")
        
        # inputs
        # data and graph temp variables
        
        self.xxxinputs = [("Classified Examples", ExampleTableWithClass, self.data, 1), ("Examples", ExampleTable, self.data, 1)]
        self.inputs = [("Examples", ExampleTable, self.data, 1)]

        # temp variables
        self.MAdata = None
        self.MAnoclass = 1 
        self.classColor = None
        self.classBrighterColor = None
        self.numberOfClasses  = 0

        self.graph.canvas().setMouseTracking(1)

        self.zoomStack = []
        self.connect(self.graph,
                     SIGNAL('plotMousePressed(const QMouseEvent&)'),
                     self.onMousePressed)
        self.connect(self.graph,
                     SIGNAL('plotMouseReleased(const QMouseEvent&)'),
                     self.onMouseReleased)

    def onMousePressed(self, e):
        if Qt.LeftButton == e.button():
            # Python semantics: self.pos = e.pos() does not work; force a copy
            self.xpos = e.pos().x()
            self.ypos = e.pos().y()
            self.graph.enableOutline(1)
            self.graph.setOutlinePen(QPen(Qt.black))
            self.graph.setOutlineStyle(Qwt.Rect)
            self.zooming = 1
            if self.zoomStack == []:
                self.zoomState = (
                    self.graph.axisScale(QwtPlot.xBottom).lBound(),
                    self.graph.axisScale(QwtPlot.xBottom).hBound(),
                    self.graph.axisScale(QwtPlot.yLeft).lBound(),
                    self.graph.axisScale(QwtPlot.yLeft).hBound(),
                    )
        elif Qt.RightButton == e.button():
            self.zooming = 0
        # fake a mouse move to show the cursor position

    # onMousePressed()

    def onMouseReleased(self, e):
        if Qt.LeftButton == e.button():
            xmin = min(self.xpos, e.pos().x())
            xmax = max(self.xpos, e.pos().x())
            ymin = min(self.ypos, e.pos().y())
            ymax = max(self.ypos, e.pos().y())
            self.graph.setOutlineStyle(Qwt.Cross)
            xmin = self.graph.invTransform(QwtPlot.xBottom, xmin)
            xmax = self.graph.invTransform(QwtPlot.xBottom, xmax)
            ymin = self.graph.invTransform(QwtPlot.yLeft, ymin)
            ymax = self.graph.invTransform(QwtPlot.yLeft, ymax)
            if xmin == xmax or ymin == ymax:
                return
            self.zoomStack.append(self.zoomState)
            self.zoomState = (xmin, xmax, ymin, ymax)
            self.graph.enableOutline(0)
        elif Qt.RightButton == e.button():
            if len(self.zoomStack):
                xmin, xmax, ymin, ymax = self.zoomStack.pop()
            else:
                self.graph.setAxisAutoScale(QwtPlot.xBottom)
                self.graph.setAxisAutoScale(QwtPlot.yLeft)
                self.graph.replot()
                return

        self.graph.setAxisScale(QwtPlot.xBottom, xmin, xmax)
        self.graph.setAxisScale(QwtPlot.yLeft, ymin, ymax)
        self.graph.replot()

    def saveToFile(self):
        qfileName = QFileDialog.getSaveFileName("graph.png","Portable Network Graphics (.PNG)\nWindows Bitmap (.BMP)\nGraphics Interchange Format (.GIF)", None, "Save to..")
        fileName = str(qfileName)
        if fileName == "": return
        (fil,ext) = os.path.splitext(fileName)
        ext = ext.replace(".","")
        ext = ext.upper()
        cl = 0
        for g in self.graphs:
            if g.isVisible():
                clfname = fil + "_" + str(cl) + "." + ext
                g.saveToFileDirect(clfname, ext)
            cl += 1

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
        selected = 0
        for i in range(qlb.count()):
            if qlb.isSelected(i):
                selected = 1
                break
        qlb.selectAll(not(selected))

    def SUAclassQLB(self):
        self.selectUnselectAll(self.classQLB)
    ##

    ## class selection (classQLB)
    def classSelectionChange(self):
        list = []
        for i in range(self.classQLB.count()):
            if self.classQLB.isSelected(i):
                list.append( 1 )
            else:
                list.append( 0 )
        self.graph.setShowClasses(list)
    ##

    def calcGraph(self):
        self.progressBarInit()
        self.graph.setData(self.MAdata, self.classColor, self.classBrighterColor, self.ShowAverageProfile, self.ShowSingleProfiles, self.progressBarSet)
        self.graph.setPointWidth(self.PointWidth)
        self.graph.setCurveWidth(self.CurveWidth)
        self.graph.setAverageCurveWidth(self.AverageCurveWidth)

        self.graph.setAxisAutoScale(QwtPlot.xBottom)
        self.graph.setAxisAutoScale(QwtPlot.yLeft)
        self.graph.replot()
        self.progressBarFinished()

    def newdata(self):
        self.classQLB.clear()
        if self.MAdata.domain.classVar.varType <> orange.VarTypes.Discrete:
            print "error, class variable not discrete:", self.MAdata.domain.classVar
        if self.MAdata <> None and self.MAdata.domain.classVar.varType == orange.VarTypes.Discrete:
            ## classQLB
            self.numberOfClasses = len(self.MAdata.domain.classVar.values)
            self.classColor = ColorPaletteHSV(self.numberOfClasses, 160)
            self.classBrighterColor = ColorPaletteHSV(self.numberOfClasses, 255)

            self.calcGraph()
            ## update graphics
            ## classQLB
            self.classQVGB.show()
            classValues = self.MAdata.domain.classVar.values.native()
            for cn in range(len(classValues)):
                self.classQLB.insertItem(ColorPixmap(self.classBrighterColor[cn]), classValues[cn])
            self.classQLB.selectAll(1)  ##or: if numberOfClasses > 0: self.classQLB.setSelected(0, 1)
            self.update()

            if self.MAnoclass:
                self.classQVGB.hide()
        else:
            self.classColor = None
            self.classBrighterColor = None
        self.graph.show()
        self.layout.activate() # this is needed to scale the widget correctly

    def data(self, MAdata):
        if not MAdata:
            self.graph.hide()
            self.classQVGB.hide()
            return
        ## if there is no class attribute, create a dummy one
        if MAdata.domain.classVar == None:
            noClass = orange.EnumVariable('noclass', values=['none'])
            noClass.getValueFrom = lambda ex, w: 0
            newDomain = orange.Domain(MAdata.domain.variables + [noClass])
            self.MAdata = MAdata.select(newDomain)
            self.MAnoclass = 1 ## remember that there is no class to display
        else:
            self.MAdata = MAdata
            self.MAnoclass = 0 ## there are classes
        self.newdata()

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
    a.setMainWidget(owdm)
    d = orange.ExampleTable('wtclassed')
    owdm.data(d)
    owdm.show()
    a.exec_loop()
    owdm.saveSettings()
