"""
<name>Epistasis Analysis</name>
<description>Epistasis analysis on microaray data</description>
<icon>icons/EpistasisAnalysis.png</icon>
<priority>1200</priority>
<prototype>1</prototype>
"""

from __future__ import absolute_import

import math

import orange, statc
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWWidget import *

from .OWDataFiles import DataFiles, ExampleSelection

##############################################################################
# parameters that determine the canvas layout

# constants for canvas drawing

canvasW = 220 # canvas width and height
canvasB = 30  # border
circleR = (canvasW - 2 * canvasB)/2.

##############################################################################
# main class

class OWEpistasisAnalysis(OWWidget):	
    settingsList = ['distype']

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'Epistasis Analysis') 
        
        self.inputs = [("Examples", ExampleTable, self.dataset, Multiple), ("Structured  Data", DataFiles, self.chipdata, 1)]
        self.outputs = [("Gene Selection", ExampleSelection), ("Examples A->B", ExampleTable), ("Examples B->A", ExampleTable), ("Examples A||B", ExampleTable)]

        self.data = []
        self.selectedFile = None
        self.A, self.B, self.AB = (None, None, None)
        self.epiA, self.epiB, self.para = (0, 0, 0)
        self.distype = 0
        #set default settings
        self.loadSettings()

        # Selection of data files
        box = OWGUI.widgetBox(self.controlArea, "Data files", addSpace=True)
#        box = QVButtonGroup("Data Files", self.controlArea)
#        self.fileLB = OWGUI.listBox(box, self, "selectedFile", "files", 
#                                    tooltip="Select example table",
##                                    callback=self.fileSelectionChanged,
#                                    selectionMode=QListWidget.SingleSelection,
#                                    )
#        self.connect(self.fileLB, SIGNAL("currentRowChanged(int)"), self.fileSelectionChanged)
        self.fileLB = QListWidget()
#        self.fileLB.setMaximumWidth(10)
        self.connect(self.fileLB, SIGNAL("highlighted(int)"), self.fileSelectionChanged)
#        self.connect(self.fileLB, SIGNAL("selected(int)"), self.setFileReferenceBySelection)
        box.layout().addWidget(self.fileLB)
        
        hbox = OWGUI.widgetBox(box, "", orientation="horizontal")
#        hbox = QHBox(box)
        self.markBtns = []
        for (i, lbl) in enumerate(["A","B","D"]):
            btn = OWGUI.button(hbox, self, lbl,
                               callback=lambda i=i: self.setMark(i),
                               disabled=1)
#            btn.setMaximumWidth(45)
#            btn.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed)
            self.markBtns.append(btn)

        # Relations, checkbox for gene selection
        self.sbox = OWGUI.widgetBox(self.controlArea, "Relations (Gene Selection)")
#        self.sbox = QVButtonGroup("Relations (Gene Selection)", self.controlArea)
        self.selChkbox = []
        self.cbinfo = (("epiA", "A -> B", "Gene B is epistatic to gene A"),
                       ("epiB", "B -> A", "Gene A is epistatic to gene B"),
                       ("para", "A || B", "Genes A and B are on parallel pathways"))
        for (i, cb) in enumerate(self.cbinfo):
            cb = OWGUI.checkBox(self.sbox, self, cb[0], cb[1],
                                tooltip=cb[2],
                                callback=self.filterSelectionChanged)
            self.selChkbox.append(cb)
        self.infochi = OWGUI.widgetLabel(self.sbox, "") #QLabel(self.sbox, '')
        self.sbox.setDisabled(1)

        OWGUI.radioButtonsInBox(self.controlArea, self, "distype", ["Average By Gene", "Average By Measurement"],
                                box="Distance", tooltips=None, callback=self.analysis)

        # scene
        self.scene = QGraphicsScene()
        self.sceneView = QGraphicsView()
        self.sceneView.setScene(self.scene)
        self.scene.setSceneRect(QRectF(0, 0, canvasW, canvasW))
        self.mainArea.layout().addWidget(self.sceneView)
        
        # canvas
#        self.canvas = QCanvas()
#        self.layout = QVBoxLayout(self.mainArea)
#        self.canvasView = QCanvasView(self.canvas, self.mainArea)
#        self.canvas.resize(canvasW, canvasW)
#        self.layout.add(self.canvasView)

        self.resize(420,350)

    ##########################################################################
    # handling of input/output signals

    def checkDomain(self, data, selection = None):
        
        cl = clo = data.domain.classVar
        if cl:
            if selection:
                cl = orange.RemoveUnusedValues(cl, selection, removeOneValued = 1)
            else:
                cl = orange.RemoveUnusedValues(cl, data, removeOneValued = 1)

        # Construct a new domain only if the class has changed
        # (ie to lesser number of values or to one value (alias None))
        if cl != clo:
            domain = orange.Domain(data.domain.attributes, cl)
            metas = data.domain.getmetas()
            for key in metas:
                domain.addmeta(key, metas[key])
            return domain
        else:
            return None
        
    def dataset(self, data, id):
        ids = [d.id for d in self.data]
        if not data:
            if id in ids:
                k = ids.index(id)
                del self.data[k]
                self.fileLB.takeItem(k)
#                self.fileLB.removeItem(k)
        else:
            # check if the same length
            if data.domain.classVar:
                domain = self.checkDomain(data)
                if domain:
                    data = orange.ExampleTable(domain, data)
            data.setattr("id", id)
            data.setattr("marker", None)
            if id in ids:
                data.id = id
                indx = ids.index(id)
                self.data[indx] = data
                self.fileLB.takeItem(indx)
                self.fileLB.insertItem(indx, self.createListItem(data))
#                self.fileLB.changeItem(self.createListItem(data), indx)
            else:
                if len(self.data) < 3:
                    data.setattr("marker", len(self.data)) #### REMOVE !!!
                self.fileLB.addItem(self.createListItem(data))
                self.data.append(data)
                
                if len(self.data) >= 3: #### REMOVE !!!
                    self.analysis()

    def chipdata(self, data):
        self.data = [] # XXX should only remove the data from the same source, use id in this rutine
        self.fileLB.clear()
        if not data:
            for i in self.canvas.allItems():
                i.setCanvas(None)
            self.canvas.update()
            return
        indx = 0
        for (strainname, ds) in data:
            for d in ds:
                self.dataset(d, indx)
                indx += 1

    ##########################################################################
    # handle events

    # mark button pressed, mark the file appropriately
    def setMark(self, mark):
        sel = self.selectedFile
        # remove this mark from previously marked file
        marks = [d.marker for d in self.data]
        if mark in marks:
            indx = marks.index(mark)
            self.data[indx].setattr("marker", None)
            self.fileLB.takeItem(indx)
            self.fileLB.insertItem(indx, self.createListItem(self.data[indx]))
#            self.fileLB.changeItem(self.createListItem(self.data[indx]), indx)
        self.data[sel].setattr("marker", mark)
        self.fileLB.takeItem(sel)
        self.fileLB.insertItem(sel, self.createListItem(self.datap[sel]))
#        self.fileLB.changeItem(self.createListItem(self.data[sel]), sel)
        self.analysis() # do the analysis (if all three selections are made)

    # user has selected a different file from the file list
    def fileSelectionChanged(self, sel):
        self.selectedFile = sel
        for btn in self.markBtns:
            btn.setEnabled(True)

    # user has selected a different file
    def filterSelectionChanged(self):
        pass
        
    def createListItem(self, data):
        pixmap = QPixmap(QSize(14, 14))
#        pixmap.resize(14,13)
        pixmap.fill(Qt.white)
        
        if data.marker is not None:
            painter = QPainter()
            painter.begin(pixmap)
            painter.setPen(Qt.black)
            painter.setBrush(Qt.white)
            painter.drawRect(0, 0, 13, 13)
            c = ['A','B','D'][data.marker]
            painter.drawText(3, 11, c)
            painter.end()
            
        listItem = QListWidgetItem(QIcon(pixmap), data.name)
#        listItem = QListBoxPixmap(pixmap)
#        listItem.setText(data.name)
        return listItem
##        return QListBoxText(text)

    ##########################################################################
    # epistasis analysis

    import math

    def analysis(self):
        markers = [d.marker for d in self.data]
        if len(filter(lambda x: x is not None, markers)) < 3:
            self.sbox.setDisabled(True)
            for item in self.scene.items():
                self.scene.removeItem(item)
#                i.setCanvas(None)
#            self.canvas.update()
            return
        self.sbox.setEnabled(True)

        pa, pb, pab = [self.data[markers.index(x)] for x in range(3)]
        
        pb = orange.ExampleTable(pa.domain, pb)
        pab = orange.ExampleTable(pa.domain, pab)
        
        dist = orange.ExamplesDistanceConstructor_Euclidean(pa, normalize=False)

        ave = [0]*3
        vote = [0]*3
        genevote = []
        for g in range(len(pa)):
            d = [dist(pb[g], pab[g]), dist(pa[g], pab[g]), dist(pa[g], pb[g])]
            voteindx = d.index(min(d))
            vote[voteindx] += 1
            genevote.append(voteindx)
            if self.distype == 1:
                for i in range(3):
                    d[i] = d[i] * d[i]
            for i in range(3):
                ave[i] += d[i]
        if self.distype == 1:
            ave = [math.sqrt(x) / len(pa) for x in ave]
        else:
            ave = [x/len(pa) for x in ave]

        # compute Chi^2 statistics,
        # update the interface (report on results)
        for i in range(3):
            self.selChkbox[i].setText(self.cbinfo[i][1] + "  (%d genes)" % vote[i])
        p = statc.chisquare([len(pa) / 3.] * 3, vote)[1]
        self.infochi.setText('Chi Square: ' + ['p = %6.4f' % p, 'p < 0.0001'][p < 0.0001])

        self.setAnalysisPlot(ave)
        self.senddata(genevote)

    def setAnalysisPlot(self, ds):
        def plotDot(x, y, w=10, z=10):
            dot = QGraphicsEllipseItem(x, y, w, w, self.scene)
#            dot = QCanvasEllipse(w, w, self.canvas)
            dot.setBrush(QBrush(Qt.black))
            dot.setZValue(z)
#            dot.setX(x); dot.setY(y); dot.setZ(z)
            dot.show()
            
        def plotLine(x0, y0, x1, y1, w=2, z=10, color=Qt.black):
            line = QGraphicsLineItem(x0, y0, x1, y1, self.scene)
#            line = QCanvasLine(self.canvas)
#            line.setPoints(x0, y0, x1, y1)
            line.setPen(QPen(color, w))
#            line.setZ(10)
            line.setZValue(10)
            line.show()
            
        def plotText(x, y, text, xoffset=0):
            t = QGraphicsSimpleTextItem("", self.scene)
#            t = QCanvasText(self.canvas)
#            t.setText(text)
            xw = t.boundingRect().width()
            t.setPos(QPointF(x+xw*xoffset, y))
            t.setZValue(20)
#            t.setX(x+xw*xoffset); t.setY(y); t.setZ(20)
            t.show()        

        for item in self.scene.items():
            self.scene.removeItem(item)

        s = sum(ds)/2.
        K = math.sqrt( s * reduce(lambda x,y: x*y, map(lambda x: s-x, ds)) )
        R = reduce(lambda x,y: x*y, ds) / (4 * K)
        scale = circleR / R
        yab = canvasB + circleR - scale * math.sqrt(R**2 - (ds[2]/2)**2)
        xa = canvasB + circleR - scale * ds[2]/2
        xb = canvasB + circleR + scale * ds[2]/2
        h = 2*K/ds[2]
        yd = yab + scale * h
        xd = xa + scale * math.sqrt(max(ds[0], ds[1])**2 - h**2)

        # plot a circle
        c = QGraphicsEllipseItem(0, 0, circleR*2, circleR*2, self.scene)
#        c = QCanvasEllipse(circleR*2, circleR*2, self.canvas)
        c.setBrush(QBrush(QColor(240,240,240)))
        c.setPos(canvasB+circleR, canvasB+circleR)
        c.setZValue(0)
#        c.setX(canvasB+circleR); c.setY(canvasB+circleR); c.setZ(0)
        c.show()

        # plot triangle dots, line, ...
        plotDot(xa, yab)
        plotDot(xb, yab)
        plotDot(xd, yd)

        plotLine(xa, yab, xb, yab, color=Qt.red)
        plotLine(xa, yab, xd, yd)
        plotLine(xb, yab, xd, yd)

        plotText(xa+(xb-xa)/2, yab-15, "%6.4f" % ds[2], xoffset=-0.5)
        ymid = yab+(yd-yab)/2
        plotText(xa+(xd-xa)/2-3, ymid, "%6.4f" % max(ds[0], ds[1]), xoffset=-1)
        plotText(xd+(xb-xd)/2+3, ymid, "%6.4f" % min(ds[0], ds[1]))

        # edge labels
        markers = [d.marker for d in self.data]
        labels = [self.data[markers.index(x)].name for x in range(3)]
        plotText(xd, yd+8, labels[2], xoffset=-0.5)
        plotText(xa, yab-20, labels[ds[:-1].index(min(ds[:-1]))], xoffset=-0.67)
        plotText(xb, yab-20, labels[ds[:-1].index(max(ds[:-1]))], xoffset=-0.33)

        self.canvas.update()

    def senddata(self, genevotes):
        def sendkeyeddata(channel, key):
            for i in range(3):
                d = datasets[i].select(genevotes, key)
                d.name = datasets[i].name
                self.send(channel, d, i)
        markers = [d.marker for d in self.data]
        datasets = [self.data[markers.index(x)] for x in range(3)]
        channels= ["Examples A->B", "Examples B->A", "Examples A||B"]
        for c in channels:
            for i in range(3): # this should be excluded, repear in heat map
                self.send(c, None, i)
        for (i,ch) in enumerate(channels):
            sendkeyeddata(ch, i)

##################################################################################################
# test script

if __name__=="__main__":
    import orange
    a = QApplication(sys.argv)
    ow = OWEpistasisAnalysis()
    a.setMainWidget(ow)

    ow.show()
    names = ['wt1', 'wt2', 'wt3', 'wt4']
##    names = [r'chipdata/pufA/pufA1.1.raw.tab', r'chipdata/yakA/yakA1.1.raw.tab', r'chipdata/yakApufA/yakApufA1.1.raw.tab', r'chipdata/yakApufA/yakApufA1.1.raw.tab']
    for i, s in enumerate(names): 
        d = orange.ExampleTable(s); d.name = s
        ow.dataset(d, i)
##    ow.dataset(None, 2)

    a.exec_loop()
    ow.saveSettings()
