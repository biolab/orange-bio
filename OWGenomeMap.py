"""
<name>Genome Map</name>
<description>Shows the locations of genes.</description>
<category>Genomics</category>
<icon>icons/GenomeMap.png</icon>
<priority>100</priority>
"""

import orange, OWGUI, math
import os.path # to find out where the local files are
from OWWidget import *
from qtcanvas import *

localdir = os.path.dirname(__file__) or "."

##############################################################################
# main class

# z coordinates for graphical objects
zchrom = 100; zchrombb=110; zgenes = 50; zticks = 40; zsel=10

class OWGenomeMap(OWWidget):
    settingsList = ["MinGeneWidth", "ShowTicks", "ColorByClass"]

    def __init__(self, parent=None, name='GenomeMap'):
        OWWidget.__init__(self, parent, name, "Shows the locations of genes on the chromosomes.")
        self.setWFlags(Qt.WResizeNoErase | Qt.WRepaintNoErase) #this works like magic.. no flicker during repaint!
        self.parent = parent        

        self.callbackDeposit = [] # deposit for OWGUI callback functions
        self.MinGeneWidth = 5
        self.ColorByClass = 1
        self.ShowTicks = 1
        self.loadSettings()

        self.classColors = []
        self.loadChromosomeDefinitions()
        self.data = None
        self.mid = orange.newmetaid() # meta id for a marker if a gene has a known position
        self.graph = ChromosomeGraph(self, self.chrom)
        
        # inputs and outputs
        self.inputs=[("Examples", ExampleTable, self.dataset, 1)]
        self.outputs = [("Examples", ExampleTable), ("Classified Examples", ExampleTableWithClass)]


        # GUI definition
        box = QVButtonGroup("Graph Options", self.controlArea)
        OWGUI.qwtHSlider(box, self, "MinGeneWidth", label='Min. mark width: ', labelWidth=80, minValue=1, maxValue=10, step=1, callback=self.graph.repaintGenes)
        self.colorByClassCB = OWGUI.checkBox(box, self, "ColorByClass", "Gene colors wrt class", callback=self.graph.repaintGenes, disabled=1)

        box=QBoxLayout(self.mainArea, QVBoxLayout.TopToBottom, 0)
        self.view = ChromosomeGraphView(self.graph, self.mainArea)
        self.view.setMinimumWidth(500)
        self.view.setMargins(0, max([x[1] for x in self.chrom]))
        box.addWidget(self.view)

    def loadChromosomeDefinitions(self):
##        import os
##        print 'PATH', os.getcwd()
##        self.chromData = orange.ExampleTable(localdir+'/test-genome.tab')
        self.chromData = orange.ExampleTable(localdir+'/Genome Maps/gene loci.tab')
        self.chrom = []
        self.geneCoordinates = {}
        id2desc = {}
        for d in self.chromData:
            if int(d['start']) == 0:
                self.chrom.append((int(d['start']), int(d['stop']), int(d['chromosome']), str(d['DDB'])))
                id2desc[int(d['chromosome'])] = len(self.chrom)-1
            else:
                self.geneCoordinates[str(d['DDB'])] = (int(d['start']), int(d['stop']), id2desc[int(d['chromosome'])] )

    def dataset(self, data):
        self.data = data

        if self.data:
            mid = self.mid
            ### XXX issue a warning if not found
            found = 0
            metas = self.data.domain.getmetas()
            for rmiIndx in metas.keys():
                if metas[rmiIndx].name == 'DDB':
                    found = 1
                    break
            if not found:
                print 'XXX warning: DDB not found in data set'
                return

            coord = self.geneCoordinates
            self.coord = [None] * len(data)
            if self.data:
                for (i,d) in enumerate(data):
    ##                rmi = str(d.getmeta(rmiIndx))    # XXX change to this once bug is removed
                    rmi = str(d['DDB'])
                    if coord.has_key(rmi):
                        self.coord[i] = coord[rmi]
                        d[mid] = 1
                    else:
                        ### XXX issue a warning
                        d[mid] = 0
                        print 'no key for', rmi
            else:
                pass

            self.colorByClassCB.setDisabled(data.domain.classVar == None)
            if data.domain.classVar:
                self.classColors = []
                for i in range(len(data.domain.classVar.values)):
                    newColor = QColor()
                    newColor.setHsv(i*360/len(data.domain.classVar.values), 255, 255)
                    self.classColors.append(newColor)

        self.graph.paint()
        
##############################################################################
# graph with chromosomes and genes

xoffset = 20; yoffset = 30; yspace = 50; ychrom = 20; ytick = 3
multiples = [(2., math.log10(2.)), (5., math.log10(5.))]
selectionColor = QColor(230,230,230)
gSelectionColor = QColor(190,190,190)
geneColor = QColor(60,60,60)
##geneColor = Qt.gray
##gSelectionColor = Qt.yellow

class ChromosomeGraph(QCanvas):
    def __init__(self, parent, chrom):
        apply(QCanvas.__init__,(self, parent, ""))
        self.parent = parent
        self.chrom = chrom
        self.selection = []

    def setMargins(self):
        view = self.parent.view
        self.bpL, self.bpR = view.margins[-1]
        self.bpW = float(self.bpR - self.bpL)
        self.ticks = []
        for c in self.chrom:
            self.ticks.append(self.getTicks(max(self.bpL, c[0]), min(self.bpR, c[1])))

    # converts bp index to canvas position
    def bp2x(self, bp):
        r = (bp - self.bpL) / self.bpW
        return int(xoffset + r * self.gwidth)

    def x2bp(self, x):
        r = (x - xoffset) / float(self.gwidth)
        return int(self.bpL + r * self.bpW)

    def paint(self):
        view = self.parent.view
        self.resize(view.width()-20, max(view.height()-5, yoffset+(len(self.chrom)+1)*yspace+yoffset))
        self.gwidth = self.width() - 2*xoffset

        # remove everything on present canvas
        for i in self.allItems():
            i.setCanvas(None)

        for (i, c) in enumerate(self.chrom):
            self.paintChrom(yoffset+i*(ychrom+yspace), max(self.bpL, c[0]), min(self.bpR, c[1]), i)
        self.paintGenes()
        self.paintSelection()
        self.update()
        
    def paintChrom(self, y, bpL, bpR, indx):
        chrom = self.chrom[indx]
        xL, xR = self.bp2x(bpL), self.bp2x(bpR)
        if xR-xL <= 0:
            return

        # paint the chromosome box
        # relW = (bpR-bpL)/float(self.bpR-self.bpL) # relative width of the displayed chromosome (0..1)
        # adjust the ticks to that

        r = QCanvasRectangle(xoffset, y, xR-xL+1, ychrom+1, self)
        r.setPen(QPen(Qt.white))
        r.y = y; r.id = indx
        r.setZ(zchrom)
        r.show()
        
        lu = QCanvasLine(self); lu.setPoints(0, 0, xR-xL, 0)
        ld = QCanvasLine(self); ld.setPoints(0, ychrom, xR-xL, ychrom)
        lines = [lu, ld]
        if bpL == chrom[0]:
            ll = QCanvasLine(self); ll.setPoints(0, 0, 0, ychrom)
            lines.append(ll)
        if bpR == chrom[1]:
            lr = QCanvasLine(self); lr.setPoints(xR-xL, 0, xR-xL, ychrom)
            lines.append(lr)
        for l in lines:
            l.setX(xoffset); l.setY(y); l.setZ(zchrombb)
            l.show()

        # paint chromosome name            
        label = QCanvasText(chrom[3], self)
        label.setX(xoffset); label.setY(y-label.boundingRect().height()-1); label.setZ(zticks)
        label.show()

        # paint the ticks
        if self.parent.ShowTicks:
            ticks = self.ticks[indx]
            for (bp, str) in ticks:
                x = self.bp2x(bp)
                tick = QCanvasLine(self); tick.setPoints(x, 0, x, ytick)
                tick.setX(0); tick.setY(y+ychrom); tick.setZ(zticks)
                tick.show()

                label = QCanvasText(str, self)
                label.setX(x - label.boundingRect().width()/2); label.setY(y+ychrom+ytick); label.setZ(zticks)
                label.show()

    # paint the genes
    def paintGenes(self):
        mid = self.parent.mid
        if not self.parent.data: return
        lborder, rborder = xoffset, self.width()-xoffset
        colorclass = self.parent.data.domain.classVar and self.parent.ColorByClass
        colors = self.parent.classColors
        for (i,d) in enumerate(self.parent.data):
            if not int(d[mid]):
                continue # position not known for this gene
            coord = self.parent.coord[i]
            if not(coord[0]>self.bpR or coord[1]<self.bpL):  # is gene in the visible area of the chromosome?
                l, r = max(coord[0], self.bpL),  min(coord[1], self.bpR)
                lp, rp = self.bp2x(l), self.bp2x(r)
                if rp-lp < self.parent.MinGeneWidth:
                    diff = int((self.parent.MinGeneWidth - (rp-lp)) / 2)
                    lp, rp = max(lp-diff, lborder), min(rp+diff, rborder)
                y = yoffset + coord[2]*(ychrom+yspace)
                r = QCanvasRectangle(lp, y, rp-lp, ychrom+1, self)
                if colorclass:
                    color = colors[int(d.getclass())]
                else:
                    color = geneColor
                r.instance = d
                r.setBrush(QBrush(color)); r.setPen(QPen(color))
                r.setZ(zgenes)
                r.show()

    def repaintGenes(self):
        for i in filter(lambda i,z=zgenes: i.z()==z, self.allItems()):
            i.setCanvas(None)
        self.paintGenes()
        self.update()

    def addSelection(self, x1, x2, cID, rect, replace=1):
        bp1 = self.x2bp(min(x1,x2)); bp2 = self.x2bp(max(x1,x2))
        if replace:
            self.selection = [(bp1, bp2, cID)]
        else:
            self.selection.append((bp1, bp2, cID))
        self.exportSelection()

    # finds which genes intersect with selection, makes and example table and sends it out
    def exportSelection(self):
        if not self.selection: return
        ex = []
        for i in filter(lambda i,z=zgenes: i.z()==zsel, self.allItems()):
            genes =  filter(lambda i,z=zgenes: i.z()==zgenes, self.collisions(i.boundingRect()))
            for g in genes:
                ex.append(g.instance)
        if len(ex):
            data = self.parent.data
            selectedData = orange.ExampleTable(data.domain, ex)

            # Reduce the number of class values, if class is defined
            cl = clo = data[0].domain.classVar
            if cl:
                cl = orange.RemoveUnusedValues(cl, ex, removeOneValued = 1)

            # Construct a new domain only if the class has changed
            # (ie to lesser number of values or to one value (alias None))
            if cl != clo:
                domain = orange.Domain(data[0].domain.attributes, cl)
                metas = data[0].domain.getmetas()
                for key in metas:
                    domain.addmeta(key, metas[key])
            else:
                domain = data[0].domain

            selectedData = orange.ExampleTable(domain, ex)
            if selectedData.domain.classVar:
                self.parent.send("Classified Examples", selectedData)
            else:
                self.parent.send("Classified Examples", None)
            self.parent.send("Examples", selectedData)


    def paintSelection(self):
        for (bp1, bp2, c) in self.selection:
            if not(bp1>self.bpR or bp2<self.bpL):
                l, r = max(bp1, self.bpL),  min(bp2, self.bpR)
                lp, rp = self.bp2x(l), self.bp2x(r)
                if rp<lp:
                    break
                r = QCanvasRectangle(lp, yoffset+c*(ychrom+yspace), rp-lp, ychrom, self)
                r.setZ(zsel)
                r.setBrush(QBrush(gSelectionColor)); r.setPen(QPen(gSelectionColor))
                r.show()
        
    def getTicks(self, lower, upper):
        ideal = (upper-lower)/2.
        if ideal<=0:
            return []
##        print 'iii %d (%d-%d)' % (ideal, upper, lower)
        log = math.log10(ideal)
        power = math.floor(log)
        fraction = log-power
        factor = 1.
        error = fraction
        for f, lf in multiples:
            e = math.fabs(fraction-lf)
            if e < error:
                error = e
                factor = f
        grid = factor * 10.**power
        digits = max(1, int(power))
        if power >= 6:
            format = '%dM'
            scale = 1000000
        elif power == 5:
            format = '%1.1fM'
            scale = 1000000
        elif power >= 3:
            format = '%dk'
            scale = 1000
        elif power == 2:
            format = '%1.1fk'
            scale = 1000
        else:
            format = '%d'
            scale = 1

        ticks = []
        t = -grid*math.floor(-lower/grid)
        while t <= upper and len(ticks) < 200:
            ticks.append((t, format % (t/scale,)))
            t = t + grid
        return ticks

##############################################################################
# the event manager

class ChromosomeGraphView(QCanvasView):
    def __init__(self, canvas, mainArea):
        apply(QCanvasView.__init__,(self,)+(canvas,mainArea))
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(1)
        self.setFocusPolicy(QWidget.ClickFocus)
        self.setFocus()
        self.margins = []
        self.selStart = None  # starting point of the current gene selection
        self.zoomStart = None # starting point for zooming
        self.shiftPressed = 0
        
    def resizeEvent(self, event):
        apply(QCanvasView.resizeEvent, (self,event))
        if self.canvas():
            self.canvas().paint()

    def setMargins(self, left, right):
        self.margins.append((left, right))
        self.canvas().setMargins()

    def retractMargins(self):
        if len(self.margins)>1:
            self.margins.pop(-1)
            self.canvas().setMargins()
            self.canvas().paint()

    def contentsMousePressEvent(self, event):
        if not self.canvas: return
        if event.button() == QMouseEvent.RightButton:
            self.retractMargins()
        elif event.button() == QMouseEvent.LeftButton:
            items = filter(lambda ci: ci.z()==zchrom, self.canvas().collisions(event.pos()))
            # items = self.canvas().collisions(event.pos())
            if items: # user pressed mouse on a chromosome
                self.addSelection = self.shiftPressed
                if not self.shiftPressed:
                    for i in filter(lambda i,z=zgenes: i.z()==zsel, self.canvas().allItems()):
                        i.setCanvas(None)
                self.chrom = items[0]
                self.selStart = event.pos().x()
                self.selRect = QCanvasRectangle(self.canvas())
                self.selRect.setY(self.chrom.y); self.selRect.setZ(zsel)
                self.selRect.setBrush(QBrush(gSelectionColor)); self.selRect.setPen(QPen(gSelectionColor))
                self.selRect.show()
                self.drawSel(self.selStart, self.selStart)
            else: # zooming
                self.zoomStart = event.pos().x()
                self.zoomRect = QCanvasRectangle(self.canvas())
                self.zoomRect.setY(0); self.zoomRect.setZ(0)
                self.zoomRect.setBrush(QBrush(selectionColor))
                self.zoomRect.setPen(QPen(selectionColor))
                self.zoomRect.show()
                self.drawZoom(self.zoomStart, self.zoomStart)

    def contentsMouseMoveEvent(self, event):
        x = event.pos().x()
        if self.zoomStart:
            self.drawZoom(self.zoomStart, x)
        elif self.selStart:
            self.drawSel(self.selStart, x)

    def contentsMouseReleaseEvent(self, event):
        x = event.pos().x()
        if self.zoomStart:
            left = self.canvas().x2bp(min(self.zoomStart, x))
            right = self.canvas().x2bp(max(self.zoomStart, x))
            if abs(right-left) < 50:
                self.zoomRect.setCanvas(None)
                self.canvas().update()
            else:
                self.setMargins(left, right)
                self.canvas().paint()
            self.zoomStart = None
        elif self.selStart:
            self.canvas().addSelection(self.selStart, x, self.chrom.id, self.selRect, replace=not self.addSelection)
            self.selStart = None

    def keyPressEvent(self, e):
        self.shiftPressed = e.key() == 4128

    def keyReleaseEvent(self, e):        
        self.shiftPressed = False            

    def drawZoom(self, left, right):
        self.zoomRect.setX(left)
        self.zoomRect.setSize(right-left, self.canvas().height())
        self.canvas().update()

    def drawSel(self, left, right):
        self.selRect.setX(left)
        self.selRect.setSize(right-left, ychrom)
        self.canvas().update()

##############################################################################
# test widget's appearance

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow=OWGenomeMap()
    a.setMainWidget(ow)
#    data = orange.ExampleTable("wtclassed.tab")
    data = orange.ExampleTable("hj.tab")
    ow.show()
    ow.dataset(data)
    a.exec_loop()
    # save settings
    ow.saveSettings()
