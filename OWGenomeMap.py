"""
<name>Genome Map</name>
<description>Shows the locations of genes.</description>
<category>Genomics</category>
<icon>icons/GenomeMap.png</icon>
<priority>200</priority>
"""

import orange, OWGUI, math, glob
import os.path # to find out where the local files are
from OWWidget import *
from qtcanvas import *

from OWGraph import ColorPaletteHSV

localdir = os.path.dirname(__file__) or "."

DEBUG = 0

##############################################################################
# main class

# z coordinates for graphical objects
zchrom = 100; zchrombb=110; zgenes = 50; zticks = 40; zsel=10

## chrom definition .tab file must have the following attributes (columns):
#   geneID
#   chromosome
#   start
#   stop
#
# chromosome definitions should be in the first lines
# chromosome definitions are those entries with column start == 0

class OWGenomeMap(OWWidget):
    settingsList = ["MinGeneWidth", "ShowTicks", "ColorByClass", "RecentGenomeMaps"]

    def __init__(self, parent=None, name='GenomeMap'):
        OWWidget.__init__(self, parent, name, "Shows the locations of genes on the chromosomes.")
        self.setWFlags(Qt.WResizeNoErase | Qt.WRepaintNoErase) #this works like magic.. no flicker during repaint!
        self.parent = parent        

        self.callbackDeposit = [] # deposit for OWGUI callback functions
        self.MinGeneWidth = 5
        self.ColorByClass = 1
        self.ShowTicks = 1
        self.RecentGenomeMaps = []
        self.GenomeMapIndx = 0
        self.genesInGenomeMapFile = {}
        self.GenomeMapLoaded = 0
        self.loadSettings()

        # received and decided based on input signal
        self.candidateGeneIDsFromSignal = [] ## list of discrete attributes present in clusterSet data signal
        self.geneIDattrIndx = 0 ## index of attribute in candidateGeneIDsFromSignal that was selected to represent the gene IDs
        self.geneIDattr = None ## self.geneIDattr is set accordingly

        # check if genome maps exist and remove those that don't
        # check that all files in directories "Genome Map" are included in the list
        self.RecentGenomeMaps = filter(os.path.exists, self.RecentGenomeMaps)
        widgetDir = os.path.dirname(os.path.abspath(__file__)) + "/"
        gmapList = glob.glob(widgetDir + 'Genome Map/*.tab')
        for f in gmapList:
            f = os.path.abspath(f)
            if f not in self.RecentGenomeMaps:
                self.RecentGenomeMaps.append( f)
        genesInRecentGenomeMapFile = {}
        #

        self.classColors = []
        self.geneCoordinates = {}
        self.data = None
        self.mid = orange.newmetaid() # meta id for a marker if a gene has a known position
        self.graph = ChromosomeGraph(self)

        # inputs and outputs
        self.inputs=[("Examples", ExampleTable, self.dataset, 1)]
        self.outputs = [("Examples", ExampleTable), ("Classified Examples", ExampleTableWithClass)]

        # GUI definition
        self.controls = QVGroupBox(self.controlArea)
        box = QVButtonGroup("Graph Options", self.controls)
        box.setMaximumSize(250, 80)
        OWGUI.qwtHSlider(box, self, "MinGeneWidth", label='Min. mark width: ', labelWidth=80, minValue=1, maxValue=10, step=1, callback=self.graph.repaintGenes)
        self.colorByClassCB = OWGUI.checkBox(box, self, "ColorByClass", "Gene colors wrt class", callback=self.graph.repaintGenes, disabled=1)

        box=QBoxLayout(self.mainArea, QVBoxLayout.TopToBottom, 0)
        self.view = ChromosomeGraphView(self.graph, self.mainArea)
        self.view.setMinimumWidth(500)
        box.addWidget(self.view)

        box = QHButtonGroup("Genome Map", self.controls)
        box.setMaximumSize(250, 50)
        self.genomeMapCombo = OWGUI.comboBox(box, self, 'GenomeMapIndx', items=[], callback=self.loadGenomeMap)
        self.genomeMapCombo.setMaximumSize(160, 20)
        self.setFilelist(self.genomeMapCombo, self.RecentGenomeMaps)
        self.genomeMapBrowse = OWGUI.button(box, self, 'Browse', callback=self.browseGenomeMap)
        self.genomeMapBrowse.setMaximumSize(50, 30)

        box = QHButtonGroup("Gene ID attribute", self.controls)
        box.setMaximumSize(250, 50)
        self.geneIDAttrCombo = OWGUI.comboBox(box, self, 'geneIDattrIndx', items=[], callback=self.geneIDchanged)
        self.geneIDAttrCombo.setMaximumSize(160, 20)
        self.setGeneIDAttributeList()

    def geneIDchanged(self):
        if len(self.candidateGeneIDsFromSignal) > self.geneIDattrIndx:
            self.geneIDAttrCombo.setCurrentItem(self.geneIDattrIndx)
            self.geneIDattr = self.candidateGeneIDsFromSignal[self.geneIDattrIndx]
        else:
            self.geneIDattr = None
        if DEBUG: print "changing geneID attribute to: " + str(self.geneIDattr)
        self.datasetChanged() ## recalculate the selected genes
        self.geneIDAttrCombo.setDisabled(len(self.candidateGeneIDsFromSignal) == 0)

    def setGeneIDAttributeList(self):
        ## refresh the list
        self.geneIDAttrCombo.clear()
        for f in self.candidateGeneIDsFromSignal:
            self.geneIDAttrCombo.insertItem(str(f.name))
        self.geneIDAttrCombo.setDisabled(len(self.candidateGeneIDsFromSignal) == 0)

    def setFilelist(self, filecombo, fileList):
        filecombo.clear()
        if fileList != []:
            for file in fileList:
                (dir, filename) = os.path.split(file)
                #leave out the path
                fnToDisp = filename
                filecombo.insertItem(fnToDisp)
            filecombo.setDisabled(False)
        else:
            filecombo.insertItem("(none)")
            filecombo.setDisabled(True)

    def loadChromosomeDefinitions(self, filename):
        self.geneCoordinates = {}
        chrom = []
        try:
            chromData = orange.ExampleTable(filename, dontCheckStored=1)
        except:
            self.graph.chrom = []
            self.GenomeMapLoaded = 0
            return
        id2desc = {}
        geneIDs = [] # all geneIDs in this file
        for d in chromData:
            geneID = str(d['geneID'])
            if int(d['start']) == 0: ## loading chromosomes definitions
                chrom.append((int(d['start']), int(d['stop']), int(d['chromosome']), geneID))
                id2desc[int(d['chromosome'])] = len(chrom)-1
            else: ## loading genes positions
                self.geneCoordinates[geneID] = (int(d['start']), int(d['stop']), id2desc[int(d['chromosome'])] )
            geneIDs.append( geneID)
        self.genesInGenomeMapFile[filename] = geneIDs ## update with new data (in case file has changed)
        self.GenomeMapLoaded = 1
        self.graph.chrom = chrom

    def repaintChromeGraph(self):
        self.view.resetMargins(0, max([10]+ [x[1] for x in self.graph.chrom]))
        self.graph.setMargins()
        self.graph.paint()

    def loadGenomeMap(self, change=1):
        if self.GenomeMapIndx < len(self.RecentGenomeMaps):
            fn = self.RecentGenomeMaps[self.GenomeMapIndx]
            if fn != "(none)":
                # remember the recent file list
                if fn in self.RecentGenomeMaps: # if already in list, remove it
                    self.RecentGenomeMaps.remove(fn)
                self.RecentGenomeMaps.insert(0, fn) # add to beginning of list
                self.setFilelist(self.genomeMapCombo, self.RecentGenomeMaps) # update combo
                self.loadChromosomeDefinitions(fn)
                if change: self.datasetChanged() ## repaint

    def browseGenomeMap(self):
        if self.RecentGenomeMaps == []:
            startfile = "."
        else:
            startfile = self.RecentGenomeMaps[0]
        filename = QFileDialog.getOpenFileName(startfile, 'Genome Map files (*.tab)\nAll files(*.*)', None, 'Genome Map File')
        fn = str(filename)
        fn = os.path.abspath(fn)
        if fn in self.RecentGenomeMaps: # if already in list, remove it
            self.RecentGenomeMaps.remove(fn)
        self.RecentGenomeMaps.insert(0, fn)
        self.GenomeMapIndx = 0
        self.loadGenomeMap()

    def findMostAppropriateGeneIDandGenomeMap(self):
        if self.data == None:
            self.candidateGeneIDsFromSignal = []
            self.geneIDattrIndx = -1
            self.geneIDattr = None
            self.setGeneIDAttributeList()
            return

        ## all discrete and string type attributes are good candidates
        self.candidateGeneIDsFromSignal = [a for a in self.data.domain.attributes + self.data.domain.getmetas().values() if a.varType == orange.VarTypes.Discrete or a.varType == orange.VarTypes.Other or a.varType == orange.VarTypes.String]
        self.setGeneIDAttributeList()
        self.geneIDAttrCombo.setDisabled(1)

        ## check if there are new genome map files present
        ## remove from geneID2genomeMapfile those not present in the RecentGenomeMaps list
        ## geneID is key, item is list of indexes in self.RecentGenomeMaps that have that geneID
        geneID2genomeMapfile = {}
        cn = 0
        for f in self.RecentGenomeMaps:
            if f not in self.genesInGenomeMapFile.keys():
                if DEBUG: print "loading", f
                try:
                    chromData = orange.ExampleTable(f, dontCheckStored=1)
                    geneIDs = [str(d['geneID']) for d in chromData] # all geneIDs in this file
                    self.genesInGenomeMapFile[f] = geneIDs # update with new data (in case file has changed)
                except:
                    self.genesInGenomeMapFile[f] = [] # update with new data (in case file has changed)
            for geneID in self.genesInGenomeMapFile[f]:
                tmpl = geneID2genomeMapfile.get(geneID, [])
                if cn not in tmpl:
                    tmpl.append(cn)
                    geneID2genomeMapfile[geneID] = tmpl
            cn += 1

        ## for each attribute look how many genesID are there, that are also present in geneID2genomeMapfile
        ## if current self.geneIDattr has count 0
        ## then select attribute with highest count
        ## else keep self.geneIDattr

        ## when best attribute selected, check if the loaded genome map is ok
        ## otherwise suggest the most appropriate genome map
        bestAttr = '' ## key is attribute, item is number of recognized geneIDs
        bestCn = 0
        bestGenomeMap = 0
        lst = self.candidateGeneIDsFromSignal
        if self.geneIDattr <> None and self.geneIDattr in self.candidateGeneIDsFromSignal: lst = [self.geneIDattr] + lst

        for attr in lst:
            vals = [ex[attr] for ex in self.data]

            ## calculate the frequency of each annotation file to which this geneID belongs to
            genomeMapFrequency = {}
            cn = 0
            for v in vals:
                v = str(v)
                i = geneID2genomeMapfile.get(v, -1) ## -1, not present
                if i <> -1:
                    for ai in i:
                        af = genomeMapFrequency.get(ai, 0)
                        genomeMapFrequency[ai] = af + 1
                    cn += 1
            if cn > bestCn or (cn > 0 and attr == self.geneIDattr):
                bestAttr = attr
                bestCn = cn
                gmfs = [(f, gmindex) for (gmindex, f) in genomeMapFrequency.items()]
                if len(gmfs) > 0:
                    gmfs.sort()
                    gmfs.reverse() ## most frequent first
                    bestGenomeMap = gmfs[0][1]
                else:
                    bestGenomeMap = 0 ## keep current
        if DEBUG: print "best attribute: " + str(bestAttr) + " with " + str(bestCn) + " gene IDs from genome map"
        if DEBUG: print "bestGenomeMap: " + str(self.RecentGenomeMaps[bestGenomeMap])

        self.geneIDattr = bestAttr
        try:
            self.geneIDattrIndx = self.candidateGeneIDsFromSignal.index(self.geneIDattr)
        except:
            self.geneIDattrIndx = 0

        ## load annotation if a better one found
        if bestGenomeMap <> 0 or not(self.GenomeMapLoaded):
            self.GenomeMapIndx = bestGenomeMap
            self.loadGenomeMap(0)
##            self.loadChromosomeDefinitions(self.RecentGenomeMaps[self.GenomeMapIndx])

        ## select the geneID, and rerun the GO term finding
        if DEBUG: print "geneID changed"
        self.geneIDchanged()

    def dataset(self, data):
        self.data = data
        self.findMostAppropriateGeneIDandGenomeMap() ## select most appropriate attribute only when first receiving the signal
        self.graph.selection = []
        self.datasetChanged()

    def datasetChanged(self):
        if self.geneIDattr == None:
            self.coord = []
            if self.data:
                for (i,d) in enumerate(self.data):
                    d[self.mid] = 0
            self.repaintChromeGraph() ## paint empty graph
            return

        ## in chrome data mark those records
        ## where the geneID matches that from geneCoordinates
        if self.data:
            ## make a self.coord the same size as input signal data
            self.coord = [None] * len(self.data)
            for (i,d) in enumerate(self.data):
                geneID = str(d[str(self.geneIDattr.name)])
                if self.geneCoordinates.has_key(geneID):
                    self.coord[i] = self.geneCoordinates[geneID]
                    d[self.mid] = 1
                else:
                    ### XXX issue a warning
                    d[self.mid] = 0
##                    print 'no key for', geneID

            ## create color map
            self.colorByClassCB.setDisabled(self.data.domain.classVar == None)
            if self.data.domain.classVar:
                self.classColors = ColorPaletteHSV(len(self.data.domain.classVar.values), 255)
            self.repaintChromeGraph()

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
    def __init__(self, parent, chrom = []):
        apply(QCanvas.__init__,(self, parent, ""))
        self.parent = parent
        self.chrom = chrom
        self.selection = []

    def setMargins(self):
        view = self.parent.view
        self.bpL, self.bpR = view.margins[-1]
        self.bpW = float(self.bpR - self.bpL)
        self.ticks = []
        ## find longest chrome
        for c in self.chrom:
            self.ticks.append(self.getTicks(max(self.bpL, c[0]), min(self.bpR, c[1]), self.bpL, self.bpR))

    # converts bp index to canvas position
    def bp2x(self, bp):
        r = (bp - self.bpL) / self.bpW
        return int(xoffset + r * self.gwidth)

    def x2bp(self, x):
        r = (x - xoffset) / float(self.gwidth)
        return int(self.bpL + r * self.bpW)

    def paint(self):
        view = self.parent.view
        self.resize(view.width()-20, max(view.height()-5, yoffset+(len(self.chrom)+1)*(ychrom+yspace) - yspace))
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
                r = QCanvasRectangle(lp, y, max(self.parent.MinGeneWidth, rp-lp), ychrom+1, self)
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
##            if cl:
##                cl = orange.RemoveUnusedValues(cl, ex, removeOneValued = 1)

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
        
    def getTicks(self, lower, upper, abslower, absupper):
        ideal = (upper-lower)/2.
        absideal = (absupper - abslower)/2
        ideal = max(absideal /4.0, ideal) ## don't display a too fine scale if the rest is quite big
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

    def resetMargins(self, left, right):
        self.margins = []
        self.setMargins(left, right)

    def setMargins(self, left, right):
        self.margins.append((max(0, left), right))
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
        if self.zoomStart <> None:
            self.drawZoom(self.zoomStart, x)
        if self.selStart <> None:
            self.drawSel(self.selStart, x)

    def contentsMouseReleaseEvent(self, event):
        x = event.pos().x()
        if self.zoomStart <> None:
            zoomStart = self.zoomStart ## remember is locally, and set it to None right away
            self.zoomStart = None
            left = self.canvas().x2bp(min(zoomStart, x))
            right = self.canvas().x2bp(max(zoomStart, x))
            if abs(right-left) < 50:
                self.zoomRect.setCanvas(None)
                self.canvas().update()
            else:
                self.setMargins(left, right)
                self.canvas().paint()
        if self.selStart <> None:
            selStart = self.selStart ## remember is locally, and set it to None right away
            self.selStart = None     ## otherwise the selection continues when other widgets process
            self.canvas().addSelection(selStart, x, self.chrom.id, self.selRect, replace=not self.addSelection)

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
