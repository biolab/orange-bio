"""
<name>Heat Map</name>
<description>Microarray Heat Map</description>
<icon>icons/HeatMap.png</icon>
<priority>50</priority>
"""

import orange, math
import OWGUI
#from string import *
from qt import *
from qtcanvas import *
from OWWidget import *
from OWOptions import *
from qwt import *
from OWChipDataFiles import ChipData

from OWChipANOVA import ANOVAResults

##############################################################################
# parameters that determine the canvas layout

c_offsetX = 10; c_offsetY = 10  # top and left border
c_spaceX = 10; c_spaceY = 10    # space btw graphical elements
c_legendHeight = 15             # height of the legend
c_averageStripeWidth = 12       # width of the stripe with averages

z_heatmap = 5                   # layer with heatmaps

##############################################################################
# main class

class OWHeatMap(OWWidget):	
    settingsList = ["CellWidth", "CellHeight", "Merge", "Gamma", "CutLow", "CutHigh", "CutEnabled", 
                    "ShowAnnotation", "LegendOnTop", "LegendOnBottom", "ShowAverageStripe", "ShowGroupLabel",
                    "MaintainArrayHeight",
                    'BShowballoon', 'BShowColumnID', 'BShowSpotIndex', 'BShowAnnotation', 'BShowGeneExpression',
                    "BSpotVar", "ShowGeneAnnotations", "ShowDataFileNames", "BAnnotationVar",
                    "SelectionType",
                    "CurrentPalette", "SortGenes"]

    def __init__(self, parent=None, name='HeatMap'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, name, 'Microarray Heat Map', FALSE, FALSE) 
        
        self.inputs = [("Examples", ExampleTable, self.dataset, 0), ("Structured Chip Data", ChipData, self.chipdata, 1)]
        self.outputs = [("Examples", ExampleTable), ("Classified Examples", ExampleTableWithClass)]

        #set default settings
        self.CellWidth = 3; self.CellHeight = 3
        self.Merge = 1; self.savedMerge = self.Merge
        self.Gamma = 1
        self.CutLow = 0; self.CutHigh = 0; self.CutEnabled = 0
        self.ShowAnnotation = 0
        self.LegendOnTop = 1           # legend stripe on top (bottom)?
        self.LegendOnBottom = 1
        self.ShowGroupLabel = 1        # show class names in case of classified data?
        self.ShowAverageStripe = 0     # show the stripe with the evarage
        self.MaintainArrayHeight = 1   # adjust cell height while changing the merge factor
        self.BShowballoon = 1          # balloon help
        self.ShowGeneAnnotations = 1   # show annotations for genes
        self.ShowDataFileNames = 1     # show the names of the data sets (only in case of multiple files)
        self.BShowColumnID = 1; self.BShowSpotIndex = 1; self.BShowAnnotation = 1; self.BShowGeneExpression = 1
        self.BSpotVar = None; self.BAnnotationVar = None  # these are names of variables
        self.BSpotIndx = None; self.BAnnotationIndx = None # these are id's of the combo boxes
        self.SortGenes = 1
        self.SelectionType = 0         # selection on a single data set
        self.setColorPalette()
        self.refFile = 0               # position index of a reference file
        self.selectedFile = None       # position of the selected file in the list box
        
        self.loadSettings()
        self.data = []
        self.maxHSize = 15; self.maxVSize = 15

        # GUI definition
        self.tabs = QTabWidget(self.controlArea, 'tabWidget')

        # SETTINGS TAB
        settingsTab = QVGroupBox(self)
        box = QVButtonGroup("Cell Size (Pixels)", settingsTab)
        OWGUI.qwtHSlider(box, self, "CellWidth", label='Width: ', labelWidth=38, minValue=1, maxValue=self.maxHSize, step=1, precision=0, callback=self.drawHeatMap)
        self.sliderVSize = OWGUI.qwtHSlider(box, self, "CellHeight", label='Height: ', labelWidth=38, minValue=1, maxValue=self.maxVSize, step=1, precision=0, callback=self.createHeatMap)
        OWGUI.qwtHSlider(settingsTab, self, "Gamma", box="Gamma", minValue=0.1, maxValue=1, step=0.1, callback=self.drawHeatMap)

        # define the color stripe to show the current palette
        colorItems = [self.createColorStripe(i) for i in range(len(self.ColorPalettes))]
        palc = OWGUI.comboBox(settingsTab, self, "CurrentPalette", box="Colors", items=None, tooltip=None, callback=self.setColor)
        for cit in colorItems:
            palc.insertItem(cit) ## because of a string cast in the comboBox constructor
        OWGUI.checkBox(settingsTab, self, "SortGenes", "Sort genes", box="Sort", callback=self.constructHeatmap)

        self.tabs.insertTab(settingsTab, "Settings")

        # FILTER TAB
        tab = QVGroupBox(self)
        box = QVButtonGroup("Threshold Values", tab)
        OWGUI.checkBox(box, self, 'CutEnabled', "Enabled", callback=self.setCutEnabled)
        self.sliderCutLow = OWGUI.qwtHSlider(box, self, 'CutLow', label='Low:', labelWidth=33, minValue=-100, maxValue=0, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.drawHeatMap)
        self.sliderCutHigh = OWGUI.qwtHSlider(box, self, 'CutHigh', label='High:', labelWidth=33, minValue=0, maxValue=100, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.drawHeatMap)
        if not self.CutEnabled:
            self.sliderCutLow.box.setDisabled(1)
            self.sliderCutHigh.box.setDisabled(1)

        box = QVButtonGroup("Merge", tab)
        OWGUI.qwtHSlider(box, self, "Merge", label='Rows:', labelWidth=33, minValue=1, maxValue=100, step=1, callback=self.mergeChanged, ticks=0)
        OWGUI.checkBox(box, self, 'MaintainArrayHeight', "Maintain array height")

        self.tabs.insertTab(tab, "Filter")

        # INFO TAB
        tab = QVGroupBox(self)

        box = QVButtonGroup("Annotation && Legends", tab)
        OWGUI.checkBox(box, self, 'LegendOnTop', 'Show legend', callback=self.drawHeatMap)
        OWGUI.checkBox(box, self, 'ShowAverageStripe', 'Stripes with averages', callback=self.drawHeatMap)
        self.geneAnnotationsCB = OWGUI.checkBox(box, self, 'ShowGeneAnnotations', 'Gene annotations', callback=self.drawHeatMap)
        
        self.annotationCombo = OWGUI.comboBox(box, self, "BAnnotationIndx", items=[], callback=lambda x='BAnnotationVar', y='BAnnotationIndx': self.setMetaID(x, y))

        box = QVButtonGroup("Balloon", tab)
        OWGUI.checkBox(box, self, 'BShowballoon', "Show balloon", callback=lambda: self.balloonInfoBox.setDisabled(not self.BShowballoon))
        box = QVButtonGroup("Balloon Info", tab)
        OWGUI.checkBox(box, self, 'BShowColumnID', "Column ID")
        self.spotIndxCB = OWGUI.checkBox(box, self, 'BShowSpotIndex', "Spot Index", callback=lambda: self.spotCombo.setDisabled(not self.BShowSpotIndex))
        self.spotCombo = OWGUI.comboBox(box, self, "BSpotIndx", items=[], callback=lambda x='BSpotVar', y='BSpotIndx': self.setMetaID(x, y))
        OWGUI.checkBox(box, self, 'BShowGeneExpression', "Gene expression")
        OWGUI.checkBox(box, self, 'BShowAnnotation', "Annotation")
        self.balloonInfoBox = box
        self.tabs.insertTab(tab, "Info")

        # FILES TAB
        self.filesTab = QVGroupBox(self)
        box = QVButtonGroup("Data Files", self.filesTab)
        self.fileLB = QListBox(box, "lb")
##        self.fileLB.setMaximumWidth(10)
        self.connect(self.fileLB, SIGNAL("highlighted(int)"), self.fileSelectionChanged)
        self.connect(self.fileLB, SIGNAL("selected(int)"), self.setFileReferenceBySelection)
        self.tabs.insertTab(self.filesTab, "Files")
        self.tabs.setTabEnabled(self.filesTab, 0)
        hbox = QHBox(box)
        self.fileUp = OWGUI.button(hbox, self, 'Up', callback=lambda i=-1: self.fileOrderChange(i), disabled=1)
        self.fileRef = OWGUI.button(hbox, self, 'Ref', self.setFileReference, disabled=1)
        self.fileDown = OWGUI.button(hbox, self, 'Down', callback=lambda i=1: self.fileOrderChange(i), disabled=1)
        for btn in [self.fileUp, self.fileRef, self.fileDown]:
            btn.setMaximumWidth(45)
            

        OWGUI.checkBox(self.filesTab, self, 'ShowDataFileNames', 'Show data file names', callback=self.drawHeatMap)
        OWGUI.radioButtonsInBox(self.filesTab, self, 'SelectionType', ['Single data set', 'Multiple data sets'], box='Selection')

        self.resize(700,400)

        # canvas with microarray
        self.layout = QVBoxLayout(self.mainArea)
        self.canvas = QCanvas()
        self.canvasView = MyCanvasView(self.canvas, self.mainArea)
        self.selection = SelectData(self, self.canvas)
        # self.canvasView = QCanvasView(self.canvas, self.mainArea)
        self.layout.add(self.canvasView)

    def createColorStripe(self, palette):
        dx = 104; dy = 18
        bmp = chr(252)*dx*2 + reduce(lambda x,y:x+y, [chr(i*250/dx) for i in range(dx)] * (dy-4)) + chr(252)*dx*2 
        image = QImage(bmp, dx, dy, 8, self.ColorPalettes[palette], 256, QImage.LittleEndian)
        pm = QPixmap()
        pm.convertFromImage(image, QPixmap.Color);
        return pm

    # set the default palettes used in the program
    # palette defines 256 colors, 250 are used for heat map, remaining 6 are extra
    # color indices for unknown is 255, underflow 253, overflow 254, white 252
    def setColorPalette(self):
        white = qRgb(255,255,255)
        gray = qRgb(200,200,200)
        self.ColorPalettes = \
          ([qRgb(255.*i/250., 255.*i/250., 255-(255.*i/250.)) for i in range(250)] + [white]*3 + [qRgb(0., 0., 255.), qRgb(255., 255., 0.), gray],
           [qRgb(0, 255.*i*2/250., 0) for i in range(125, 0, -1)] + [qRgb(255.*i*2/250., 0, 0) for i in range(125)] + [white]*3 + [qRgb(0, 255., 0), qRgb(255., 0, 0), gray],
           [qRgb(255.*i/250., 0, 0) for i in range(250)] + [white]*3 + [qRgb(0., 0, 0), qRgb(255., 0, 0), gray])
        self.SelectionColors = [QColor(0,0,0), QColor(255,255,128), QColor(0,255,255)]
        self.CurrentPalette = 0
        
    # any time the data changes, the two combo boxes showing meta attributes
    # have to be adjusted
    def setMetaCombo(self, cb, value, enabled=1, default=None):
        cb.clear()
        if len(self.meta)==0:
            cb.setDisabled(True)
            self.spotIndxCB.setDisabled(1); self.geneAnnotationsCB.setDisabled(1)
            return (None, None)
        cb.setDisabled(not enabled)
        self.spotIndxCB.setEnabled(1); self.geneAnnotationsCB.setEnabled(1)
        for m in self.meta:
            cb.insertItem(m)
        
        if not (value in self.meta):
            if default in self.meta:
                value = default
            else:
                value = None

        if value in self.meta:
            cb.setCurrentItem(self.meta.index(value))
            indx = self.meta.index(value)
        else:
            cb.setCurrentItem(0)
            value = self.meta[0]; indx = 0
        return (value, indx)

    def setMetaID(self, val, valIndx):
        setattr(self, val, self.meta[getattr(self, valIndx)])
        if val=='BAnnotationVar':
            self.drawHeatMap()

    def setMetaCombos(self):
        self.meta = [m.name for m in self.data[0].domain.getmetas().values()]
        self.BSpotVar, self.BSpotIndx = self.setMetaCombo(self.spotCombo, self.BSpotVar, enabled=self.BShowSpotIndex, default='RMI')
        self.BAnnotationVar, self.BAnnotationIndx = self.setMetaCombo(self.annotationCombo, self.BAnnotationVar, enabled=self.BShowAnnotation, default='xannotation')

    ##########################################################################
    # handling of input/output signals

    def dataset(self, data, id):
        print 'GOT', id
        ids = [d.id for d in self.data]
        if not data:
            if id in ids:
                k = ids.index(id)
                del self.data[k]
                self.fileLB.removeItem(k)
                if self.refFile == k:
                    self.refFile = 0
                    if len(self.data):
                        self.fileLB.changeItem(self.createListItem(self.data[0].name, self.refFile), self.refFile)
        else:
            # check if the same length
            if data.domain.classVar:
                domain = self.checkDomain(data)
                if domain:
                    data = orange.ExampleTable(domain, data)
            data.setattr("id", id)
            if id in ids:
                data.id = id
                indx = ids.index(id)
                self.data[indx] = data
                self.fileLB.changeItem(self.createListItem(data.name, indx), indx)
            else:
                self.fileLB.insertItem(self.createListItem(data.name, len(self.data)))
                self.data.append(data)
#                self.fileLB.insertItem(data.name)

            if len(self.data) > 1:
                self.tabs.setTabEnabled(self.filesTab, 1)
            else:
                self.tabs.setTabEnabled(self.filesTab, 0)
            self.setMetaCombos() # set the two combo widgets according to the data
            
        self.send('Classified Examples', None)
        self.send('Examples', None)
        self.constructHeatmap()
        self.canvas.update()

    def chipdata(self, data):
        print 'CHIPDATA'
        self.data = [] # XXX should only remove the data from the same source, use id in this rutine
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
        self.createHeatMap()
        self.canvas.update()
        
    def constructHeatmap(self):
        if len(self.data):
            self.heatmapconstructor = [None] * len(self.data)
            if self.SortGenes:
                self.heatmapconstructor[self.refFile] = orange.HeatmapConstructor(self.data[self.refFile])
            else:
                self.heatmapconstructor[self.refFile] = orange.HeatmapConstructor(self.data[self.refFile], None)
            for i in range(len(self.data)):
                if i <> self.refFile:
                    self.heatmapconstructor[i] = orange.HeatmapConstructor(self.data[i], self.heatmapconstructor[self.refFile])
        else:
            self.heatmapconstructor = []
        self.createHeatMap()

    # remove unused values from the class of the data set
    def checkDomain(self, data, selection = None):
        # Reduce the number of class values, if class is defined
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

    # send out the data for selected rows, rows = [(group, from, to), ...]
    def sendData(self, rows):
        ex = []
        for (g,s,e) in rows:
            hm = self.heatmaps[0][g]
            ex += hm.examples[hm.exampleIndices[s] : hm.exampleIndices[e+1]]

        # Reduce the number of class values, if class is defined
        newdomain = self.checkDomain(self.data[0], selection=ex)
        if not newdomain:
            newdomain = self.data[0].domain
        selectedData = orange.ExampleTable(newdomain, ex)
        if selectedData.domain.classVar:
            self.send("Classified Examples", selectedData)
        else:
            self.send("Classified Examples", None)
        self.send("Examples", selectedData, 1)

    ##########################################################################
    # callback functions

    def setColor(self):
        if self.CurrentPalette == len(self.ColorPalettes):
            self.CurrentPalette = 0
            # put a code here that allows to define ones own colors
        else:
            pm = self.createColorStripe(self.CurrentPalette)
            self.drawHeatMap()

    def setCutEnabled(self):
        self.sliderCutLow.box.setDisabled(not self.CutEnabled)
        self.sliderCutHigh.box.setDisabled(not self.CutEnabled)
        self.drawHeatMap()

    def mergeChanged(self):
        self.oldMerge = self.savedMerge
        if self.MaintainArrayHeight and self.oldMerge <> self.Merge:
            k = self.Merge / self.oldMerge
            l = max(1, min(self.CellHeight * k, self.maxVSize))
            if l <> self.CellHeight:
                self.CellHeight = l
                self.sliderVSize.setValue(self.CellHeight)

        self.createHeatMap()
        self.savedMerge = self.Merge

    def fileOrderChange(self, chg):
        if chg==-1 and self.selectedFile>0:
            switchFiles(self.selectedFile, self.selectedFile-1)
        if chg==1  and self.selectedFile < len(self.data - 1):
            switchFiles(self.selectedFile, self.selectedFile+1)
        

    ##########################################################################
    # drawing

    def drawLegend(self, x, y, width, height, palette):
        legend = self.heatmapconstructor[0].getLegend(width, height, self.Gamma)

        lo = self.CutEnabled and self.CutLow   or self.lowerBound
        hi = self.CutEnabled and self.CutHigh  or self.upperBound

        t = QCanvasText("%3.1f" % lo, self.canvas)
        t.setX(x); t.setY(y)
        t.show()
        t = QCanvasText("%3.1f" % hi, self.canvas)
        t.setX(x+width-t.boundingRect().width()); t.setY(y)
        t.show()
        y += t.boundingRect().height()+1
        self.legendItem = ImageItem(legend, self.canvas, width, height, palette, x=x, y=y)
        return y + c_legendHeight + c_spaceY

    def drawFileName(self, label, x, y, width):
        t = QCanvasText(label, self.canvas)
        t.setX(x); t.setY(y)
        t.show()
        line = QCanvasLine(self.canvas)
        line.setPoints(0, 0, width, 0)
        y += t.boundingRect().height()
        line.setX(x); line.setY(y)
        line.show()
        return y + 5

    def drawGroupLabel(self, label, x, y, width):
        t = QCanvasText(label, self.canvas)
        t.setX(x); t.setY(y)
        t.show()
        return y + t.boundingRect().height() + 1

    def drawGeneAnnotation(self, x, y, group):
        # determine the appropriate font width for annotation
        # this part is ugly, we need to do it computationally
        font = QFont()
        dummy = QCanvasText("dummy", self.canvas)
        last = 2; offset = 0
        for fsize in range(2,9):
            font.setPointSize(fsize)
            dummy.setFont(font)
            if dummy.boundingRect().height() > self.CellHeight:
                break
            offset = (self.CellHeight - dummy.boundingRect().height())/2
            last = fsize
        font.setPointSize(last)
        y += offset

        # annotate
        hm = self.heatmaps[0][group]
        if self.BAnnotationVar:
            for (row, indices) in enumerate(hm.exampleIndices[:-1]):
                t = QCanvasText(str(hm.examples[hm.exampleIndices[row]][self.BAnnotationVar]), self.canvas)
                t.setFont(font)
                t.setX(x); t.setY(y)
                t.show()
                y += self.CellHeight

    def drawHeatMap(self):
        # remove everything from current canvas
        for i in self.canvas.allItems():
            i.setCanvas(None)
        if not len(self.data):
            return

        lo = self.CutEnabled and self.CutLow   or self.lowerBound
        hi = self.CutEnabled and self.CutHigh  or self.upperBound

        self.canvasView.heatmapParameters(self, self.CellWidth, self.CellHeight) # needed for event handling

        palette = self.ColorPalettes[self.CurrentPalette]
        groups = (not self.data[0].domain.classVar and 1) or len(self.data[0].domain.classVar.values) # mercy! (just had to do this)

        self.bmps = []; self.heights = []; self.imgStart = []; self.imgEnd = []
        for (i,hm) in enumerate(self.heatmaps):
            bmpl = []
            for g in range(groups):
                bmp, self.imageWidth, imageHeight = hm[g].getBitmap(int(self.CellWidth), int(self.CellHeight), lo, hi, self.Gamma)
                bmpl.append(bmp)
                if not i: self.heights.append(imageHeight)
            self.bmps.append(bmpl)

        self.canvas.resize(2000, 2000) # this needs adjustment
        x = c_offsetX; y0 = c_offsetY

        self.legend = self.heatmapconstructor[0].getLegend(self.imageWidth, c_legendHeight, self.Gamma)
        if self.LegendOnTop:
            y0 = self.drawLegend(x, y0, self.imageWidth, c_legendHeight, palette)

        for i in range(len(self.data)):
            y = y0; y1 = y0
            if self.ShowDataFileNames and len(self.data)>1:
                y1 = self.drawFileName(self.data[i].name, x, y, self.imageWidth+self.ShowAverageStripe*(c_averageStripeWidth + c_spaceX))
            x0 = x
            # plot the heatmap (and group label)
            ycoord = []
            y = y1; x += self.ShowAverageStripe * (c_averageStripeWidth + c_spaceX)
            for g in range(groups):
              if self.heights[g]:
                if self.ShowGroupLabel and groups>1:
                    y = self.drawGroupLabel(self.data[i][0].domain.classVar.values[g], x, y, self.imageWidth)
                if not i: self.imgStart.append(y)
                ycoord.append(y)
                image = ImageItem(self.bmps[i][g], self.canvas, self.imageWidth, self.heights[g], palette, x=x, y=y, z=z_heatmap)
                image.hm = self.heatmaps[i][g] # needed for event handling
                image.height = self.heights[g]; image.width = self.imageWidth
                if not i: self.imgEnd.append(y+self.heights[g]-1)
                y += self.heights[g] + c_spaceY

            x = x0
            # plot stripe with averages
            if self.ShowAverageStripe:
                for g in range(groups):
                    avg, avgWidth, avgHeight = self.heatmaps[i][g].getAverages(c_averageStripeWidth, int(self.CellHeight), lo, hi, self.Gamma)
                    ImageItem(avg, self.canvas, avgWidth, avgHeight, palette, x=x, y=ycoord[g])
            x += self.imageWidth + c_offsetX + self.ShowAverageStripe * (c_averageStripeWidth + c_spaceX)

        # plot the gene annotation
        for g in range(groups):
            if self.ShowGeneAnnotations and self.CellHeight>4:
                self.drawGeneAnnotation(x, ycoord[g], g)

        self.selection.redraw()        
        self.canvas.update()
        
    def createHeatMap(self):
        if len(self.data):
            merge = min(self.Merge, float(len(self.data[0])))
            squeeze = 1. / merge
            self.lowerBound = 1000; self.upperBound = -1000 # CHANGE!!!
            self.heatmaps = []
            for (i, hmc) in enumerate(self.heatmapconstructor):
                hm, lb, ub = hmc(squeeze)
                self.heatmaps.append(hm)
                self.lowerBound = min(self.lowerBound, lb)
                self.upperBound = max(self.upperBound, ub)

            self.sliderCutLow.setRange(self.lowerBound, 0, 0.1)
            self.sliderCutHigh.setRange(1e-10, self.upperBound, 0.1)
            self.CutLow = max(self.CutLow, self.lowerBound)
            self.CutHigh = min(self.CutHigh, self.upperBound)
            self.sliderCutLow.setValue(self.CutLow)
            self.sliderCutHigh.setValue(self.CutHigh)
            self.selection.remove()
        self.drawHeatMap()

    ##########################################################################
    # file list management

    # rel = -1 for up, or 1 for down
    def fileOrderChange(self, rel):
        sel = self.selectedFile
        data = self.data
        data[sel], data[sel+rel] = (data[sel+rel], data[sel])
        if sel == self.refFile:
            self.refFile += rel
        elif sel + rel == self.refFile:
            self.refFile += -rel
        # i got lazy here
        self.fileLB.clear()
        for i in range(len(data)):
            self.fileLB.insertItem(self.createListItem(data[i].name, i))
        self.fileLB.setSelected(sel + rel, 1)
        self.constructHeatmap()

    def setFileReferenceBySelection(self, sel):
        self.fileSelectionChanged(sel)
        self.setFileReference()
        
    def setFileReference(self):
        sel = self.selectedFile
        self.fileLB.changeItem(self.createListItem(self.data[self.refFile].name, -1), self.refFile)
        self.refFile = sel
        self.fileLB.changeItem(self.createListItem(self.data[sel].name, sel), sel)
        self.constructHeatmap()
##        self.canvas.update()

    def fileSelectionChanged(self, sel):
        # self.fileRef.setDisabled(sel==0)
        self.selectedFile = sel
        self.fileDown.setEnabled(sel < len(self.data)-1)
        self.fileUp.setEnabled(sel>0)
        self.fileRef.setEnabled(sel <> self.refFile)

    def createListItem(self, text, position):
        pixmap = QPixmap()
        pixmap.resize(14,13)
        pixmap.fill(Qt.white)
        
        if position == self.refFile:
            painter = QPainter()
            painter.begin(pixmap)
            painter.setPen(Qt.black)
            painter.setBrush(Qt.black)
            painter.drawRect(3, 3, 8, 8)
            painter.end()
            
        listItem = QListBoxPixmap(pixmap)
        listItem.setText(text)
        return listItem
##        return QListBoxText(text)


##################################################################################################
##################################################################################################
# color palette dialog

class OWDisplayColorOptions(OWOptions):
    def __init__(self, parent=None, name=None, master=None):
        OWOptions.__init__(self, "Color Palette", "OrangeWidgetsIcon.png", parent, name)
        self.master = master
        pms = [master.createColorStripe(i) for i in range(len(master.ColorPalettes))]
        OWGUI.radioButtonsInBox(self.top, master, "Color Palettes", pms, "CurrentPalette", tooltips=None, callback=self.master.setCurrentPalette)

##################################################################################################
# new canvas items

class ImageItem(QCanvasRectangle):
    def __init__(self, bitmap, canvas, width, height, palette, depth=8, numColors=256, x=0, y=0, z=0):
        QCanvasRectangle.__init__(self, canvas)
        self.image = QImage(bitmap, width, height, depth, palette, numColors, QImage.LittleEndian)
        self.image.bitmap = bitmap # this is tricky: bitmap should not be freed, else we get mess. hence, we store it in the object
        #self.pixmap = QPixmap()
        #self.pixmap.convertFromImage(image, QPixmap.Color)
        self.canvas = canvas
        self.setSize(width, height)
##        self.setSize(int(image.width()), int(image.height()))
        self.setX(x); self.setY(y); self.setZ(z)
        self.show()

    def drawShape(self, painter):
        painter.drawImage(self.x(), self.y(), self.image, 0, 0, -1, -1)
        
##################################################################################################
# mouse event handler

v_sel_width = 2

class MyCanvasView(QCanvasView):
    def __init__(self, *args):
        apply(QCanvasView.__init__,(self,) + args)
        self.canvas = args[0]
        self.clicked = False
        self.viewport().setMouseTracking(True)
        self.setFocusPolicy(QWidget.ClickFocus)
        self.setFocus()
        self.shiftPressed = False

    def heatmapParameters(self, master, cellWidth, cellHeight):
        self.master = master
        self.dx, self.dy = cellWidth, cellHeight
        self.selector = QCanvasRectangle(0, 0, self.dx + 2 * v_sel_width - 1, self.dy + 2 * v_sel_width - 1, self.canvas)
#        self.selector.setPen(QPen(Qt.black, v_sel_width))
        self.selector.setPen(QPen(self.master.SelectionColors[self.master.CurrentPalette], v_sel_width))
        self.selector.setZ(10)
        self.bubble = BubbleInfo(self.canvas)

    def contentsMouseMoveEvent(self, event):
        # handling of selection
        if self.clicked:
            self.master.selection(self.clicked, event.y())

        # balloon handling
        try:
            if self.master <> None and not self.master.BShowballoon: return
        except:
            return
        items = filter(lambda ci: ci.z()==z_heatmap, self.canvas.collisions(event.pos()))
        if len(items) == 0: # mouse over nothing special
            self.selector.hide()
            self.bubble.hide()
            self.canvas.update()
        else:
            item = items[0]
            hm = item.hm
            x, y = event.x() - item.x(), event.y() - item.y()
            if x<0 or y<0 or x>item.width-1 or y>item.height-1: 
                self.selector.hide()
                self.canvas.update()
                return
            col, row = int(x / self.dx), int(y / self.dy)
            # hm.getCellIntensity(row, col), hm.getRowIntensity(row)
            ex = hm.examples[hm.exampleIndices[row] : hm.exampleIndices[row+1]]
            self.selector.setX(item.x()+col*self.dx-v_sel_width+1)
            self.selector.setY(item.y()+row*self.dy-v_sel_width+1)
            self.selector.show()

            # bubble, construct head
            if hm.getCellIntensity(row, col):
                head = "%6.4f" % hm.getCellIntensity(row, col)
            else:
                head = "Missing Data"
            if self.master.BShowColumnID:
                head += "\n"+ex[0].domain.attributes[col].name
            self.bubble.head.setText(head)
            # bubble, construct body
            body = None
            if (self.master.BShowSpotIndex and self.master.BSpotVar) or self.master.BShowAnnotation or self.master.BShowGeneExpression:
                for (i, e) in enumerate(ex):
                    if i>5:
                        body += "\n... (%d more)" % (len(ex)-5)
                        break
                    else:
                        s = []
                        if self.master.BShowSpotIndex and self.master.BSpotVar:
                            s.append(str(e[self.master.BSpotVar]))
                        if self.master.BShowGeneExpression:
                            s.append(str(e[col]))
                        if self.master.BShowAnnotation and self.master.BAnnotationVar:
                            s.append(str(e[self.master.BAnnotationVar]))
                    if body: body += "\n"
                    else: body=""
                    body += reduce(lambda x,y: x + ' | ' + y, s)
            self.bubble.body.setText(body)
            self.bubble.show()
            self.bubble.move(event.x()+20, event.y()+20)
            
            self.canvas.update()

    def keyPressEvent(self, e):
        self.shiftPressed = e.key() == 4128

    def keyReleaseEvent(self, e):        
        self.shiftPressed = False

    def contentsMousePressEvent(self, event):
        # self.viewport().setMouseTracking(False)
        self.clicked = event.y()
        self.master.selection.start(self.clicked, event.y(), self.shiftPressed)

    def contentsMouseReleaseEvent(self, event):
        self.clicked = False
        self.canvas.update()
        self.master.selection.release()

##################################################################################################
# data selection

class SelectData:
    def __init__(self, master, canvas):
        self.canvas = canvas
        self.master = master
        self.add = FALSE
        self.squares = []; self.rows = []    # cumulative, used for multiple selection

    # removes the selection and relate information
    def remove(self):
        for r in self.squares:
            r.setCanvas(None)
        self.squares = []; self.rows = []

    # starts the selection, called after the first click
    def start(self, p1, p2, add=False):
        if not add:
            self.remove()
        self.cSquares = []; self.cRows = [] # current selection
        self.add = add
        self.__call__(p1, p2)
 
    # called during dragging (extending the selection)
    def __call__(self, p1, p2):
        for r in self.cSquares:
            r.setCanvas(None)
        y1 = min(p1, p2); y2 = max(p1, p2)
        self.cRows = self.findSelectionRows([(y1,y2)])
        self.cSquares = self.draw(self.cRows)

    # merges two selection lists, account for any overlaping or directly neighboring selections
    def mergeRows(self, pl1, pl2):
        new = []
        l = {}
        for (g, p1, p2) in pl1+pl2:
            if l.has_key((g,p1)): l[(g,p1)] += 1
            else: l[(g,p1)] = 1
            if l.has_key((g,p2)): l[(g,p2)] += -1
            else: l[(g,p2)] = -1
        kk = l.keys()
        kk.sort()
        current = 0 # sum of ups and downs, if positive then selection, else, empty
        intermediate = FALSE; # true if previous end was only one point before, indicates join
        for k in kk:
            if current == 0 and not intermediate:
                start = k
            current += l[k]
            if current == 0:
                if l.has_key((k[0],k[1]+1)):
                    intermediate = TRUE
                else:
                    intermediate = FALSE
                    new += [(start[0], start[1], k[1])]
        return new

    # mouse release, end of selection
    # if shift click, then merge with previous, else only remember current
    def release(self):
        if self.add:
            newrows = self.mergeRows(self.rows, self.cRows)
            self.remove()
            for r in self.cSquares:
                r.setCanvas(None)
            self.rows = newrows
            squares = self.draw(self.rows)
            self.squares = squares
        else:
            self.rows = self.cRows
            self.squares = self.cSquares
        if self.rows:
            self.master.sendData(self.rows)

    def findSelectionRows(self, points):
        start = self.master.imgStart; end = self.master.imgEnd
        rows = []
        for (y1, y2) in points: # this could be optimized, since points are ordered
            for i in range(len(start)):
                if y2<start[i]:
                    break
                if y1<end[i]:
                    if y1>start[i]:
                        a = int((y1-start[i])/self.master.CellHeight)
                    else:
                        a = 0
                    if y2<end[i]:
                        b = int((y2-start[i])/self.master.CellHeight)
                    else:
                        b = int((end[i]-start[i])/self.master.CellHeight)
                    rows.append((i,a,b))
        return rows

    # draws rectangles around selected points
    def draw(self, rows):
        start = self.master.imgStart; end = self.master.imgEnd
        x = c_offsetX + self.master.ShowAverageStripe * (c_averageStripeWidth + c_spaceX)
        lines = []
        for (g, r1, r2) in rows:
            y1 = start[g] + r1 * self.master.CellHeight
            y2 = start[g] + (r2+1) * self.master.CellHeight - 1
            r = QCanvasRectangle(x-v_sel_width+1, y1-v_sel_width+1, self.master.imageWidth+2*v_sel_width-1, y2-y1+v_sel_width+2, self.canvas)
            r.setPen(QPen(self.master.SelectionColors[self.master.CurrentPalette], v_sel_width))
            #r.setPen(QPen(Qt.red, v_sel_width))
            r.setZ(10)
            r.show()
            lines.append(r)
        self.canvas.update()
        return lines

    def redraw(self):
        for r in self.squares:
            r.setCanvas(None)
        if self.rows:
            self.squares = self.draw(self.rows)

##################################################################################################
# bubble info class

bubbleBorder = 4

class BubbleInfo(QCanvasRectangle):
    def __init__(self, *args):
        apply(QCanvasRectangle.__init__, (self,) + args)
        self.canvas = args[0]
        self.setBrush(QBrush(Qt.white))
        #self.setPen(QPen(Qt.black, v_sel_width))
        self.bubbleShadow = QCanvasRectangle(self.canvas)
        self.bubbleShadow.setBrush(QBrush(Qt.black))
        self.bubbleShadow.setPen(QPen(Qt.black))
        self.head = QCanvasText(self.canvas)
        self.line = QCanvasLine(self.canvas)
        self.body = QCanvasText(self.canvas)
        self.items = [self.head, self.line, self.body]
        self.setZ(110)
        self.bubbleShadow.setZ(109)
        for i in self.items:
            i.setZ(111)

    def move(self, x, y):
        QCanvasRectangle.move(self, x, y)
        self.setX(x); self.setY(y)
        self.bubbleShadow.move(x+5, y+5)
        for item in self.items:
            item.setX(x + bubbleBorder)
        w = max(100, self.head.boundingRect().width() + 2 * bubbleBorder, self.body.boundingRect().width() + 2 * bubbleBorder)
        y += 2
        self.head.setY(y)
        y += self.head.boundingRect().height()
        self.line.setPoints(0,0,w,0)
        self.line.setX(x); self.line.setY(y)
        y += 2
        self.body.setY(y)
        h = 2 * (2 + (self.body.text()<>None)) + self.head.boundingRect().height() + (self.body.text()<>None) * self.body.boundingRect().height()
        self.setSize(w,h)
        self.bubbleShadow.setSize(w,h)
        
    def show(self):
        QCanvasRectangle.show(self)
        self.bubbleShadow.show()
        self.head.show()
        if self.body.text():
            self.line.show()
            self.body.show()

    def hide(self):
        QCanvasRectangle.hide(self)
        self.bubbleShadow.hide()
        for item in self.items:
            item.hide()

##################################################################################################
# test script

if __name__=="__main__":
    import orange
    a = QApplication(sys.argv)
    ow = OWHeatMap()
    a.setMainWidget(ow)

    ow.show()
##    d = orange.ExampleTable('wt'); d.name = 'wt'
##    d = orange.ExampleTable('wt-nometa'); d.name = 'wt'
##    ow.dataset(d, 2)
##    d = orange.ExampleTable('wt-nometa'); d.name = 'wt'
##    ow.dataset(d, 1)
##    ow.dataset(None, 1)
##    ow.dataset(None, 2)

    names = ['wt1', 'wt2', 'wt3', 'wt4']
    for i, s in enumerate(names): 
        d = orange.ExampleTable(s); d.name = s
        ow.dataset(d, i)
    ow.dataset(None, 2)

    a.exec_loop()
    ow.saveSettings()
