"""
<name>Heat Map</name>
<description>Microarray Heat Map</description>
<category>Genomics</category>
<icon>icons/ClassificationTreeViewer2D.png</icon>
<priority>1010</priority>
"""

import orange, math
import OWGUI, OData
#from string import *
from qt import *
from qtcanvas import *
from OWWidget import *
from OWOptions import *
from qwt import *

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
                    "BSpotVar", "BAnnotationVar",
                    "CurrentPalette"]

    def __init__(self, parent=None, name='OWHeatMap'):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, name, 'Microarray Heat Map', FALSE, FALSE) 
        
        self.inputs = [("Examples", ExampleTableWithClass, self.microarrayData, 1)]
        self.outputs = [("Examples", ExampleTable)]

        #set default settings
        self.CellWidth = 3; self.CellHeight = 3
        self.Merge = 1; self.savedMerge = self.Merge
        self.Gamma = 1
        self.CutLow = 0; self.CutHigh = 0; self.CutEnabled = 1
        self.ShowAnnotation = 0
        self.LegendOnTop = 1           # legend stripe on top (bottom)?
        self.LegendOnBottom = 1
        self.ShowGroupLabel = 1        # show class names in case of classified data?
        self.ShowAverageStripe = 1     # show the stripe with the evarage
        self.MaintainArrayHeight = 1   # adjust cell height while changing the merge factor
        self.BShowballoon = 1           # balloon help
        self.BShowColumnID = 1; self.BShowSpotIndex = 1; self.BShowAnnotation = 1; self.BShowGeneExpression = 1
        self.BSpotVar = None; self.BAnnotationVar = None
        self.setColorPalette()
        
        self.loadSettings()
        self.data = None
        self.maxHSize = 10; self.maxVSize = 10

        # GUI definition
        self.tabs = QTabWidget(self.controlArea, 'tabWidget')

        # SETTINGS TAB
        settingsTab = QVGroupBox(self)
        box = QVButtonGroup("Cell Size (Pixels)", settingsTab)
        OWGUI.qwtHSlider(box, self, "CellWidth", label='Width: ', labelWidth=38, minValue=1, maxValue=self.maxHSize, step=1, precision=0, callback=self.createHeatMap)
        self.sliderVSize = OWGUI.qwtHSlider(box, self, "CellHeight", label='Height: ', labelWidth=38, minValue=1, maxValue=self.maxVSize, step=1, precision=0, callback=self.createHeatMap)
        OWGUI.qwtHSlider(settingsTab, self, "Gamma", box="Gamma", minValue=0.1, maxValue=1, step=0.1, callback=self.createHeatMap)

        # define the color stripe to show the current palette
        colorItems = [self.createColorStripe(i) for i in range(len(self.ColorPalettes))] + ["Custom ..."]
        OWGUI.comboBox(settingsTab, self, "CurrentPalette", box="Colors", items=colorItems, tooltip=None, callback=self.setColor)

        box = QVButtonGroup("Annotations", settingsTab)
        OWGUI.checkOnly(box, self, 'Legend (top)', 'LegendOnTop', callback=self.createHeatMap)
        OWGUI.checkOnly(box, self, 'Stripes with averages', 'ShowAverageStripe', callback=self.createHeatMap)
        
        self.tabs.insertTab(settingsTab, "Settings")

        # FILTER TAB
        tab = QVGroupBox(self)
        box = QVButtonGroup("Treshold Values", tab)
        OWGUI.checkOnly(box, self, "Enabled", 'CutEnabled', callback=self.setCutEnabled)
        self.sliderCutLow = OWGUI.qwtHSlider(box, self, 'CutLow', label='Low:', labelWidth=33, minValue=-100, maxValue=0, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.createHeatMap)
        self.sliderCutHigh = OWGUI.qwtHSlider(box, self, 'CutHigh', label='High:', labelWidth=33, minValue=0, maxValue=100, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.createHeatMap)
        if not self.CutEnabled:
            self.sliderCutLow.box.setDisabled(1)
            self.sliderCutHigh.box.setDisabled(1)

        box = QVButtonGroup("Merge", tab)
        OWGUI.qwtHSlider(box, self, "Merge", label='Rows:', labelWidth=33, minValue=1, maxValue=100, step=1, callback=self.mergeChanged, ticks=0)
        OWGUI.checkOnly(box, self, "Maintain array height", 'MaintainArrayHeight')
        
            
        self.tabs.insertTab(tab, "Filter")

        # INFO TAB
        tab = QVGroupBox(self)
        box = QVButtonGroup("Balloon", tab)
        OWGUI.checkOnly(box, self, "Show balloon", 'BShowballoon', callback=lambda: self.balloonInfoBox.setDisabled(not self.BShowballoon))
        box = QVButtonGroup("Balloon Info", tab)
        OWGUI.checkOnly(box, self, "Column ID", 'BShowColumnID')
        OWGUI.checkOnly(box, self, "Spot Index", 'BShowSpotIndex', callback=lambda: self.spotCombo.setDisabled(not self.BShowSpotIndex))
        self.spotCombo = OWGUI.comboBox(box, self, "BSpotVar", items=[])
        OWGUI.checkOnly(box, self, "Annotation", 'BShowAnnotation', callback=lambda: self.annotationCombo.setDisabled(not self.BShowAnnotation))
        self.annotationCombo = OWGUI.comboBox(box, self, "BAnnotationVar", items=[])
        OWGUI.checkOnly(box, self, "Gene expression", 'BShowGeneExpression')
        self.balloonInfoBox = box
        self.tabs.insertTab(tab, "Info")
        
        self.resize(400,400)

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
        self.ColorPalettes = \
          ([qRgb(255.*i/250., 255.*i/250., 255-(255.*i/250.)) for i in range(250)] + [white]*6,
           [qRgb(0, 255.*i*2/250., 0) for i in range(125, 0, -1)] + [qRgb(255.*i*2/250., 0, 0) for i in range(125)] + [white]*6,
           [qRgb(255.*i/250., 0, 0) for i in range(250)] + [white]*6)
        self.CurrentPalette = 0
        
    # any time the data changes, the two combo boxes showing meta attributes
    # have to be adjusted
    def setMetaCombo(self, cb, value, enabled):
        cb.clear()
        self.meta = [m.name for m in self.data.domain.getmetas().values()]
        if len(self.meta)==0:
            cb.setDisabled(True)
            return
        for m in self.meta:
            cb.insertItem(m)
        if value in self.meta:
            cb.setCurrentItem(self.meta.index(value))
        else:
            cb.setCurrentItem(0)
        cb.setDisabled(not enabled)

    def setMetaCombos(self):
        self.setMetaCombo(self.spotCombo, 'RMI', self.BShowSpotIndex)
        self.setMetaCombo(self.annotationCombo, 'annotation', self.BShowAnnotation)

    ##########################################################################
    # handling of input signals

    def microarrayData(self, data):
        self.data = data

        # figure out the max in min value of expression
        bs = orange.DomainBasicAttrStat(self.data)
        minVal = min([bs[x].min for x in data.domain.attributes])
        maxVal = max([bs[x].max for x in data.domain.attributes])
        #self.sliderCutLow.setScale(minVal, 0)
        self.sliderCutLow.setRange(minVal, 0, 0.1)
        self.sliderCutHigh.setRange(0, maxVal, 0.1)
        changed = 0
        if minVal > self.CutLow:
            self.CutLow = minValue
            changed = 1
        if maxVal < self.CutLow:
            self.CutLow = maxValue
            changed = 1
        self.sliderCutLow.setValue(self.CutLow)
        self.sliderCutHigh.setValue(self.CutHigh)
        
        self.heatmapconstructor = orange.HeatmapConstructor(self.data)
        self.createHeatMap()
        self.setMetaCombos() # set the two combo widgets according to the data
        
    ##########################################################################
    # callback functions

    def setColor(self):
        if self.CurrentPalette == len(self.ColorPalettes):
            self.CurrentPalette = 0
            # put a code here that allows to define ones own colors
        else:
            pm = self.createColorStripe(self.CurrentPalette)
            self.createHeatMap()

    def setCutEnabled(self):
        self.sliderCutLow.box.setDisabled(not self.CutEnabled)
        self.sliderCutHigh.box.setDisabled(not self.CutEnabled)
        self.createHeatMap()

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

    ##########################################################################
    # drawing

    def drawLegend(self, x, y, width, height, palette):
        legend = self.heatmapconstructor.getLegend(width, height, self.Gamma)

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

    def drawGroupLabel(self, label, x, y, width):
        t = QCanvasText(label, self.canvas)
        t.setX(x); t.setY(y)
        t.show()
        return y + t.boundingRect().height() + 1

    def drawHeatMap(self):
        lo = self.CutEnabled and self.CutLow   or self.lowerBound
        hi = self.CutEnabled and self.CutHigh  or self.upperBound

        # remove everything from current canvas
        for i in self.canvas.allItems():
            i.setCanvas(None)
        self.canvasView.heatmapParameters(self, self.CellWidth, self.CellHeight) # needed for event handling

        palette = self.ColorPalettes[self.CurrentPalette]
        groups = (not data.domain.classVar and 1) or len(data.domain.classVar.values) # mercy! (just had to do this)

        self.bmps = []; self.heights = []; self.imgStart = []; self.imgEnd = []
        for g in range(groups):
            bmp, self.imageWidth, imageHeight = self.heatmaps[g].getBitmap(int(self.CellWidth), int(self.CellHeight), lo, hi, self.Gamma)
            self.bmps.append(bmp)
            self.heights.append(imageHeight)

        self.canvas.resize(2000, 2000) # this needs adjustment

        x = c_offsetX; y = c_offsetY
        if self.ShowAverageStripe:
            x += c_averageStripeWidth + c_spaceX

        self.legend = self.heatmapconstructor.getLegend(self.imageWidth, c_legendHeight, self.Gamma)
        if self.LegendOnTop:
            y = self.drawLegend(x, y, self.imageWidth, c_legendHeight, palette)

        for g in range(groups):
            if self.ShowGroupLabel and groups>1:
                y = self.drawGroupLabel(data.domain.classVar.values[g], x, y, self.imageWidth)
            self.imgStart.append(y)
            image = ImageItem(self.bmps[g], self.canvas, self.imageWidth, self.heights[g], palette, x=x, y=y, z=z_heatmap)
            image.hm = self.heatmaps[g] # needed for event handling
            image.height = self.heights[g]; image.width = self.imageWidth
            if self.ShowAverageStripe:
                avg, avgWidth, avgHeight = self.heatmaps[g].getAverages(c_averageStripeWidth, int(self.CellHeight), lo, hi, self.Gamma)
                ImageItem(avg, self.canvas, avgWidth, avgHeight, palette, x=c_offsetX, y=y)
            self.imgEnd.append(y+self.heights[g]-1)
            y += self.heights[g] + c_spaceY

        print self.imgEnd
        for i in range(len(self.imgStart)):
            l = QCanvasLine(self.canvas)
            l.setPoints(2, self.imgStart[i], 2, self.imgEnd[i])
            print 'line %d->%d ' % (self.imgStart[i], self.imgEnd[i]),
            l.show()
        print
        self.canvas.update()
        
    def createHeatMap(self):
        merge = min(self.Merge, float(len(self.data)))
        squeeze = 1. / merge
        self.heatmaps, self.lowerBound, self.upperBound = self.heatmapconstructor(squeeze)
        self.drawHeatMap()

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
        self.selector.setPen(QPen(Qt.black, v_sel_width))
        self.selector.setZ(10)
        self.bubble = BubbleInfo(self.canvas)

    def contentsMouseMoveEvent(self, event):
        # handling of selection
        if self.clicked:
            self.master.selection(self.clicked, event.y(), self.shiftPressed)


        # balloon handling
        if not self.master.BShowballoon: return
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
            head = "%6.4f" % hm.getCellIntensity(row, col)
            if self.master.BShowColumnID:
                head += "\n"+ex[0].domain.attributes[col].name
            self.bubble.head.setText(head)
            # bubble, construct body
            body = None
            if self.master.BShowSpotIndex or self.master.BShowAnnotation or self.master.BShowGeneExpression:
                for (i, e) in enumerate(ex):
                    if i>5:
                        body += "\n..."
                        break
                    else:
                        s = []
                        if self.master.BShowSpotIndex:
                            s.append(str(e[self.master.meta[self.master.BSpotVar]]))
                        if self.master.BShowGeneExpression:
                            s.append(str(e[col]))
                        if self.master.BShowAnnotation:
                            s.append(str(e[self.master.meta[self.master.BAnnotationVar]]))
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
        self.master.selection(self.clicked, event.y(), self.shiftPressed)

    def contentsMouseReleaseEvent(self, event):
        self.clicked = False
        self.canvas.update()
        self.master.selection.release()

##################################################################################################
# data selection

class SelectData:
    def __init__(self, master, canvas):
        self.canvas = canvas
        self.currentLines = None
        self.master = master
        self.squares = []
        self.points = []

    def __call__(self, p1, p2, add=False):
        if not add and self.squares:
            for r in self.squares:
                r.setCanvas(None)
            self.squares = []
            self.points = []
        points = self.findSelectionPoints(p1, p2)
        self.points = self.mergePoints(self.points, points)
        self.squares = self.drawOne(points)

    def mergePoints(self, pl1, pl2):
        print '->', pl1, '|||', pl2
        
        return pl1

    def release(self):
        pass
        
    def findSelectionPoints(self, p1, p2):
        y1 = min(p1, p2); y2 = max(p1, p2)
        start = self.master.imgStart; end = self.master.imgEnd
        s = []
        for i in range(len(start)):
            if y2<start[i]:
                break
            if y1<end[i]:
                if y1>start[i]:
                    a = start[i] + int((y1-start[i])/self.master.CellHeight) * self.master.CellHeight
                else:
                    a = start[i]
                if y2<end[i]:
                    b = start[i] + (int((y2-start[i])/self.master.CellHeight)+1) * self.master.CellHeight-1
                else:
                    b = end[i]
                s.append((a,b))
        return s

    def drawOne(self, points):
        x = c_offsetX + self.master.ShowAverageStripe * (c_averageStripeWidth + c_spaceX)
        lines = []
        for (y1,y2) in points:
            r = QCanvasRectangle(x-v_sel_width+1, y1-v_sel_width+1, self.master.imageWidth+2*v_sel_width-1, y2-y1+v_sel_width+2, self.canvas)
            r.setPen(QPen(Qt.black, v_sel_width))
            r.setZ(10)
            r.show()
            lines.append(r)
        self.canvas.update()
        return lines

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

##    data = orange.ExampleTable('wt-large')
##    data = orange.ExampleTable('wt')
    data = orange.ExampleTable('wtclassed')
    ow.microarrayData(data)
    ow.show()
    a.exec_loop()
    ow.saveSettings()
