"""
<name>Heat Map</name>
<description>Heatmap-based data visualization.</description>
<contact>Ales Erjavec, Blaz Zupan, Janez Demsar</contact>
<icon>icons/HeatMap.png</icon>
<priority>50</priority>
"""

import orange, math
import orangene
import OWGUI
from OWWidget import *
from OWDlgs import OWChooseImageSizeDlg
from ColorPalette import signedPalette
from OWClustering import HierarchicalClusterItem
import OWColorPalette
try:
    from OWDataFiles import DataFiles
except Exception:
    class DataFiles(object):
        pass

import warnings
warnings.filterwarnings("ignore", "'strain'", orange.AttributeWarning)

# from OWChipANOVA import ANOVAResults

##############################################################################
# parameters that determine the canvas layout

c_offsetX = 10; c_offsetY = 10  # top and left border
c_spaceY = 10                   # space btw graphical elements
c_spaceAverageX = 5             # space btw stripe with average and microarray
c_legendHeight = 15             # height of the legend
c_averageStripeWidth = 12       # width of the stripe with averages

z_heatmap = 5                   # layer with heatmaps

##############################################################################
# main class

class OWHeatMap(OWWidget):	
    settingsList = ["CellWidth", "CellHeight", "SpaceX", "Merge",
                    "Gamma", "CutLow", "CutHigh", "CutEnabled", 
                    "ShowAnnotation", "LegendOnTop", "LegendOnBottom",
                    "ShowAverageStripe", "ShowGroupLabel",
                    "MaintainArrayHeight",
                    "BShowballoon", "BShowColumnID", "BShowSpotIndex",
                    "BShowAnnotation", 'BShowGeneExpression',
                    "BSpotVar", "ShowGeneAnnotations",
                    "ShowDataFileNames", "BAnnotationVar",
                    "SelectionType",
                    "CurrentPalette", "SortGenes", "colorSettings", "selectedSchemaIndex",
                    "palette", "ShowColumnLabels", "ColumnLabelPosition"]

    def __init__(self, parent=None, signalManager = None):
        self.callbackDeposit = [] # deposit for OWGUI callback functions
        OWWidget.__init__(self, parent, signalManager, 'HeatMap', TRUE)
        
        self.inputs = [("Structured Data", DataFiles, self.chipdata, Single + NonDefault), ("Examples", ExampleTable, self.dataset, Default + Multiple)]
        self.outputs = [("Structured Data", DataFiles, Single + NonDefault), ("Examples", ExampleTable, Default)]

        #set default settings
        self.CellWidth = 3; self.CellHeight = 3
        self.SpaceX = 10
        self.Merge = 1; self.savedMerge = self.Merge
        self.Gamma = 1
        self.CutLow = 0; self.CutHigh = 0; self.CutEnabled = 0
        self.ShowAnnotation = 0
        self.LegendOnTop = 0           # legend stripe on top (bottom)?
        self.LegendOnBottom = 1
        self.ShowGroupLabel = 1        # show class names in case of classified data?
        self.ShowAverageStripe = 0     # show the stripe with the evarage
        self.MaintainArrayHeight = 0   # adjust cell height while changing the merge factor
        self.BShowballoon = 1          # balloon help
        self.ShowGeneAnnotations = 1   # show annotations for genes
        self.ShowColumnLabels = 1
        self.ColumnLabelPosition = 0
        self.ShowDataFileNames = 1     # show the names of the data sets (only in case of multiple files)
        self.BShowColumnID = 1; self.BShowSpotIndex = 1; self.BShowAnnotation = 1; self.BShowGeneExpression = 1
        self.BSpotVar = None; self.BAnnotationVar = None  # these are names of variables
        self.BSpotIndx = None; self.BAnnotationIndx = None # these are id's of the combo boxes
        self.SortGenes = 1
        self.ShowClustering = 1
        self.SelectionType = 0         # selection on a single data set
        self.setColorPalette()
        self.refFile = 0               # position index of a reference file
        self.selectedFile = None       # position of the selected file in the list box

        self.colorSettings =None
        self.selectedSchemaIndex = 0

        self.palette = self.ColorPalettes[0]
        
        self.loadSettings()
        self.data = []
        self.maxHSize = 30; self.maxVSize = 15


        # GUI definition
        self.connect(self.graphButton, SIGNAL("clicked()"), self.saveFig)
        self.tabs = OWGUI.tabWidget(self.controlArea) #QTabWidget(self.controlArea, 'tabWidget')

        # SETTINGS TAB
        settingsTab = OWGUI.createTabPage(self.tabs, "Settings") #QVGroupBox(self)
        box = OWGUI.widgetBox(settingsTab, "Cell Size (Pixels)", addSpace=True) #QVButtonGroup("Cell Size (Pixels)", settingsTab)
        OWGUI.qwtHSlider(box, self, "CellWidth", label='Width: ', labelWidth=38, minValue=1, maxValue=self.maxHSize, step=1, precision=0, callback=self.drawHeatMap)
        self.sliderVSize = OWGUI.qwtHSlider(box, self, "CellHeight", label='Height: ', labelWidth=38, minValue=1, maxValue=self.maxVSize, step=1, precision=0, callback=self.createHeatMap)
        OWGUI.qwtHSlider(box, self, "SpaceX", label='Space: ', labelWidth=38, minValue=0, maxValue=50, step=2, precision=0, callback=self.drawHeatMap)
        OWGUI.qwtHSlider(settingsTab, self, "Gamma", box="Gamma", minValue=0.1, maxValue=1, step=0.1, callback=self.drawHeatMap)
        OWGUI.separator(settingsTab)

        # define the color stripe to show the current palette
        box = OWGUI.widgetBox(settingsTab, "Color", orientation = "horizontal")
        self.colorCombo = OWColorPalette.PaletteSelectorComboBox(self)
        try:
            self.colorCombo.setPalettes("palette", self.createColorDialog())
        except Exception, ex:
            print >> sys.stderr, ex, "Error loading saved color palettes!\nCreating new default palette!"
            self.colorSettings = None
            self.colorCombo.setPalettes("palette", self.createColorDialog())
        self.colorCombo.setCurrentIndex(self.selectedSchemaIndex)
        self.setColor(self.selectedSchemaIndex, update=False)
        self.connect(self.colorCombo, SIGNAL("activated(int)"), self.setColor)
        box.layout().addWidget(self.colorCombo, 2)
        button = OWGUI.button(box, self, "Edit colors", callback=self.openColorDialog, tooltip="Edit the heatmap color palette")
        
        OWGUI.separator(settingsTab)

##        OWGUI.checkBox(settingsTab, self, "SortGenes", "Sort genes", box="Sort", callback=self.constructHeatmap)
        OWGUI.comboBox(settingsTab, self, "SortGenes", "Sort genes", items=["No sorting", "Sort genes", "Clustering", "Clustering with leaf ordering"], callback=self.constructHeatmap)
        OWGUI.rubber(settingsTab)
        
        # FILTER TAB
        tab = OWGUI.createTabPage(self.tabs, "Filter") #QVGroupBox(self)
        box = OWGUI.widgetBox(tab, "Threshold Values", addSpace=True) #QVButtonGroup("Threshold Values", tab)
        OWGUI.checkBox(box, self, 'CutEnabled', "Enabled", callback=self.setCutEnabled)
        self.sliderCutLow = OWGUI.qwtHSlider(box, self, 'CutLow', label='Low:', labelWidth=33, minValue=-100, maxValue=0, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.drawHeatMap)
        self.sliderCutHigh = OWGUI.qwtHSlider(box, self, 'CutHigh', label='High:', labelWidth=33, minValue=0, maxValue=100, step=0.1, precision=1, ticks=0, maxWidth=80, callback=self.drawHeatMap)
        if not self.CutEnabled:
            self.sliderCutLow.box.setDisabled(1)
            self.sliderCutHigh.box.setDisabled(1)

        box = OWGUI.widgetBox(tab, "Merge", addSpace=True) #QVButtonGroup("Merge", tab)
##        OWGUI.qwtHSlider(box, self, "Merge", label='Rows:', labelWidth=33, minValue=1, maxValue=500, step=1, callback=self.mergeChanged, precision=0, ticks=0)
        OWGUI.spin(box, self, "Merge", min=1, max=500, step=1, label='Rows:', callback=self.mergeChanged, callbackOnReturn=True)
        OWGUI.checkBox(box, self, 'MaintainArrayHeight', "Maintain array height")
        OWGUI.rubber(tab)

        # INFO TAB
        tab = OWGUI.createTabPage(self.tabs, "Info") #QVGroupBox(self)

        box = OWGUI.widgetBox(tab,'Annotation && Legends') #QVButtonGroup("Annotation && Legends", tab)
        OWGUI.checkBox(box, self, 'LegendOnTop', 'Show legend', callback=self.drawHeatMap)
        OWGUI.checkBox(box, self, 'ShowAverageStripe', 'Stripes with averages', callback=self.drawHeatMap)
        self.geneAnnotationsCB = OWGUI.checkBox(box, self, 'ShowGeneAnnotations', 'Gene annotations', callback=self.drawHeatMap)
        
        self.annotationCombo = OWGUI.comboBox(box, self, "BAnnotationIndx", items=[], callback=lambda x='BAnnotationVar', y='BAnnotationIndx': self.setMetaID(x, y))

        box = OWGUI.widgetBox(tab, 'Column Labels')
        columnLabelCB = OWGUI.checkBox(box, self, "ShowColumnLabels", "Display column labels", callback=self.drawHeatMap)
        comboBox = OWGUI.comboBox(OWGUI.indentedBox(box), self, "ColumnLabelPosition", "Position", items=["Top", "Bottom"], callback=self.drawHeatMap)
        columnLabelCB.disables.append(comboBox.box)
        columnLabelCB.makeConsistent()
        
        box = OWGUI.widgetBox(tab, "Ballon") #QVButtonGroup("Balloon", tab)
        OWGUI.checkBox(box, self, 'BShowballoon', "Show balloon", \
            callback=lambda: self.balloonInfoBox.setDisabled(not self.BShowballoon))
        box = OWGUI.widgetBox(tab, "Ballon info") #QVButtonGroup("Balloon Info", tab)
        OWGUI.checkBox(box, self, 'BShowColumnID', "Column ID")
        self.spotIndxCB = OWGUI.checkBox(box, self, 'BShowSpotIndex', "Spot Index", \
            callback=lambda: self.spotCombo.setDisabled(not self.BShowSpotIndex))
        self.spotCombo = OWGUI.comboBox(box, self, "BSpotIndx", items=[], \
            callback=lambda x='BSpotVar', y='BSpotIndx': self.setMetaID(x, y))
        OWGUI.checkBox(box, self, 'BShowGeneExpression', "Gene expression")
        OWGUI.checkBox(box, self, 'BShowAnnotation', "Annotation")
        self.balloonInfoBox = box
        OWGUI.rubber(tab)

        # FILES TAB
        self.filesTab = OWGUI.createTabPage(self.tabs, "Files")
        box = OWGUI.widgetBox(self.filesTab, "Data Files")
        self.fileLB = QListWidget(box)
        box.layout().addWidget(self.fileLB)
##        self.fileLB.setMaximumWidth(10)
        self.connect(self.fileLB, SIGNAL("highlighted(int)"), self.fileSelectionChanged)
        self.connect(self.fileLB, SIGNAL("selected(int)"), self.setFileReferenceBySelection)
##        self.tabs.insertTab(self.filesTab, "Files")
        self.tabs.setTabEnabled(self.tabs.indexOf(self.filesTab), 0)
        hbox = OWGUI.widgetBox(box, orientation="horizontal") #QHBox(box)
        self.fileUp = OWGUI.button(hbox, self, 'Up', \
            callback=lambda i=-1: self.fileOrderChange(i), disabled=1)
        self.fileRef = OWGUI.button(hbox, self, 'Ref', self.setFileReference, disabled=1)
        self.fileDown = OWGUI.button(hbox, self, 'Down', \
            callback=lambda i=1: self.fileOrderChange(i), disabled=1)
        for btn in [self.fileUp, self.fileRef, self.fileDown]:
            btn.setMaximumWidth(45)

        OWGUI.checkBox(self.filesTab, self, 'ShowDataFileNames', 'Show data file names',
           callback=self.drawHeatMap)
        OWGUI.radioButtonsInBox(self.filesTab, self, 'SelectionType', \
           ['Single data set', 'Multiple data sets'], box='Selection', callback=self.removeSelection)
        OWGUI.rubber(self.filesTab)

        self.resize(700,400)

        # canvas with microarray
        self.scene = HeatMapGraphicsScene()
        self.sceneView = MyGraphicsView(self.scene, self.mainArea)
        self.selection = SelectData(self, self.scene)
        self.currentHighlightedCluster = None
        self.selectedClusters = []
        self.mainArea.layout().addWidget(self.sceneView)

    def createColorStripe(self, palette):
        dx = 104; dy = 18
        bmp = chr(252)*dx*2 + reduce(lambda x,y:x+y, \
           [chr(i*250/dx) for i in range(dx)] * (dy-4)) + chr(252)*dx*2 
##        image = QImage(bmp, dx, dy, 8, self.ColorPalettes[palette], 256, QImage.LittleEndian)
        image = QImage(bmp, dx, dy, QImage.Format_Indexed8)# self.ColorPalettes[palette], 256, QImage.LittleEndian)
        image.setColorTable(signedPalette(self.ColorPalettes[palette]))

        pm = QPixmap.fromImage(image, Qt.AutoColor);
        return pm

    # set the default palettes used in the program
    # palette defines 256 colors, 250 are used for heat map, remaining 6 are extra
    # color indices for unknown is 255, underflow 253, overflow 254, white 252
    def setColorPalette(self):
        white = qRgb(255,255,255)
        gray = qRgb(200,200,200)
        self.ColorPalettes = \
          ([qRgb(255.*i/250., 255.*i/250., 255-(255.*i/250.)) \
            for i in range(250)] + [white]*3 + [qRgb(0., 0., 255.), qRgb(255., 255., 0.), gray],
           [qRgb(0, 255.*i*2/250., 0) for i in range(125, 0, -1)] \
           + [qRgb(255.*i*2/250., 0, 0) for i in range(125)] + [white]*3 \
           + [qRgb(0, 255., 0), qRgb(255., 0, 0), gray],
           [qRgb(255.*i/250., 0, 0) for i in range(250)] + [white]*3 \
           + [qRgb(0., 0, 0), qRgb(255., 0, 0), gray])
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
            cb.addItem(m)
        
        if not (value in self.meta):
            if default in self.meta:
                value = default
            else:
                value = None

        if value in self.meta:
            cb.setCurrentIndex(self.meta.index(value))
            indx = self.meta.index(value)
        else:
            cb.setCurrentIndex(0)
            value = self.meta[0]; indx = 0
        return (value, indx)

    def setMetaID(self, val, valIndx):
        setattr(self, val, self.meta[getattr(self, valIndx)])
        if val=='BAnnotationVar':
            self.drawHeatMap()

    def setMetaCombos(self):
        self.meta = [m.name for m in self.data[0].domain.getmetas().values()]
        self.BSpotVar, self.BSpotIndx = self.setMetaCombo(self.spotCombo, self.BSpotVar, \
            enabled=self.BShowSpotIndex, default='RMI')
        self.BAnnotationVar, self.BAnnotationIndx = self.setMetaCombo(self.annotationCombo, \
            self.BAnnotationVar, enabled=self.BShowAnnotation, default='xannotation')

    def saveFig(self):
        sizeDlg = OWChooseImageSizeDlg(self.scene)
        sizeDlg.exec_()

    ##########################################################################
    # handling of input/output signals

    def dataset(self, data, id, blockUpdate=0):
        ids = [d.id for d in self.data]
        if not data:
            if id in ids:
                k = ids.index(id)
                del self.data[k]
                self.fileLB.takeItem(k)
                if self.refFile == k:
                    self.refFile = 0
                    if len(self.data):
                        self.fileLB.changeItem(self.createListItem(self.data[0].name, self.refFile), \
                                               self.refFile)
        else:
            # check if the same length
            if data.domain.classVar:
                domain = self.checkDomain(data)
                if domain:
                    data = orange.ExampleTable(domain, data)
            data.setattr("id", id)
            if id in ids:
                indx = ids.index(id)
                self.data[indx] = data
##                self.fileLB.changeItem(self.createListItem(data.name, indx), indx)
                self.fileLB.takeItem(indx)
                self.fileLB.insertItem(indx, self.createListItem(data.name, indx))
            else:
                self.fileLB.addItem(self.createListItem(data.name, len(self.data)))
                self.data.append(data)

            if len(self.data) > 1:
                self.tabs.setTabEnabled(self.tabs.indexOf(self.filesTab), True)
            else:
                self.tabs.setTabEnabled(self.tabs.indexOf(self.filesTab), False)
            self.setMetaCombos() # set the two combo widgets according to the data

        self.unorderedData = None
        self.groupClusters = None
        if not blockUpdate:
            self.send('Examples', None)
            self.send('Structured Data', None)
            self.constructHeatmap()
            self.scene.update()

    def chipdata(self, data):
        self.data = [] # XXX should only remove the data from the same source, use id in this rutine
        self.fileLB.clear()
        self.refFile = 0
        if not data:
            for i in self.scene.items():
                self.scene.removeItem(i)
            self.scene.update()
            return
        indx = 0
        for (strainname, ds) in data:
            for d in ds:
                self.dataset(d, indx, blockUpdate=1)
                indx += 1
#        self.createHeatMap()

        pb = OWGUI.ProgressBar(self, iterations=len(self.data))
        self.constructHeatmap(callback=pb.advance)
        self.scene.update()
        pb.finish()

    def orderClustering(self, data):
        import orngClustering
        self.progressBarInit()
        clusterRoots = []
        orderedData = []
        mapping = []
        progressCallback = lambda value, caller=None: self.progressBarSet(value)
        if data.domain.classVar and data.domain.classVar.values:
            valuesCount = len(data.domain.classVar.values)
            for i, val in enumerate(data.domain.classVar.values):
                self.progressBarSet(100.0*i/valuesCount)
                tmpData = orange.ExampleTable([ex for ex in data if ex.getclass()==val])
                root = orngClustering.hierarchicalClustering(tmpData, progressCallback=lambda value: progressCallback((100.0*i + value)/valuesCount), order=self.SortGenes==3)
                orderedData.extend([tmpData[i] for i in root.mapping])
                mapping.extend([i+len(mapping) for i in root.mapping])
                clusterRoots.append(root)
            
        else:
            root = orngClustering.hierarchicalClustering(data, progressCallback=progressCallback, order=self.SortGenes==3)
            orderedData.extend([data[i] for i in root.mapping])
            mapping = list(root.mapping)
            clusterRoots.append(root)

        self.progressBarFinished()
        return orange.ExampleTable(orderedData), clusterRoots, mapping
            
        
    def constructHeatmap(self, callback=None):
        if len(self.data):
            self.heatmapconstructor = [None] * len(self.data)
            self.unorderedData = self.data if not self.unorderedData else self.unorderedData
            self.groupClusters = []
            self.attrCluster = None
            self.mapping = None
            sortData = lambda data, mapping: orange.ExampleTable(data.domain, [data[i] for i in mapping])
            if self.SortGenes:
                if self.SortGenes > 1: ## cluster sort
                    refData, self.groupClusters , self.mapping = self.orderClustering(self.unorderedData[self.refFile])
                    sortedData = sortData(self.data[self.refFile], self.mapping)
                    self.heatmapconstructor[self.refFile] = \
                        orangene.HeatmapConstructor(sortedData, None)
                    self.heatmapconstructor[self.refFile].setattr("_sortedData", sortedData)
                else:
                    self.heatmapconstructor[self.refFile] = \
                        orangene.HeatmapConstructor(self.data[self.refFile])
            else:
                self.heatmapconstructor[self.refFile] = \
                    orangene.HeatmapConstructor(self.data[self.refFile], None)
            if callback: callback()

            for i in range(len(self.data)):
                if i <> self.refFile:
                    if self.mapping:
                        self.heatmapconstructor[i] = orangene.HeatmapConstructor(sortData(self.data[i],self.mapping),
                            self.heatmapconstructor[self.refFile])
                    else:                        
                        self.heatmapconstructor[i] = orangene.HeatmapConstructor(self.data[i],
                            self.heatmapconstructor[self.refFile])
                    if callback: callback()
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
    def prepareData(self, rows, indx):
        ex = []
        for (g,s,e) in rows:
            hm = self.heatmaps[indx][g]
            ex += hm.examples[hm.exampleIndices[s] : hm.exampleIndices[e+1]]

        # Reduce the number of class values, if class is defined
        newdomain = self.checkDomain(self.data[indx], selection=ex)
        if not newdomain:
            newdomain = self.data[indx].domain
        selectedData = orange.ExampleTable(newdomain, ex)
        return selectedData

    def sendOne(self, data, indx):
        self.send("Examples", data)
    
    def sendData(self, rows, indxs):
        indxs.sort()

        newdata = [None] * (max(indxs) + 1)
        for i in indxs:
            newdata[i] = self.prepareData(rows, i)
            newdata[i].name = self.data[i].name
            if not hasattr(self.data[i], "strain"):
                newdata[i].strain = "NoName (%d)" % i
            else:
                newdata[i].strain = self.data[i].strain
            
            self.sendOne(newdata[i], i)
            
        groups = {}
        for i in indxs:
            groups[newdata[i].strain] = []
        for i in indxs:
            groups[newdata[i].strain].append(i)
        strains = groups.keys()
        strains.sort()
        datafiles = []
        for s in strains:
            datafiles.append( (s, [newdata[i] for i in groups[s]]) )
        datafiles
        self.send("Structured Data", datafiles)
        
    ##########################################################################
    # callback functions

    def getGammaCorrectedPalette(self):
        return [QColor(*self.contPalette.getRGB(float(i)/250, gamma=self.Gamma)).rgb() for i in range(250)] + self.palette[-6:]

    def setColor(self, index, dialog=None, update=True):
        self.selectedSchemaIndex = index
        if not dialog:
            dialog = self.createColorDialog()

        self.colorCombo.setPalettes("palette", dialog)
        self.colorCombo.setCurrentIndex(self.selectedSchemaIndex)
        self.contPalette = palette = dialog.getExtendedContinuousPalette("palette")
        unknown = dialog.getColor("unknown").rgb()
        underflow = dialog.getColor("underflow").rgb()
        overflow = dialog.getColor("overflow").rgb()
        self.palette = [QColor(*palette.getRGB(float(i)/250, gamma=self.Gamma)).rgb() for i in range(250)] + [qRgb(255, 255, 255)]*3 +[underflow, overflow, unknown]

        if update:        
            self.drawHeatMap()
        
    def openColorDialog(self):
        dialog = self.createColorDialog()
        if dialog.exec_():
            self.colorSettings = dialog.getColorSchemas()
            self.selectedSchemaIndex = dialog.selectedSchemaIndex
            self.colorCombo.setCurrentIndex(self.selectedSchemaIndex)
            self.setColor(self.selectedSchemaIndex, dialog)

    def createColorDialog(self):
        c = OWColorPalette.ColorPaletteDlg(self, "Color Palette")
        c.createExtendedContinuousPalette("palette", "Continuous Palette", initialColor1=QColor(Qt.blue), initialColor2=QColor(255, 255, 0).rgb(), extendedPassThroughColors = ((Qt.red, 1), (Qt.darkYellow, 1), (Qt.black, 1), (Qt.magenta, 1), (Qt.green, 1)))
        box = c.createBox("otherColors", "Other Colors")
        c.createColorButton(box, "unknown", "Unknown", Qt.gray)
        box.layout().addSpacing(5)
        c.createColorButton(box, "overflow", "Overflow", Qt.black)
        box.layout().addSpacing(5)
        c.createColorButton(box, "underflow", "Underflow", Qt.white)
        c.setColorSchemas(self.colorSettings, self.selectedSchemaIndex)
        return c


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
        

    # ########################################################################
    # drawing

    def drawLegend(self, x, y, width, height, palette):
        legend = self.heatmapconstructor[0].getLegend(width, height, 1.0) #self.Gamma)

        lo = self.CutEnabled and self.CutLow   or self.lowerBound
        hi = self.CutEnabled and self.CutHigh  or self.upperBound

        t = QGraphicsSimpleTextItem("%3.1f" % lo, None, self.scene) #QCanvasText("%3.1f" % lo, self.canvas)
        t.setPos(x, y) #setX(x); t.setY(y)
        t.show()
        t = QGraphicsSimpleTextItem("%3.1f" % hi, None, self.scene) #QCanvasText("%3.1f" % hi, self.canvas)
        t.setPos(x+width-t.boundingRect().width(), y) #setX(x+width-t.boundingRect().width()); t.setY(y)
        t.show()
        y += t.boundingRect().height()+1
        self.legendItem = ImageItem(legend, self.scene, width, height, palette, x=x, y=y)
        return y + c_legendHeight + c_spaceY

    def drawFileName(self, label, x, y, width):
        t = QGraphicsSimpleTextItem(label, None, self.scene)
        t.setPos(x, y) #setX(x); t.setY(y)
        t.show()
        line = QGraphicsLineItem(None, self.scene)
        line.setPoints(0, 0, width, 0)
        y += t.boundingRect().height()
        line.setPos(x, y) #setX(x); line.setY(y)
        line.show()
        return y + 5

    def drawGroupLabel(self, label, x, y, width):
        t = QGraphicsSimpleTextItem(label, None, self.scene)
        t.setPos(x, y) #(x); t.setY(y)
        t.show()
        return y + t.boundingRect().height() + 1

    def drawGeneAnnotation(self, x, y, group):
        # determine the appropriate font width for annotation
        # this part is ugly, we need to do it computationally
        font = QFont()
##        dummy = QGraphicsSimpleTextItem("dummy", None)
##        last = 2; offset = 0
##        for fsize in range(2,9):
##            font.setPointSize(fsize)
##            dummy.setFont(font)
##            if dummy.boundingRect().height() > self.CellHeight:
##                break
##            offset = (self.CellHeight - dummy.boundingRect().height())/2
##            last = fsize
##        font.setPointSize(last)
##        y += offset
        font.setPixelSize(max(self.CellHeight - 1, 1))

        # annotate
        hm = self.heatmaps[0][group]
        if self.BAnnotationVar:
            for (row, indices) in enumerate(hm.exampleIndices[:-1]):
                t = QGraphicsSimpleTextItem(str(hm.examples[hm.exampleIndices[row]][self.BAnnotationVar]), None, self.scene)
                t.setFont(font)
                t.setPos(x, y)
                t.show()
                y += self.CellHeight

    def drawColumnLabels(self, x, y, heatmap):
        font = QFont()
        font.setPixelSize(min(max(self.CellWidth - 1, 1), 11))
        t = QGraphicsSimpleTextItem("Dummy123", None, self.scene)
        t.setFont(font)
        if t.boundingRect().height() < self.CellWidth:
            x += (self.CellWidth - t.boundingRect().height()) / 2
        self.scene.removeItem(t)

        maxY = y
        items = []
        if self.ShowColumnLabels:
            angle = -90 #-90 if self.ColumnLabelPosition == 0 else 90
            for attr in heatmap[0].examples.domain.attributes:
                t = QGraphicsSimpleTextItem(str(attr.name), None, self.scene)
                t.setFont(font)
                t.setPos(x, y)
                t.show()
                x += self.CellWidth
                maxY = max(y + t.boundingRect().width(), maxY)
                t.rotate(angle)
                items.append(t)

        for item in items:
            if self.ColumnLabelPosition == 0:
                item.setPos(item.x(), maxY)
            else:
                item.setPos(item.x(), item.y() + item.boundingRect().width())

        return maxY                

    def drawHeatMap(self):
        # remove everything from current canvas
        for i in self.scene.items():
            self.scene.removeItem(i)
        if not len(self.data):
            return

        lo = self.CutEnabled and self.CutLow   or self.lowerBound
        hi = self.CutEnabled and self.CutHigh  or self.upperBound

##        self.sceneView.heatmapParameters(self, self.CellWidth, self.CellHeight) # needed for event handling
        self.scene.heatmapParameters(self, self.CellWidth, self.CellHeight) # needed for event handling

##        palette = self.ColorPalettes[self.CurrentPalette]
        palette = self.getGammaCorrectedPalette() if self.Gamma !=0 else self.palette
        groups = (not self.data[0].domain.classVar and 1) or \
                 len(self.data[0].domain.classVar.values) # mercy! (just had to do this)

        self.bmps = []; self.heights = []; self.widths = []; self.imgStart = []; self.imgEnd = []
        for (i,hm) in enumerate(self.heatmaps):
            bmpl = []
            for g in range(groups):
                bmp, self.imageWidth, imageHeight = hm[g].getBitmap(int(self.CellWidth), \
                    int(self.CellHeight), lo, hi, 1.0) #self.Gamma)
                bmpl.append(bmp)
                if not i: self.heights.append(imageHeight)
            self.bmps.append(bmpl)
            self.widths.append(self.imageWidth)

        totalHeight = max(reduce(lambda x,y:x+y, self.heights) + 500, 2000)
        self.scene.setSceneRect(0, 0, 2000, totalHeight) # this needs adjustment
        x = c_offsetX; y0 = c_offsetY

        self.legend = self.heatmapconstructor[0].getLegend(self.imageWidth, c_legendHeight, 1.0) #self.Gamma)
        if self.LegendOnTop:
            y0 = self.drawLegend(x, y0, self.imageWidth, c_legendHeight, palette)

        self.heatmapPositionsX = [] # start and end positions of heatmaps
        for i in range(len(self.data)):
            y = y0; y1 = y0
            if self.ShowDataFileNames and len(self.data)>1:
                y1 = self.drawFileName(self.data[i].name, x, y, \
                    self.imageWidth+self.ShowAverageStripe*(c_averageStripeWidth + c_spaceAverageX))
            x0 = x                    
            # plot the heatmap (and group label)
            showClusters = (i == 0 and self.groupClusters and self.ShowClustering)
            ycoord = []
            y = y1; x += self.ShowAverageStripe * (c_averageStripeWidth + c_spaceAverageX)

            if self.ColumnLabelPosition == 0:
                y = self.drawColumnLabels(x, y, self.heatmaps[i]) + 2
            
            self.heatmapPositionsX.append((x, x + self.widths[i]-1))
            for g in range(groups):
              if self.heights[g]:
                if self.ShowGroupLabel and groups>1:
                    y = self.drawGroupLabel(self.data[i][0].domain.classVar.values[g], x, y, self.imageWidth)          
                if not i: self.imgStart.append(y)
                ycoord.append(y)
                if showClusters:
                    item = HierarchicalClusterItem(self.groupClusters[g], None, self.scene)
                    item.setTransform(QTransform().scale(100.0/item.rect().height(), self.heights[g]/float(len(item.cluster))).\
                                    rotate(90).translate(0, -item.rect().height()))
                    item.setPos(-100, y+self.CellHeight/2.0)
                    
                image = ImageItem(self.bmps[i][g], self.scene, self.imageWidth, \
                                  self.heights[g], palette, x=x, y=y, z=z_heatmap)
                image.hm = self.heatmaps[i][g] # needed for event handling
                image.height = self.heights[g]; image.width = self.imageWidth
                if not i: self.imgEnd.append(y+self.heights[g]-1)
                y += self.heights[g] + c_spaceY
            
            if self.ColumnLabelPosition == 1:
                self.drawColumnLabels(x, y - c_spaceY + 2, self.heatmaps[i])
            x = x0
            # plot stripe with averages
            if self.ShowAverageStripe:
                for g in range(groups):
                    avg, avgWidth, avgHeight = self.heatmaps[i][g].getAverages(c_averageStripeWidth, \
                        int(self.CellHeight), lo, hi, 1.0) # self.Gamma)
                    ImageItem(avg, self.scene, avgWidth, avgHeight, palette, x=x, y=ycoord[g])
            x += self.imageWidth + self.SpaceX + self.ShowAverageStripe * \
                (c_averageStripeWidth + c_spaceAverageX)

        # plot the gene annotation
        for g in range(groups):
            if self.ShowGeneAnnotations and self.CellHeight>4:
                self.drawGeneAnnotation(x, ycoord[g], g)

        self.selection.redraw()
        self.scene.setSceneRect(self.scene.itemsBoundingRect())
        self.scene.update()
        
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
            self.fileLB.addItem(self.createListItem(data[i].name, i))
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

    def fileSelectionChanged(self, sel):
        # self.fileRef.setDisabled(sel==0)
        self.selectedFile = sel
        self.fileDown.setEnabled(sel < len(self.data)-1)
        self.fileUp.setEnabled(sel>0)
        self.fileRef.setEnabled(sel <> self.refFile)

    def createListItem(self, text, position):
        pixmap = QPixmap(14, 13)
        pixmap.fill(Qt.white)
        
        if position == self.refFile:
            painter = QPainter()
            painter.begin(pixmap)
            painter.setPen(Qt.black)
            painter.setBrush(Qt.black)
            painter.drawRect(3, 3, 8, 8)
            painter.end()
            
        listItem = QListWidgetItem(QIcon(pixmap), text) #QListBoxPixmap(pixmap)
##        listItem.setText(text)
        return listItem

    # remove gene selection (when changing some of the options)
    def removeSelection(self):
        if self.selection:
            self.selection.remove()
            self.scene.update()

##################################################################################################
# new canvas items

class ImageItem(QGraphicsRectItem):
    def __init__(self, bitmap, scene, width, height, palette, depth=8, numColors=256, x=0, y=0, z=0):
        QGraphicsRectItem.__init__(self, None, scene)
        self.image = QImage(bitmap, width, height, QImage.Format_Indexed8)
        self.image.bitmap = bitmap # this is tricky: bitmap should not be freed,
                                   # else we get mess. hence, we store it in the object
        self.image.setColorTable(signedPalette(palette))
        #self.pixmap = QPixmap()
        #self.pixmap.convertFromImage(image, QPixmap.Color)
        self.setRect(0, 0, width, height)
        self.setPos(x, y) #setX(x); self.setY(y); self.setZ(z)
        self.setZValue(z)
        self.show()

    def paint(self, painter, options, widget=None):
        x, y, w, h = options.exposedRect.x(), options.exposedRect.y(), options.exposedRect.width(), options.exposedRect.height()
        painter.drawImage(x, y, self.image, x, y, w, h)


        
# ################################################################################################
# mouse event handler

v_sel_width = 2

class HeatMapGraphicsScene(QGraphicsScene):
    def __init__(self, *args):
        QGraphicsScene.__init__(self, *args)
        self.clicked = False
        self.shiftPressed = False
        self.currentHighlightedCluster = None
        self.selectedClusters = []

    def heatmapParameters(self, master, cellWidth, cellHeight):
        self.master = master
        self.dx, self.dy = cellWidth, cellHeight
        self.selector = QGraphicsRectItem(0, 0, \
            self.dx + 2 * v_sel_width - 1, self.dy + 2 * v_sel_width - 1, None, self)
        self.selector.setPen(QPen(self.master.SelectionColors[self.master.CurrentPalette], v_sel_width))
        self.selector.setZValue(10)

    def mouseMoveEvent(self, event):
        QGraphicsScene.mouseMoveEvent(self, event)
        # handling of selection
        if self.clicked:
            self.master.selection(self.clicked, (event.scenePos().x(), event.scenePos().y()))

        # balloon handling
        try:
            if self.master <> None and not self.master.BShowballoon: return
        except:
            return
        item = self.itemAt(event.scenePos())
        if self.currentHighlightedCluster and self.currentHighlightedCluster != item:
            self.currentHighlightedCluster.setHighlight(False)
        if isinstance(item, HierarchicalClusterItem):
            item.setHighlight(True)
            self.currentHighlightedCluster = item
        items = filter(lambda ci: ci.zValue()==z_heatmap, self.items(event.scenePos()))
        if len(items) == 0: # mouse over nothing special
            self.selector.hide()
            self.update()
        else:
            item = items[0]
            hm = item.hm
            x, y = event.scenePos().x() - item.x(), event.scenePos().y() - item.y()
            if x<0 or y<0 or x>item.width-1 or y>item.height-1: 
                self.selector.hide()
                self.scene.update()
                return
            col, row = int(x / self.dx), int(y / self.dy)
            # hm.getCellIntensity(row, col), hm.getRowIntensity(row)
            ex = hm.examples[hm.exampleIndices[row] : hm.exampleIndices[row+1]]
            self.selector.setPos(item.x()+col*self.dx-v_sel_width+1, item.y()+row*self.dy-v_sel_width+1)
            self.selector.show()

            # bubble, construct head
            if hm.getCellIntensity(row, col)!=None:
                head = "%6.4f" % hm.getCellIntensity(row, col)
            else:
                head = "Missing Data"
            if self.master.BShowColumnID:
                head += "\n"+ex[0].domain.attributes[col].name
            # bubble, construct body
            body = None
            if (self.master.BShowSpotIndex and self.master.BSpotVar) or \
                    self.master.BShowAnnotation or self.master.BShowGeneExpression:
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

            QToolTip.showText(QPoint(event.screenPos().x(), event.screenPos().y()), "")
            QToolTip.showText(QPoint(event.screenPos().x(), event.screenPos().y()), head + body)

    def keyPressEvent(self, e):
        if e.key() == 4128:
            self.shiftPressed = True
        else:
            QGraphicsScene.keyPressEvent(self, e)

    def keyReleaseEvent(self, e):        
        if e.key() == 4128:
            self.shiftPressed = False
        else:
            QGraphicsScene.keyReleaseEvent(self, e)

    def mousePressEvent(self, event):
        # self.viewport().setMouseTracking(False)
        item = self.itemAt(event.scenePos())
        if isinstance(item ,HierarchicalClusterItem):
            leaves = list(item)
            print len(leaves)
            first, last = leaves[0], leaves[-1]
            first_x, first_y = first.mapToScene(first.rect().topLeft()).x(), first.mapToScene(first.rect().topLeft()).y()
            last_x, last_y = last.mapToScene(last.rect().topLeft()).x(), last.mapToScene(last.rect().topLeft()).y()
            self.master.selection.start((first_x, first_y), (first_x, first_y), self.shiftPressed)
            self.master.selection((first_x, first_y), (last_x, last_y))
            self.master.selection.release()
            return
        self.clicked = (event.scenePos().x(), event.scenePos().y())
        if not self.master.selection.start(self.clicked, self.clicked, self.shiftPressed):
            self.clicked = None

    def mouseReleaseEvent(self, event):
        if self.clicked:
            self.clicked = False
            self.update()
            self.master.selection.release()
            

class MyGraphicsView(QGraphicsView):
    def __init__(self, *args):
        QGraphicsView.__init__(self, *args)
        self.clicked = False
        self.viewport().setMouseTracking(True)
        self.setFocusPolicy(Qt.ClickFocus)
        self.setFocus()
        self.shiftPressed = False

    def keyPressEvent(self, e):
        if e.key() == Qt.Key_Shift:
            self.scene().shiftPressed = True
        else:
            QGraphicsView.keyPressEvent(self, e)

    def keyReleaseEvent(self, e):        
        if e.key() == Qt.Key_Shift:
            self.scene().shiftPressed = False
        else:
            QGraphicsView.keyReleaseEvent(self, e)

# ################################################################################################
# data selection

def interval_position(x, l, forced=None):
    if forced:
        for (i,(min,max)) in enumerate(l):
            if x>=min and x<=max:
                return i
        return None
    else:
        for (i,(min,max)) in enumerate(l):
            if x>=min and x<=max:
                return i
            if x<min:
                if i>forced:
                    if i-1>0: return i-1
                    return 0
                else:
                    return i
        return i

class SelectData(object):
    def __init__(self, master, scene):
        self.scene = scene
        self.master = master
        self.add = False
        self.squares = []; self.rows = []    # cumulative, used for multiple selection
        self.cRows = []; self.cSquares = []  # workaround for when release(self) is called
                        # before any clicking at all has occured (maybe caused by a Qt bug?)
        self.startIndx = None; self.indxs = []

    # removes the selection and relate information
    def remove(self):
        for r in self.squares:
            self.scene.removeItem(r)
        self.squares = []; self.rows = []

    # starts the selection, called after the first click
    def start(self, p1, p2, add=False):
        indx = interval_position(p1[0], self.master.heatmapPositionsX)
        if indx == None:
            self.remove()
            self.indx = []; self.startIndx = None
            return False
        
        if not add:
            self.remove()
            self.startIndx = indx
            self.indxs = [indx]
        else:
            if self.master.SelectionType==0:
                if self.startIndx<>None and self.startIndx<>indx:
                    self.remove()
                    self.indx = [indx]
                    self.startIndx = indx                    
            else:
                self.startIndx = indx
                if indx not in self.indxs:
                    self.indxs.append(indx)

        self.cSquares = []; self.cRows = [] # current selection
        self.add = add
        self.__call__(p1, p2)
        # determines which microarray was clicked
        return True
 
    # called during dragging (extending the selection)
    def __call__(self, p1, p2):
        for r in self.cSquares:
            self.scene.removeItem(r)
        y1 = min(p1[1], p2[1]); y2 = max(p1[1], p2[1])
        if self.master.SelectionType == 1:
            indx = interval_position(p2[0], self.master.heatmapPositionsX)
            if indx==None:
                indx = interval_position(p2[0], self.master.heatmapPositionsX, forced=self.startIndx)
            irange = range(min(indx, self.startIndx), max(indx, self.startIndx)+1)
            for i in irange:
                if indx not in self.indxs:
                    self.indxs.append(indx)

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
                self.scene.removeItem(r)
            self.rows = newrows
            squares = self.draw(self.rows)
            self.squares = squares
        else:
            self.rows = self.cRows
            self.squares = self.cSquares
        if self.rows:
            self.master.sendData(self.rows, self.indxs)

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

    # draws rectangles around selected points for a indx-th microarray
    def drawOne(self, rows, indx):
        lines = []
        start = self.master.imgStart; end = self.master.imgEnd
        self.master.heatmapPositionsX[indx][0]
        for (g, r1, r2) in rows:
            y1 = start[g] + r1 * self.master.CellHeight
            y2 = start[g] + (r2+1) * self.master.CellHeight - 1
            r = QGraphicsRectItem(self.master.heatmapPositionsX[indx][0]-v_sel_width+1, y1-v_sel_width+1, \
                    self.master.widths[indx]+2*v_sel_width-1, y2-y1+v_sel_width+2, None, self.scene)
            r.setPen(QPen(self.master.SelectionColors[self.master.CurrentPalette], v_sel_width))
            #r.setPen(QPen(Qt.red, v_sel_width))
            r.setZValue(10)
            r.show()
            lines.append(r)
        return lines
            
    # draws rectangles around selected points
    def draw(self, rows):
        if self.master.SelectionType==0: # selection on single data set
            lines = self.drawOne(rows, self.startIndx)
        else: # selection on multiple data set
            lines = []
            for i in self.indxs:
                l = self.drawOne(rows, i)
                lines += l
        self.scene.update()
        return lines

    def redraw(self):
        for r in self.squares:
            self.scene.removeItem(r)
        if self.rows:
            self.squares = self.draw(self.rows)


##################################################################################################
# test script

if __name__=="__main__":
    import OWDataFiles, orngSignalManager
    signalManager = orngSignalManager.SignalManager(0)
    a=QApplication(sys.argv)
    ow=OWHeatMap(signalManager = signalManager)
    signalManager.addWidget(ow)
    a.setMainWidget(ow)
    ow.show()
    ds = OWDataFiles.OWDataFiles(signalManager = signalManager, loaddata=0)
    signalManager.addWidget(ds)
    ds.loadData("potato.sub100")
#    ds.loadData("yakApufA")
#    ds.loadData("smallchipdata")
#    ds.loadData("testchip")
    signalManager.setFreeze(1)
    signalManager.addLink(ds, ow, 'Structured Data', 'Structured Data', 1)
    signalManager.setFreeze(0)
    a.exec_loop()
    ow.saveSettings()
    
##    d = orange.ExampleTable('wt'); d.name = 'wt'
##    d = orange.ExampleTable('wt-nometa'); d.name = 'wt'
##    ow.dataset(d, 2)
##    d = orange.ExampleTable('wt-nometa'); d.name = 'wt'
##    ow.dataset(d, 1)
##    ow.dataset(None, 1)
##    ow.dataset(None, 2)

##    names = ['wt1', 'wt2', 'wt3', 'wt4']
##    for i, s in enumerate(names): 
##        d = orange.ExampleTable(s); d.name = s
##        ow.dataset(d, i)
##    ow.dataset(None, 2)

#    a.exec_loop()
#    ow.saveSettings()
