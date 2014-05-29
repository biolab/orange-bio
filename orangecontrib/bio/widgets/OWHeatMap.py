"""
<name>Heat Map</name>
<description>Heatmap visualization.</description>
<contact>Ales Erjavec, Blaz Zupan, Janez Demsar</contact>
<icon>icons/Heatmap.svg</icon>
<priority>1040</priority>
"""

from __future__ import absolute_import

from collections import defaultdict
import itertools
import math

import orange
import orangene
from Orange.orng import orngClustering
from Orange.OrangeWidgets import OWColorPalette
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.ColorPalette import signedPalette
from Orange.OrangeWidgets.OWClustering import HierarchicalClusterItem, DendrogramWidget, DendrogramItem
from Orange.OrangeWidgets.OWDlgs import OWChooseImageSizeDlg
from Orange.OrangeWidgets.OWWidget import *

import functools

NAME = "Heat Map"
DESCRIPTION = "Heatmap visualization."
ICON = "icons/Heatmap.svg"
PRIORITY = 1040

INPUTS = [("Examples", Orange.data.Table, "set_dataset")]
OUTPUTS = [("Examples", Orange.data.Table, Single + Default)]

REPLACES = ["_bioinformatics.widgets.OWHeatMap.OWHeatMap"]


DEBUG = False

#BUG: OWHeatMap does not support heatmaps which need a image which is larger than maxint X maxint pixels!
#and no warning appears

import warnings
warnings.filterwarnings("ignore", "'strain'", orange.AttributeWarning)

def split_domain(domain, split_label):
    """ Split the domain based on values of `split_label` value.
    """
    groups = defaultdict(list)
    for attr in domain.attributes:
        groups[attr.attributes.get(split_label)].append(attr)
        
    attr_values = [attr.attributes.get(split_label) for attr in domain.attributes]
    
    domains = []
    for value, attrs in groups.items():
        group_domain = orange.Domain(attrs, domain.class_var)
        group_domain.add_metas(domain.get_metas())
        domains.append((value, group_domain))
        
    if domains:
        assert(all(len(dom) == len(domains[0][1]) for _, dom in domains))
        
    return sorted(domains, key=lambda t: attr_values.index(t[0]))
    
def vstack_by_subdomain(data, sub_domains):
    domain = sub_domains[0]
    newtable = orange.ExampleTable(domain)
    
    for sub_dom in sub_domains:
        for ex in data:
            vals = [ex[a].native() for a in sub_dom]
            newtable.append(orange.Example(domain, vals))
    
    return newtable
    
    
def select_by_class(data, class_):
    indices = select_by_class_indices(data, class_)
    return data.select(indices)

def select_by_class_indices(data, class_):
    return [1 if class_ == ex.getclass() else 0 for ex in data]

def group_by_unordered(iterable, key):
    groups = defaultdict(list)
    for item in iterable:
        groups[key(item)].append(item)
    return groups.items()

class PP_callback(object):

    def __init__(self, parts=None, progress_callback=None):
        self.progress_callback = progress_callback
        self.parts = parts
        self.part = -1

    def __call__(self, value):
        return self.progress_callback(value / self.parts + 100.0 * self.part / self.parts)

    def newpart(self, k=1):
        self.part += k
        self.__call__(0)
 

def hierarchical_cluster_ordering_examples(data, group_domains=None, opt_order=False, pp_callback=None):
       
    classVar = data.domain.classVar
    if classVar and isinstance(classVar, orange.EnumVariable):
        class_data = [select_by_class_indices(data, val) for val in data.domain.classVar.values]
    else:
        class_data = [[1] * len(data)]
 
    def indices_map(indices):
        map = zip(range(len(indices)), indices)
        map = [i for i, test in map if test]
        return dict(enumerate(map))
    
    data_ordering = []
    data_clusters = []
    for i, indices in enumerate(class_data):
        pp_callback.newpart()
        sub_data = data.select(indices)
        cluster = orngClustering.hierarchicalClustering(sub_data, order=opt_order, progressCallback=pp_callback)
        ind_map = indices_map(indices)
        data_ordering.append([ind_map[m] for m in cluster.mapping])
        data_clusters.append(cluster)
        
    return data_ordering, data_clusters 

def hierarchical_cluster_ordering_attributes(data, group_domains=None, opt_order=False, pp_callback=None):

    pp_callback.newpart()

    if group_domains is not None and len(group_domains) > 1:
        stacked = vstack_by_subdomain(data, group_domains)
    else:
        stacked = data    

    attr_cluster = orngClustering.hierarchicalClustering_attributes(stacked, order=opt_order, progressCallback=pp_callback)
    attr_ordering = list(attr_cluster.mapping)
       
    return attr_ordering, attr_cluster

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

class ExampleTableContextHandler(ContextHandler):
    def match(self, context, imperfect, examples):
        return context.checksum == examples.checksum() and 2
    
    def findOrCreateContext(self, widget, examples):
        context, isNew = ContextHandler.findOrCreateContext(self, widget, examples)
        if not context:
            return None, False
        context.checksum = examples.checksum()
        return context, isNew
    
from .OWGenotypeDistances import SetContextHandler

class OWHeatMap(OWWidget):
    contextHandlers = {"": DomainContextHandler("", ["CellWidth", "CellHeight"]),
                       "Selection": ExampleTableContextHandler("Selection", contextDataVersion=2)}
    
    settingsList = ["CellWidth", "CellHeight", "SpaceX", "Merge",
                    "Gamma", "CutLow", "CutHigh", "CutEnabled", 
                    "ShowAnnotation", "LegendOnTop", "LegendOnBottom",
                    "ShowAverageStripe", "ShowGroupLabel",
                    "MaintainArrayHeight",
                    "GShowToolTip", "BShowColumnID", "BShowSpotIndex",
                    "BShowAnnotation", 'BShowGeneExpression',
                    "BSpotVar", "ShowGeneAnnotations",
                    "ShowDataFileNames", "BAnnotationVar",
                    "SelectionType",
                    "SortExamples", "SortAttributes",
                    "CurrentPalette", "colorSettings", "selectedSchemaIndex",
                    "palette", "ShowColumnLabels", "ColumnLabelPosition"]

    def __init__(self, parent=None, signalManager = None):
        OWWidget.__init__(self, parent, signalManager, 'HeatMap', TRUE)
        
        self.inputs = [("Examples", ExampleTable, self.set_dataset)]
        self.outputs = [("Examples", ExampleTable, Default)]

        #set default settings
        self.CellWidth = 3; self.CellHeight = 3
        self.SpaceX = 10
        self.Merge = 1; self.savedMerge = self.Merge
        self.Gamma = 1
        self.CutLow = self.CutHigh = self.CutEnabled = 0
        self.ShowAnnotation = 0
        self.LegendOnTop = 0           # legend stripe on top (bottom)?
        self.LegendOnBottom = 1
        self.ShowLegend = 1
        self.ShowGroupLabel = 1        # show class names in case of classified data?
        self.ShowAverageStripe = 0     # show the stripe with the evarage
        self.MaintainArrayHeight = 0   # adjust cell height while changing the merge factor
        self.GShowToolTip = 1          # balloon help
        self.ShowGeneAnnotations = 1   # show annotations for genes
        self.ShowColumnLabels = 1
        self.ColumnLabelPosition = 0
        self.ShowDataFileNames = 1     # show the names of the data sets (only in case of multiple files)
        self.BShowColumnID = 1; self.BShowSpotIndex = 1; self.BShowAnnotation = 1; self.BShowGeneExpression = 1
        self.BSpotVar = None; self.BAnnotationVar = None  # these are names of variables
        self.BSpotIndx = None; self.BAnnotationIndx = None # these are id's of the combo boxes
        self.ShowClustering = 1
        self.SelectionType = 0         # selection on a single data set
        self.setColorPalette()
        self.refFile = 0               # position index of a reference file
        self.selectedFile = None       # position of the selected file in the list box

        self.colorSettings =None
        self.selectedSchemaIndex = 0
        self.auto_commit = True
        self.selection_changed_flag = False

        self.palette = self.ColorPalettes[0]

        self.SortExamples = 0
        self.SortAttributes = 0

        self.loadSettings()
        self.data = []
        self.maxHSize = 30; self.maxVSize = 15

        

        # GUI definition
        self.connect(self.graphButton, SIGNAL("clicked()"), self.saveFig)
        self.tabs = OWGUI.tabWidget(self.controlArea)

        # SETTINGS TAB
        settingsTab = OWGUI.createTabPage(self.tabs, "Settings")
        box = OWGUI.widgetBox(settingsTab, "Cell Size (Pixels)", addSpace=True)
        OWGUI.qwtHSlider(box, self, "CellWidth", label='Width: ',
                         labelWidth=38, minValue=1, maxValue=self.maxHSize,
                         step=1, precision=0, callback=self.update_cell_size)
        
#        OWGUI.hSlider(box, self, "CellWidth", label="Width:", minValue=1,
#                      maxValue=self.maxHSize, step=1, ticks=5,
#                      callback=self.update_cell_size,
#                      tooltip="Width of each heatmap cell.")
        
        self.sliderVSize = OWGUI.qwtHSlider(box, self, "CellHeight",
                                label='Height: ', labelWidth=38,
                                minValue=1, maxValue=self.maxVSize,
                                step=1, precision=0,
                                callback=self.update_cell_size)
#        self.sliderVSize = OWGUI.hSlider(box, self, "CellHeight",
#                      label='Height:', 
#                      minValue=1, maxValue=50, ticks=5,
#                      callback = self.update_cell_size)
        
        OWGUI.qwtHSlider(box, self, "SpaceX", label='Space: ', 
                         labelWidth=38, minValue=0, maxValue=50,
                         step=2, precision=0,
                         callback=self.update_grid_spacing)
#                         callback=self.drawHeatMap)
        
        OWGUI.qwtHSlider(settingsTab, self, "Gamma", box="Gamma",
                         minValue=0.1, maxValue=1, step=0.1,
                         callback=self.update_color_schema)
#                         callback=self.drawHeatMap)
        
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
        button = OWGUI.button(box, self, "Edit colors", callback=self.openColorDialog, tooltip="Edit the heatmap color palette", debuggingEnabled=0)
        
        OWGUI.separator(settingsTab)
        
        
        boxy = OWGUI.widgetBox(settingsTab, "Sorting", orientation = "vertical", addSpace=True)

        # For attributes
        OWGUI.comboBox(boxy, self, "SortAttributes", 
                        items=["No sorting", "Clustering",
                               "Clustering with leaf ordering"], label='X axis',
                               callback=self.update_sorting_attributes)
 
        # For examples
        OWGUI.comboBox(boxy, self, "SortExamples",
                        items=["No sorting", "Sort examples", "Clustering",
                               "Clustering with leaf ordering"], label='Y axis',
                               callback=self.update_sorting_examples)
       
        OWGUI.rubber(settingsTab)
        
        # FILTER TAB
        tab = OWGUI.createTabPage(self.tabs, "Filter")
        box = OWGUI.widgetBox(tab, "Threshold Values", addSpace=True)
        OWGUI.checkBox(box, self, 'CutEnabled', "Enabled",
                       callback=self.update_thresholds)
        
        self.sliderCutLow = OWGUI.qwtHSlider(box, self, 'CutLow', label='Low:', 
                            labelWidth=40, minValue=-100, maxValue=0, step=0.1,
                            precision=1, ticks=0, maxWidth=80,
                            callback=self.update_thresholds)
        
        self.sliderCutHigh = OWGUI.qwtHSlider(box, self, 'CutHigh', label='High:', 
                            labelWidth=40, minValue=0, maxValue=100, step=0.1,
                            precision=1, ticks=0, maxWidth=80,
                            callback=self.update_thresholds)
        
        if not self.CutEnabled:
            self.sliderCutLow.box.setDisabled(1)
            self.sliderCutHigh.box.setDisabled(1)

        box = OWGUI.widgetBox(tab, "Merge", addSpace=True)
##        OWGUI.qwtHSlider(box, self, "Merge", label='Rows:', labelWidth=33, minValue=1, maxValue=500, step=1, callback=self.mergeChanged, precision=0, ticks=0)
        OWGUI.spin(box, self, "Merge", min=1, max=500, step=1, label='Rows:',
#                   callback=self.mergeChanged,
                   callback=self.on_merge_changed,
                   callbackOnReturn=True)
        OWGUI.checkBox(box, self, 'MaintainArrayHeight', "Maintain array height")
        OWGUI.rubber(tab)

        # INFO TAB
        tab = OWGUI.createTabPage(self.tabs, "Info")

        box = OWGUI.widgetBox(tab,'Annotation && Legends')
        OWGUI.checkBox(box, self, 'ShowLegend', 'Show legend', 
                       callback=self.update_legend)
        OWGUI.checkBox(box, self, 'ShowAverageStripe', 'Stripes with averages', 
                       callback=self.update_averages_stripe)
        self.geneAnnotationsCB = OWGUI.checkBox(box, self, 'ShowGeneAnnotations', 'Gene annotations', 
                                                callback=self.update_annotations)
        
        self.annotationCombo = OWGUI.comboBox(box, self, "BAnnotationIndx", items=[],
                                              callback=self.update_annotations) 
                                              #callback=lambda x='BAnnotationVar', y='BAnnotationIndx': self.setMetaID(x, y))

        box = OWGUI.widgetBox(tab, 'Column Labels')
        columnLabelCB = OWGUI.checkBox(box, self, "ShowColumnLabels", "Display column labels",
                                       callback=self.update_column_annotations)
        posbox = OWGUI.widgetBox(OWGUI.indentedBox(box), "Position", flat=True)
        comboBox = OWGUI.comboBox(posbox, self, "ColumnLabelPosition", 
                                  items=["Top", "Bottom"], 
                                  callback=self.update_column_annotations)
        columnLabelCB.disables.append(comboBox.box)
        columnLabelCB.makeConsistent()
        
        box = OWGUI.widgetBox(tab, "Tool Tips")
        cb = OWGUI.checkBox(box, self, 'GShowToolTip', "Show tool tips")
        box = OWGUI.widgetBox(OWGUI.indentedBox(box), "Tool Tip Info")
        box.setFlat(True)
        OWGUI.checkBox(box, self, 'BShowColumnID', "Column ID")
        self.spotIndxCB = OWGUI.checkBox(box, self, 'BShowSpotIndex', "Spot Index", \
            callback=lambda: self.spotCombo.setDisabled(not self.BShowSpotIndex))
        self.spotCombo = OWGUI.comboBox(box, self, "BSpotIndx", items=[], \
            callback=lambda x='BSpotVar', y='BSpotIndx': self.setMetaID(x, y))
        OWGUI.checkBox(box, self, 'BShowGeneExpression', "Gene expression")
        OWGUI.checkBox(box, self, 'BShowAnnotation', "Annotation")
        self.toolTipInfoBox = box
        cb.disables.append(box)
        cb.makeConsistent()
        OWGUI.rubber(tab)

        # SPLIT TAB
        self.splitTab = OWGUI.createTabPage(self.tabs, "Split Data")
        box = OWGUI.widgetBox(self.splitTab, "Split By")
        self.splitLB = QListWidget(box)
        box.layout().addWidget(self.splitLB)
        self.connect(self.splitLB, SIGNAL("itemSelectionChanged()"), self.split_changed)

        # Scene with microarray
        self.heatmap_scene = self.scene = HeatmapScene()
        self.selection_manager = HeatmapSelectionManager(self)
        self.connect(self.selection_manager, SIGNAL("selection_changed()"), self.on_selection_changed)
        self.connect(self.selection_manager, SIGNAL("selection_finished()"), self.on_selection_finished)
        self.heatmap_scene.set_selection_manager(self.selection_manager)
        item = QGraphicsRectItem(0, 0, 10, 10, None, self.heatmap_scene)
        self.heatmap_scene.itemsBoundingRect()
        self.heatmap_scene.removeItem(item)
        
        self.sceneView = QGraphicsView(self.scene)
        self.currentHighlightedCluster = None
        self.selectedClusters = []
        self.mainArea.layout().addWidget(self.sceneView)
        self.heatmap_scene.widget = None
        self.heatmap_widget_grid = [[]]
        self.attr_annotation_widgets = []
        self.attr_dendrogram_widgets = []
        self.gene_annotation_widgets = []
        self.gene_dendrogram_widgets = []
        
        self.heatmaps = []
        
        self.selection_rects = []
        self.selected_rows = []

        self.attr_cluster = None
        self.data_clusters = []
        self.sorted_data = None
        
        self._ordering_cache = {}
        
        self.resize(800,400)

    def createColorStripe(self, palette):
        dx = 104; dy = 18
        bmp = chr(252)*dx*2 + reduce(lambda x,y:x+y, \
           [chr(i*250/dx) for i in range(dx)] * (dy-4)) + chr(252)*dx*2 
        
        image = QImage(bmp, dx, dy, QImage.Format_Indexed8)
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
            self.update_color_schema()
#            self.drawHeatMap()
        
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
#        if val=='BAnnotationVar':
#            self.drawHeatMap()

    def setMetaCombos(self):
        self.meta = [m.name for m in self.data.domain.getmetas().values()]
        self.BSpotVar, self.BSpotIndx = self.setMetaCombo(self.spotCombo, self.BSpotVar, \
            enabled=self.BShowSpotIndex, default='RMI')
        self.BAnnotationVar, self.BAnnotationIndx = self.setMetaCombo(self.annotationCombo, \
            self.BAnnotationVar, enabled=self.BShowAnnotation, default='xannotation')
        
    def set_meta_combos(self):
        self.spotCombo.clear()
        self.annotationCombo.clear()
        
        self.meta = self.data.domain.getmetas().values()
        names = [m.name for m in self.meta]
        
        self.spotCombo.addItems(names)
        self.annotationCombo.addItems(names)
        enabled = bool(self.meta)
        
        self.spotIndxCB.setEnabled(enabled)
        self.geneAnnotationsCB.setEnabled(enabled)
        
        self.spotCombo.setEnabled(enabled and self.BShowSpotIndex)
        self.annotationCombo.setEnabled(enabled and self.BShowAnnotation)
        
        self.BSpotIndx = 0
        self.BSpotVar = self.meta[0] if self.meta else None
        self.BAnnotationIndx = 0
        self.BAnnotationVar = self.meta[0] if self.meta else None
        
    def get_candidate_splits(self):
        """ Return candidate labels on which we can split the data. 
        """
        if self.data is not None:
            groups = defaultdict(list)
            for attr in self.data.domain.attributes:
                for item in attr.attributes.items():
                    groups[item].append(attr)
                
            by_keys = defaultdict(list)
            for (key, value), attrs in groups.items():
                by_keys[key].append(attrs)
           
            # Find the keys for which all values have the same number of attributes.
            candidates = []
            for key, groups in by_keys.items():
                count = len(groups[0])
                if all(len(attrs) == count for attrs in groups) and len(groups) > 1 and count > 1:
                    candidates.append(key)
                    
            return candidates
        else:
            return []
            
    def set_split_labels(self):
        """ Set the list view in Split tab.
        """
        self.splitLB.addItems(self.get_candidate_splits())
        
    def selected_split_label(self):
        item = self.splitLB.currentItem()
        return str(item.text()) if item else None
        
    def clear(self):
        self.data = None
        self.sorted_data = None
        self.spotCombo.clear()
        self.annotationCombo.clear()
        self.splitLB.clear()
        self.meta = []
        self.clear_scene()
        self.selected_rows = []

    def clear_scene(self):
        self.selection_manager.set_heatmap_widgets([[]])
        self.heatmap_scene.clear()
        self.heatmap_scene.widget = None
        self.heatmap_widget_grid = [[]]
        self.attr_annotation_widgets = []
        self.attr_dendrogram_widgets = []
        self.gene_annotation_widgets = []
        self.gene_dendrogram_widgets = []
        
        self.selection_rects = []
        
    def saveFig(self):
        sizeDlg = OWChooseImageSizeDlg(self.scene, parent=self)
        sizeDlg.exec_()

    ##########################################################################
    # handling of input/output signals
        
    def set_dataset(self, data=None, id=None):
        self.closeContext("Selection")
        self._ordering_cache.clear()
        
        self.clear()
        self.data = data
        if data is not None:
#            self.setMetaCombos()
            self.set_meta_combos()
            self.set_split_labels()
            
        self.unorderedData = None
        self.groupClusters = None
        
    def handleNewSignals(self):
        if self.data:
            self.update_heatmaps()
        else:
            self.clear()

        if self.data:
            self.openContext("Selection", self.data)

        if self.selected_rows:
            self.commit()
        else:
            self.send('Examples', None)

    def construct_heatmaps(self, data, split_label=None):
        if split_label is not None:
            groups = split_domain(data.domain, split_label)
        else:
            groups = [("", data.domain)]
            
        group_domains = [dom for _, dom in groups]

        attr_ordering = range(len(group_domains[0][1].attributes))
        attr_cluster = None
        data_ordering = []
        data_clusters = [None]
        sorted_data = data

        self.progressBarInit()

        progress_parts = 0

        args_key_examples = tuple(tuple(d) for d in group_domains), self.SortExamples == 3, "data"
        args_key_attributes = tuple(tuple(d) for d in group_domains), self.SortAttributes == 2, "attributes"

        if self.SortExamples > 1 and args_key_examples not in self._ordering_cache:
            classVar = data.domain.class_var
            if classVar and isinstance(classVar, orange.EnumVariable):
                progress_parts += len(classVar.values)
            else:
                progress_parts += 1

        if self.SortAttributes > 0 and args_key_attributes not in self._ordering_cache:
            progress_parts += 1

        progress_bar = PP_callback(progress_callback=self.progressBarSet, parts=progress_parts)

        # rows
        if self.SortExamples > 1:

            cluster_ordering = self._ordering_cache.get(args_key_examples, None)
            if cluster_ordering is None:

                # Rows separately
                data_ordering, data_clusters = \
                        hierarchical_cluster_ordering_examples(data, group_domains,
                                      opt_order=self.SortExamples == 3,
                                      pp_callback=progress_bar)

                # Cache the clusters
                self._ordering_cache[args_key_examples] = (data_ordering, data_clusters)
            else:
                data_ordering, data_clusters = cluster_ordering
            
            sorted_data = [data[i] for i in itertools.chain(*data_ordering)]
        
        # columns
        if self.SortAttributes > 0:

            cluster_ordering = self._ordering_cache.get(args_key_attributes, None)
            if cluster_ordering is None:

                # Columns separately
                attr_ordering, attr_cluster = \
                        hierarchical_cluster_ordering_attributes(data, group_domains,
                                      opt_order=self.SortAttributes == 2,
                                      pp_callback=progress_bar)

                # Cache the clusters
                self._ordering_cache[args_key_attributes] = (attr_ordering, attr_cluster)
            else:
                attr_ordering, attr_cluster = cluster_ordering

        self.progressBarFinished()

        self.heatmapconstructor = []
        self._group_data = []
        
        for name, group_domain in groups:
            if attr_ordering != sorted(attr_ordering):
                domain = orange.Domain([group_domain[i] for i in attr_ordering], group_domain.classVar)
                domain.addmetas(group_domain.getmetas())
                group_domain = domain
                
            group_data = orange.ExampleTable(group_domain, sorted_data)
            self._group_data.append((group_data, group_domain)) # Crashes at accessing the heatmap.examples[0] without this 
            if self.SortExamples == 1:
                hc = orangene.HeatmapConstructor(group_data)
            else:
                hc = orangene.HeatmapConstructor(group_data, None)
            
            self.heatmapconstructor.append(hc)
            
        self.attr_cluster = attr_cluster
        self.data_clusters = data_clusters
        self.sorted_data = sorted_data
        self.group_domains = groups
             
    def create_heatmaps(self, constructors):
        self.lowerBound = 1000
        self.upperBound = -1000
        squeeze = 1.0 / self.Merge
        self.heatmaps = []
        for hmc in constructors:
            hm, lb, ub = hmc(squeeze)
            
            self.lowerBound = min(self.lowerBound, lb)
            self.upperBound = max(self.upperBound, ub)
                
            self.heatmaps.append(hm)
            
        for cluster, heatmap in zip(self.data_clusters, self.heatmaps[0]):
            if cluster is not None:
                cluster._heatmap = heatmap
            
        self.sliderCutLow.setRange(self.lowerBound, 0, 0.1)
        self.sliderCutHigh.setRange(1e-10, self.upperBound, 0.1)
        self.CutLow = max(self.CutLow, self.lowerBound)
        self.CutHigh = min(self.CutHigh, self.upperBound)
        self.sliderCutLow.setValue(self.CutLow)
        self.sliderCutHigh.setValue(self.CutHigh)
            
    def point_size_hint(self, height):
        font = QFont(self.font())
        font.setPointSize(height)
        fix = 0
        while QFontMetrics(font).lineSpacing() > height and height - fix > 1:
            fix += 1
            font.setPointSize(height - fix)
        return height - fix
    
    def construct_heatmaps_scene(self, heatmaps, data, attr_cluster=None, data_clusters=None):
        self.heatmap_scene.clear()
        widget = GridWidget()
        self.heatmap_scene.addItem(widget)
        layout = QGraphicsGridLayout()
        layout.setSpacing(self.SpaceX)
        widget.setLayout(layout)
        
        classVar = data.domain.classVar
        if classVar and isinstance(classVar, orange.EnumVariable):
            classes = classVar.values
        else:
            classes = [None]
        
        if self.CutEnabled:
            cut_low, cut_high = self.CutLow, self.CutHigh
        else:
            cut_low, cut_high = self.lowerBound, self.upperBound
            
        palette = self.getGammaCorrectedPalette() if self.Gamma !=0 else self.palette
        
        class_dendrograms = []
        attr_dendrograms = []
        heatmap_widgets = []
        attr_annotation_widgets = []
        gene_annotation_widgets = []
        attr_annotation_widgets_top = []
        attr_annotation_widgets_bottom = []
        
        # Dendrograms on the left side
        if data_clusters and any(data_clusters):
            for i, cluster in enumerate(data_clusters):
                class_dendrogram = DendrogramWidget(cluster, parent=widget, orientation=Qt.Vertical)
                class_dendrogram.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
                left, top, right, bottom = class_dendrogram.layout().getContentsMargins()
                class_dendrogram.layout().setContentsMargins(left, self.CellHeight / 2.0 / self.Merge, 0.0, self.CellHeight / 2.0 / self.Merge)
                class_dendrogram.setMinimumWidth(100)
                class_dendrogram.setMaximumWidth(100)
                
                layout.addItem(class_dendrogram, i*2 + 5, 0)
                class_dendrograms.append(class_dendrogram)
            
        # Class labels    
        for i, class_ in enumerate(classes):
            if class_ is not None:
                item = GtI(class_, widget)
                item = GraphicsSimpleTextLayoutItem(item, parent=widget)
                item.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
                layout.addItem(item, i*2 + 4, 2)
                layout.setRowSpacing(i*2 + 4, 2)
                layout.setAlignment(item, Qt.AlignLeft | Qt.AlignVCenter)
                
        font = QFont()
        font.setPointSize(self.point_size_hint(self.CellHeight))
        
        class_row_labels = [map(str, hm.exampleIndices) for hm in heatmaps[0]]
        group_column_labels = [[a.name for a in hm[0].examples.domain.attributes] for hm in heatmaps]
        
        # Gene annotations on the right side
        for i, labels in enumerate(class_row_labels):
            list = GraphicsSimpleTextList(labels, parent=widget, orientation=Qt.Vertical)
            list.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
            list.setFont(font)
            list.setContentsMargins(0.0, 0.0, 0.0, 0.0)
            list.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
            
            layout.addItem(list, i*2 + 5, len(self.heatmaps) + 2)
            layout.setAlignment(list, Qt.AlignLeft)
            gene_annotation_widgets.append(list)
            
        font = QFont()
        font.setPointSizeF(self.point_size_hint(self.CellWidth))
        
        if self.ShowAverageStripe:
            stripe_offset = c_averageStripeWidth + 2
        else:
            stripe_offset = 0
            
        for column, (hm, labels, group) in enumerate(zip(heatmaps, group_column_labels, self.group_domains)):
            column_heatmap_widgets = []
            
            # Top group label
            if len(heatmaps) > 1:
                item = GtI(group[0], widget)
                item = GraphicsSimpleTextLayoutItem(item, parent=widget)
                layout.addItem(item, 1, column + 2)
                layout.setRowSpacing(1, 2)
                layout.setRowMaximumHeight(1, item.geometry().height())
                layout.setAlignment(item, Qt.AlignLeft | Qt.AlignVCenter)
            else:
                layout.setRowMaximumHeight(1, 0)
                layout.setRowSpacing(1, 0)
                
            # Top dendrogram
            if attr_cluster is not None:
                attr_dendrogram = DendrogramWidget(attr_cluster, parent=widget, orientation=Qt.Horizontal)
                attr_dendrogram.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
                attr_dendrogram.translate(0.0, attr_dendrogram.size().height())
                attr_dendrogram.scale(1.0, -1.0)
                
                left, top, right, bottom = attr_dendrogram.layout().getContentsMargins()
                attr_dendrogram.layout().setContentsMargins(stripe_offset + self.CellWidth / 2.0, 0.0, self.CellWidth / 2.0, bottom)
                attr_dendrogram.setMinimumHeight(100)
                attr_dendrogram.setMaximumHeight(100)
                
                layout.addItem(attr_dendrogram, 2, column + 2)
                layout.setRowMaximumHeight(2, 100)
                attr_dendrograms.append(attr_dendrogram)
            
            
            # Heatmap widget for each class 
            for i, (class_, chm) in enumerate(zip(classes, hm)): 
                hm_widget = GraphicsHeatmapWidget(heatmap=chm, parent=widget)
                hm_widget.set_cell_size(int(self.CellWidth), int(self.CellHeight))
                hm_widget.set_cuts(cut_low, cut_high)
                hm_widget.set_color_table(palette)
                hm_widget.set_show_averages(self.ShowAverageStripe)
                hm_widget.cell_tool_tip = lambda row, col, hm=hm_widget: self.cell_tool_tip(hm, row, col)
                layout.addItem(hm_widget, i*2 + 5, column + 2)
                column_heatmap_widgets.append(hm_widget)
            heatmap_widgets.append(column_heatmap_widgets)
            
            # Top attr annotations
            list = GraphicsSimpleTextList(labels, parent=widget, orientation=Qt.Horizontal)
            list.setAlignment(Qt.AlignBottom | Qt.AlignLeft)
            
            list.setFont(font)
            list.layout().setContentsMargins(stripe_offset, 0, 0, 0)
            list.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            
            layout.addItem(list, 3, column + 2, Qt.AlignBottom | Qt.AlignLeft)
            attr_annotation_widgets.append(list)
            attr_annotation_widgets_top.append(list)
            
            # Bottom attr annotations
            list = GraphicsSimpleTextList(labels, parent=widget, orientation=Qt.Horizontal)
            list.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
            
            list.setFont(font)
            list.layout().setContentsMargins(stripe_offset, 0, 0, 0)
            list.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            
            layout.addItem(list, len(hm)*2 + 5, column + 2)
            attr_annotation_widgets.append(list)
            attr_annotation_widgets_bottom.append(list)
            
            # Legend
            if column == 0:
                item = GraphicsLegendWidget(self.heatmapconstructor[0], self.lowerBound, self.upperBound, parent=widget)
                item.set_color_table(palette)
                item.setVisible(self.ShowLegend)
                layout.addItem(item, 0, 2, 1, len(self.heatmaps) + 1)
                layout.setRowSpacing(0, 2)
#                layout.setRowMaximumHeight(0, item.geometry().height())
            
        self.heatmap_scene.addItem(widget)
        self.heatmap_scene.widget = widget
        self.heatmap_widget_grid = heatmap_widgets
        self.gene_annotation_widgets = gene_annotation_widgets
        self.attr_annotation_widgets = attr_annotation_widgets
        self.attr_annotation_widgets_top = attr_annotation_widgets_top
        self.attr_annotation_widgets_bottom = attr_annotation_widgets_bottom
        self.attr_dendrogram_widgets = attr_dendrograms
        
        self.update_annotations()
        self.update_column_annotations()
        
        self.fix_grid_layout()
        
        self.selection_manager.set_heatmap_widgets(heatmap_widgets)
        
    def fix_grid_layout(self):
        """ Fix grid layout when cell size changes or average
        stripes are shown/hiddens.
        """
        if self.heatmap_scene.widget:
            layout = self.heatmap_scene.widget.layout()
            layout.invalidate()
            layout.activate()
            
            for i, hw in enumerate(self.heatmap_widget_grid[0]):
                max_h = hw.size().height()
                layout.setRowMaximumHeight(i*2 + 5, max_h)
                
            for i, hw in enumerate(self.heatmap_widget_grid):
                max_w = hw[0].size().width()
                layout.setColumnMaximumWidth(i + 2, max_w)
                if self.attr_dendrogram_widgets:
                    dendrogram = self.attr_dendrogram_widgets[i]
#                    dendrogram.resize(max_w, -1)
                    dendrogram.setMaximumWidth(max_w)
                self.attr_annotation_widgets_top[i].setMaximumWidth(max_w)
                self.attr_annotation_widgets_bottom[i].setMaximumWidth(max_w)
                
#            for i, (hw, dend) in enumerate(zip(self.heatmap_widget_grid, self.attr_dendrogram_widgets)):
#                max_w = hw[0].size().width()
            
#            self.update_widget_margins()    
            self.heatmap_scene.widget.resize(self.heatmap_scene.widget.sizeHint(Qt.PreferredSize))
            
            self.on_selection_changed()
            self.update_scene_rect()
        
    def update_scene_rect(self):
        rect = QRectF()
        for item in self.heatmap_scene.items():
            rect |= item.sceneBoundingRect()
        self.heatmap_scene.setSceneRect(rect)
            
    def heatmap_widgets(self):
        """ Iterate over heatmap widgets.
        """
        for item in self.heatmap_scene.items():
            if isinstance(item, GraphicsHeatmapWidget):
                yield item
                
    def label_widgets(self):
        """ Iterate over GraphicsSimpleTextList widgets.
        """
        for item in self.heatmap_scene.items():
            if isinstance(item, GraphicsSimpleTextList):
                yield item
                
    def dendrogram_widgets(self):
        """ Iterate over dendrogram widgets
        """
        for item in self.heatmap_scene.items():
            if isinstance(item, DendrogramWidget):
                yield item
                
    def legend_widgets(self):
        for item in self.heatmap_scene.items():
            if isinstance(item, GraphicsLegendWidget):
                yield item
                
    def update_cell_size(self):
        """ Update cell sizes (by user request - height/width sliders)
        """ 
        for heatmap in self.heatmap_widgets():
            heatmap.set_cell_size(self.CellWidth, self.CellHeight)
            
        hor_font = QFont(self.font())
        hor_font.setPointSize(self.point_size_hint(self.CellWidth))
        vert_font = QFont(self.font())
        vert_font.setPointSize(self.point_size_hint(self.CellHeight))
        
        # Also update the annotation items font.
        for labels in self.label_widgets():
            if labels.orientation == Qt.Vertical:
                labels.setFont(vert_font)
            else:
                labels.setFont(hor_font)

        self.update_widget_margins()
        ## To hide the annotations if cell sizes to small.
        self.update_annotations()
        self.update_column_annotations()
            
        self.fix_grid_layout()
        
    def update_widget_margins(self):
        """ Update dendrogram and text list widgets margins to incude the
        space for average stripe.
        """
        if self.ShowAverageStripe:
            stripe_offset = c_averageStripeWidth + 2
        else:
            stripe_offset = 0
        right = self.CellWidth / 2.0
        
        top = self.CellHeight / 2.0 / self.Merge
        bottom = self.CellHeight / 2.0 / self.Merge
        
        for dendrogram in self.dendrogram_widgets():
            layout = dendrogram.layout() 
            if dendrogram.orientation == Qt.Horizontal:
#                index = self.attr_dendrogram_widgets.index(dendrogram)
#                heatmap = self.heatmap_widget_grid[index][0]
#                h_w = heatmap.size().width()
#                d_w = dendrogram.size().width()
#                right_ = d_w - stripe_offset - self.CellWidth - h_w
                _, top_h, _, bottom_h = layout.getContentsMargins()
                layout.setContentsMargins(stripe_offset + self.CellWidth / 2.0, top_h, right, bottom_h)
            else:
                left_v, _, right_v, _ = layout.getContentsMargins()
                layout.setContentsMargins(left_v, top, right_v, bottom)
                
        for widget in self.label_widgets():
            layout = widget.layout()
            if widget.orientation == Qt.Horizontal:
                left_h, top, right, bottom = layout.getContentsMargins()
                layout.setContentsMargins(stripe_offset, top, right, bottom)
        
    def update_averages_stripe(self):
        """ Update the visibility of the averages stripe.
        """
        if self.data:
            for widget in self.heatmap_widgets():
                widget.set_show_averages(self.ShowAverageStripe)
                
            self.update_widget_margins()
            self.fix_grid_layout()
            
    def update_grid_spacing(self):
        """ Update layout spacing.
        """
        if self.scene.widget:
            layout = self.scene.widget.layout()
            layout.setSpacing(self.SpaceX)
            self.fix_grid_layout()
        
    def update_color_schema(self):
        palette = self.getGammaCorrectedPalette() if self.Gamma !=0 else self.palette
        for heatmap in self.heatmap_widgets():
            heatmap.set_color_table(palette)
            
        for legend in self.legend_widgets():
            legend.set_color_table(palette)
            
    def update_thresholds(self):
        self.sliderCutLow.box.setDisabled(not self.CutEnabled)
        self.sliderCutHigh.box.setDisabled(not self.CutEnabled)
            
        if self.data:
            if self.CutEnabled:
                low, high = self.CutLow, self.CutHigh
            else:
                low, high = self.lowerBound, self.upperBound
            for heatmap in self.heatmap_widgets():
                heatmap.set_cuts(low, high)
    
    def update_sorting(self):
        if self.data:
            self.update_heatmaps()
        
    def update_sorting_examples(self):
        if self.data:
            self.update_heatmaps()

    def update_sorting_attributes(self):
        if self.data:
            self.update_heatmaps()

    def update_legend(self):
        for item in self.heatmap_scene.items():
            if isinstance(item, GraphicsLegendWidget):
                item.setVisible(self.ShowLegend)
        
    def update_annotations(self):
        if self.data:
            if self.meta:
                attr = self.meta[self.BAnnotationIndx]
            else:
                attr = None
            
            show = self.ShowGeneAnnotations and attr and self.Merge == 1
            show = show and self.CellHeight > 3
            for list_widget, hm in zip(self.gene_annotation_widgets, self.heatmap_widget_grid[0]):
                list_widget.setVisible(bool(show))
                if show:
                    hm = hm.heatmap
                    examples = hm.examples
                    indices = hm.exampleIndices[:-1]
                    labels = [str(examples[i][attr]) for i in indices]
                    list_widget.set_labels(labels)

    def update_column_annotations(self):
        if self.data:
            show = self.CellWidth > 3
            show_top = self.ShowColumnLabels and self.ColumnLabelPosition == 0 and show
            show_bottom = self.ShowColumnLabels and self.ColumnLabelPosition == 1 and show
            
            for list_widget in self.attr_annotation_widgets_top:
                list_widget.setVisible(show_top)
                
            layout = self.heatmap_scene.widget.layout()
            layout.setRowMaximumHeight(3,  -1 if show_top else 0)
            layout.setRowSpacing(3, -1 if show_top else 0)
                
            for list_widget in self.attr_annotation_widgets_bottom:
                list_widget.setVisible(show_bottom)
                
            layout.setRowMaximumHeight(len(self.heatmap_widget_grid[0]) + 4, -1 if show_top else 0)
                
            self.fix_grid_layout()
            
    def update_heatmaps(self):
        if self.data:
            self.construct_heatmaps(self.data, self.selected_split_label())
            self.create_heatmaps(self.heatmapconstructor)
            self.clear_scene()
            self.construct_heatmaps_scene(self.heatmaps, self.data,
                                          attr_cluster=self.attr_cluster,
                                          data_clusters=self.data_clusters)
        else:
            self.clear()
        
    def update_heatmaps_stage2(self):
        if self.data:
            self.create_heatmaps(self.heatmapconstructor)
            self.clear_scene()
            self.construct_heatmaps_scene(self.heatmaps, self.data,
                                          attr_cluster=self.attr_cluster,
                                          data_clusters=self.data_clusters)
        
    def cell_tool_tip(self, heatmap_widget, row, column):
        if not self.GShowToolTip:
            return ""
        hm = heatmap_widget.heatmap
        examples = hm.examples[hm.exampleIndices[row] : hm.exampleIndices[row+1]]
        domain = hm.examples.domain
        if hm.getCellIntensity(row, column) != None:
            head = "%6.4f" % hm.getCellIntensity(row, column)
        else:
            head = "Missing Data"
        if self.BShowColumnID:
            head += "\n" + domain.attributes[column].name
        # tool tip, construct body
        body = ""
        if (self.BShowSpotIndex and self.BSpotVar) or \
                (self.BShowAnnotation and self.BAnnotationVar) or \
                 self.BShowGeneExpression:
            for (i, e) in enumerate(examples):
                if i > 5:
                    body += "\n... (%d more)" % (len(examples) - 5)
                    break
                else:
                    s = []
                    if self.BShowSpotIndex and self.BSpotVar:
                        s.append(str(e[self.BSpotVar]))
                    if self.BShowGeneExpression:
                        s.append(str(e[column]))
                    if self.BShowAnnotation and self.BAnnotationVar:
                        s.append(str(e[self.BAnnotationVar]))
            
                body += "\n"
                body += " | ".join(s)
        return head + body
    
    def on_merge_changed(self):
        self.oldMerge = self.savedMerge
        if self.MaintainArrayHeight and self.oldMerge != self.Merge:
            k = float(self.Merge) / self.oldMerge
            l = max(1, min(int(self.CellHeight * k), self.maxVSize))
            if l != self.CellHeight:
                self.CellHeight = l
                self.sliderVSize.setValue(self.CellHeight)

        self.update_heatmaps_stage2()
        self.savedMerge = self.Merge
            
    def on_selection_changed(self):
        for item in self.selection_rects:
            item.hide()
            item.setParentItem(None)
            item.update()
            self.heatmap_scene.removeItem(item)
        self.selection_rects = []
        self.selection_manager.update_selection_rects()
        rects = self.selection_manager.selection_rects
        for rect in rects:
            item = QGraphicsRectItem(rect, None, self.heatmap_scene)
            item.setPen(QPen(Qt.black, 2))
            self.selection_rects.append(item)
            
    def on_selection_finished(self):
        self.selected_rows = self.selection_manager.selections
        self.commit_if()
        
    def commit_if(self):
        if self.auto_commit:
            self.commit()
        else:
            self.selection_changed_flag = True
        
    def commit(self):
        data = None
        if self.sorted_data:
            if self.selected_rows:
                
                #obtain examples directly from the heatmap, so their order does not matter

                rd = self.selection_manager.rows_to_heatmaps()
                hr = self.selection_manager._heatmap_ranges

                examples = []
                for row in self.selected_rows:
                    h = rd[row][0]
                    begin,_ = hr[h]
                    hm = h.heatmap
                    examples.extend(hm.examples[hm.exampleIndices[row-begin] : hm.exampleIndices[row+1-begin]])
                    
                data = orange.ExampleTable(examples, name=self.data.name)

            else:
                data = None
        
        self.send("Examples", data)
        self.selection_changed_flag
            
    ## handle saved selections 
    def settingsFromWidgetCallbackSelection(self, handler, context):
        context.selection = self.selection_manager.selections

    def settingsToWidgetCallbackSelection(self, handler, context):
        selection = getattr(context, "selection", None)
        if selection:
            self.selection_manager.select_rows(selection)
            self.selected_rows = selection

    def split_changed(self):
        if self.data:
            self.clear_scene()
            self.construct_heatmaps(self.data, self.selected_split_label())
            self.create_heatmaps(self.heatmapconstructor)
            self.construct_heatmaps_scene(self.heatmaps, self.data,
                                          attr_cluster=self.attr_cluster,
                                          data_clusters=self.data_clusters)
        

class GraphicsPixmapLayoutItem(QGraphicsLayoutItem):
    """ A layout item wraping a QGraphicsPixmapItem
    """
    def __init__(self, pixmap_item, parent=None):
        QGraphicsLayoutItem.__init__(self, parent)
        self.pixmap_item = pixmap_item
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
    def setGeometry(self, rect):
        QGraphicsLayoutItem.setGeometry(self, rect)
        self.pixmap_item.setPos(rect.topLeft())
        
    def sizeHint(self, which, constraint=QSizeF()):
        return QSizeF(self.pixmap_item.pixmap().size())
    
    def setPixmap(self, pixmap):
        self.pixmap_item.setPixmap(pixmap)
        self.updateGeometry()
        
        
class GraphicsHeatmapWidget(QGraphicsWidget):
    def __init__(self, heatmap=None, parent=None, scene=None):
        QGraphicsWidget.__init__(self, parent)
        self.setAcceptHoverEvents(True)
        layout = QGraphicsLinearLayout(Qt.Horizontal)
        layout.setContentsMargins(0, 0, 0, 0)
        item = QGraphicsPixmapItem(self)
        item.setShapeMode(QGraphicsPixmapItem.BoundingRectShape)
        self.heatmap_item = GraphicsPixmapLayoutItem(item, self)
        
        item = QGraphicsPixmapItem(self)
        item.setShapeMode(QGraphicsPixmapItem.BoundingRectShape)
        self.averages_item = GraphicsPixmapLayoutItem(item, self)
        
        layout.addItem(self.averages_item)
        layout.addItem(self.heatmap_item)
        layout.setItemSpacing(0, 2)
        
        self.setLayout(layout)
        
        self.heatmap = None
        self.show_averages = True
        self._pixmap_args = None
        self.color_table = None
        self.selection_manager = None
        self.set_cell_size(4, 4)
        self.set_cuts(0, 255)
        self.set_heatmap(heatmap)
        
        if scene is not None:
            scene.addItem(self)
            
    def clear(self):
        """ Clear the current heatmap.
        """
        self.heatmap = None
        self._pixmap_args = None
        self.heatmap_item.setPixmap(QPixmap())
        self.averages_item.setPixmap(QPixmap())
        self.show_averages = True
        self.layout().invalidate()
            
    def set_heatmap(self, heatmap):
        """ Set the heatmap for display.
        """
        self.clear()
        self.heatmap = heatmap
        self.update()
        
    def set_cell_size(self, width, height):
        self.cell_width = width
        self.cell_height = height
        self.update()
        
    def set_cuts(self, low, high):
        self.cut_low = low
        self.cut_high = high
        self.update()
        
    def set_show_averages(self, show):
        self.show_averages = show
        self._pixmap_args = None
        self.update()

    def set_color_table(self, color_table):
        if qVersion() <= "4.5":
            self.color_table = signedPalette(color_table)
        else:
            self.color_table = color_table
        
        self._pixmap_args = None
        self.update()
            
    def _update_pixmap(self):
        """ Update the pixmap if its construction arguments changed.
        """
        if self.heatmap:
            args = (int(self.cell_width), int(self.cell_height),
                    self.cut_low, self.cut_high, 1.0)
            
            if args != self._pixmap_args:
                bitmap, width, height = self.heatmap.getBitmap(*args)
                image = QImage(bitmap, width, height, QImage.Format_Indexed8)
                color_table = self.color_table
                if color_table:
                    image.setColorTable(color_table)
                self.pixmap = QPixmap.fromImage(image)
                
                bitmap, width, height = self.heatmap.getAverages(*((c_averageStripeWidth,) + args[1:]))
                image = QImage(bitmap, width, height, QImage.Format_Indexed8)
                if self.color_table:
                    image.setColorTable(color_table)
                self.averages_pixmap = QPixmap.fromImage(image)
                
                self._pixmap_args = args
                
                self.layout().invalidate()
        else:
            self.averages_pixmap = None
            self.pixmap = None
            
        self.heatmap_item.setPixmap(self.pixmap or QPixmap())
        if self.show_averages and self.averages_pixmap:
            self.averages_item.setPixmap(self.averages_pixmap)
        else:
            self.averages_item.setPixmap(QPixmap())
            
    def update(self):
        self._update_pixmap()
        QGraphicsWidget.update(self)
        
    def set_selection_manager(self, manager):
        self.selection_manager = manager
        
    def cell_at(self, pos):
        """ Return the cell row, column from a point `pos` in local
        coordinates.
        
        """
        pos = self.mapToItem(self.heatmap_item.pixmap_item, pos)
        x, y = pos.x(), pos.y()
        def clamp(i, m):
            return int(min(max(i, 0), m-1))
        return (clamp(math.floor(y / self.cell_height), self.heatmap.height),
                clamp(math.floor(x / self.cell_width), self.heatmap.width))
    
    def cell_rect(self, row, column):
        """ Return a QRectF in local coordinates containing the cell
        at `row` and `column`.
        
        """
        top = QPointF(column * self.cell_width, row * self.cell_height)
        top = self.mapFromItem(self.heatmap_item.pixmap_item, top)
        size = QSizeF(self.cell_width, self.cell_height)
        return QRectF(top, size)

    def row_rect(self, row):
        """ Return a QRectF in local coordinates containing the entire row.
        """
        rect = self.cell_rect(row, 0).united(self.cell_rect(row, self.heatmap.width - 1))
        rect.setLeft(0) # To include the average stripe if show.
        return rect
    
    def cell_tool_tip(self, row, column):
        hm = self.heatmap
        start = int(hm.exampleIndices[row])
        end = int(hm.exampleIndices[row + 1])
        examples = [hm.examples[start]]
        attr = hm.examples.domain[column]
        val = "%i, %i: %f" % (row, column, float(examples[0][attr]))
        return val
    
    def hoverEnterEvent(self, event):
        row, col = self.cell_at(event.pos())
    
    def hoverMoveEvent(self, event):
        pos = event.pos()
        row, column = self.cell_at(pos)
        tooltip = self.cell_tool_tip(row, column)
        QToolTip.showText(event.screenPos(), tooltip)
        return QGraphicsWidget.hoverMoveEvent(self, event)
    
    def hoverLeaveEvent(self, event):
        row, col = self.cell_at(event.pos())
    
    if DEBUG:
        def paint(self, painter, option, widget=0):
            rect =  self.geometry()
            rect.translate(-self.pos())
            painter.drawRect(rect.adjusted(-1, -1, 1, 1))
    
    
class GridWidget(QGraphicsWidget):
    def __init__(self, parent=None):
        QGraphicsWidget.__init__(self, parent)
        
    if DEBUG:
        def paint(self, painter, option, widget=0):
            rect =  self.geometry()
            rect.translate(-self.pos())
            painter.drawRect(rect)
            
            
class HeatmapScene(QGraphicsScene):
    """ A Graphics Scene with heatmap widgets.
    """
    def __init__(self, parent=None):
        QGraphicsScene.__init__(self, parent)
        self.selection_manager = HeatmapSelectionManager()
        
    def set_selection_manager(self, manager):
        self.selection_manager = manager
        
    def _items(self, pos=None, cls=object):
        if pos is not None:
            items = self.items(QRectF(pos, QSizeF(3, 3)).translated(-1.5, -1.5))
        else:
            items = self.items()
            
        for item in items:
            if isinstance(item, cls):
                yield item
            
    def heatmap_at_pos(self, pos):
        items  = list(self._items(pos, GraphicsHeatmapWidget))
        if items:
            return items[0]
        else:
            return None
        
    def dendrogram_at_pos(self, pos):
        items  = list(self._items(pos, DendrogramItem))
        if items:
            return items[0]
        else:
            return None
        
    def heatmap_widgets(self):
        return self._items(None, GraphicsHeatmapWidget)
        
    def select_from_dendrogram(self, dendrogram, clear=True):
        """ Select all heatmap rows which belong to the dendrogram.
        """
        dendrogram_widget = dendrogram.parentWidget()
        anchors = list(dendrogram_widget.leaf_anchors())
        cluster = dendrogram.cluster
        start, end = anchors[cluster.first], anchors[cluster.last - 1]
        start, end = dendrogram_widget.mapToScene(start), dendrogram_widget.mapToScene(end)
        # Find a heatmap widget containing start and end y coordinates.
        
        heatmap = None
        for hm in self.heatmap_widgets():
            b_rect = hm.sceneBoundingRect()
            if b_rect.contains(QPointF(b_rect.center().x(), start.y())):
                heatmap = hm
                break
            
        if dendrogram:
            b_rect = hm.boundingRect()
            start, end = hm.mapFromScene(start), hm.mapFromScene(end)
            start, _ = hm.cell_at(QPointF(b_rect.center().x(), start.y()))
            end, _ = hm.cell_at(QPointF(b_rect.center().x(), end.y()))
            self.selection_manager.selection_add(start, end, hm, clear=clear)
        return
        
    def mousePressEvent(self, event):
        pos = event.scenePos()
        heatmap = self.heatmap_at_pos(pos)
        if heatmap and event.button() & Qt.LeftButton:
            row, _ = heatmap.cell_at(heatmap.mapFromScene(pos))
            self.selection_manager.selection_start(heatmap, event)
            
        dendrogram = self.dendrogram_at_pos(pos)
        if dendrogram and event.button() & Qt.LeftButton:
            if dendrogram.orientation == Qt.Vertical:
                self.select_from_dendrogram(dendrogram, clear=not event.modifiers() & Qt.ControlModifier)
            return 
        
        return QGraphicsScene.mousePressEvent(self, event)
    
    def mouseMoveEvent(self, event):
        pos = event.scenePos()
        heatmap = self.heatmap_at_pos(pos)
        if heatmap and event.buttons() & Qt.LeftButton:
            row, _ = heatmap.cell_at(heatmap.mapFromScene(pos))
            self.selection_manager.selection_update(heatmap, event)
            
        dendrogram = self.dendrogram_at_pos(pos)
        if dendrogram and dendrogram.orientation == Qt.Horizontal: # Filter mouse move events
            return
            
        return QGraphicsScene.mouseMoveEvent(self, event)
    
    def mouseReleaseEvent(self, event):
        pos = event.scenePos()
        heatmap = self.heatmap_at_pos(pos)
        if heatmap:
            row, _ = heatmap.cell_at(heatmap.mapFromScene(pos))
            self.selection_manager.selection_finish(heatmap, event)
        
        dendrogram = self.dendrogram_at_pos(pos)
        if dendrogram and dendrogram.orientation == Qt.Horizontal: # Filter mouse events
            return
        
        return QGraphicsScene.mouseReleaseEvent(self, event)
    
    def mouseDoubleClickEvent(self, event):
        pos = event.scenePos()
        dendrogram = self.dendrogram_at_pos(pos)
        if dendrogram: # Filter mouse events
            return
        return QGraphicsScene.mouseDoubleClickEvent(self, event)
        
        
class GtI(QGraphicsSimpleTextItem):
    if DEBUG:
        def paint(self, painter, option, widget =0):
            QGraphicsSimpleTextItem.paint(self, painter, option, widget)
            painter.drawRect(self.boundingRect())
            
    
class GraphicsSimpleTextLayoutItem(QGraphicsLayoutItem):
    """ A Graphics layout item wrapping a QGraphicsSimpleTextItem alowing it 
    to be managed by a layout.
    """
    def __init__(self, text_item, orientation=Qt.Horizontal, parent=None):
        QGraphicsLayoutItem.__init__(self, parent)
        self.orientation = orientation
        self.text_item = text_item
        if orientation == Qt.Vertical:
            self.text_item.rotate(-90)
            self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        else:
            self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        
    def setGeometry(self, rect):
        QGraphicsLayoutItem.setGeometry(self, rect)
        if self.orientation == Qt.Horizontal:
            self.text_item.setPos(rect.topLeft())
        else:
            self.text_item.setPos(rect.bottomLeft())
        
    def sizeHint(self, which, constraint=QSizeF()):
        if which in [Qt.PreferredSize]:
            size = self.text_item.boundingRect().size()
            if self.orientation == Qt.Horizontal:
                return size
            else:
                return QSizeF(size.height(), size.width())
        else:
            return QSizeF()
    
    def setFont(self, font):
        self.text_item.setFont(font)
        self.updateGeometry()
        
    def setText(self, text):
        self.text_item.setText(text)
        self.updateGeometry()
        
        
class GraphicsSimpleTextList(QGraphicsWidget):
    """ A simple text list widget.
    """
    def __init__(self, labels=[], orientation=Qt.Vertical, parent=None, scene=None):
        QGraphicsWidget.__init__(self, parent)
        layout = QGraphicsLinearLayout(orientation)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.setLayout(layout)
        self.orientation = orientation
        self.alignment = Qt.AlignCenter
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.set_labels(labels)
        
        if scene is not None:
            scene.addItem(self)
        
    def clear(self):
        """ Remove all text items.
        """
        layout = self.layout()
        for i in reversed(range(layout.count())):
            item = layout.itemAt(i)
            item.text_item.setParentItem(None)
            if self.scene():
                self.scene().removeItem(item.text_item)
            layout.removeAt(i)
        
        self.label_items = []
        self.updateGeometry()
        
    def set_labels(self, labels):
        """ Set the text labels to show in the widget.
        """
        self.clear()
        orientation = Qt.Horizontal if self.orientation == Qt.Vertical else Qt.Vertical
        for text in labels:
#            item = QGraphicsSimpleTextItem(text, self)
            item = GtI(text, self)
            item.setFont(self.font())
            item.setToolTip(text)
            item = GraphicsSimpleTextLayoutItem(item, orientation, parent=self)
            self.layout().addItem(item)
            self.layout().setAlignment(item, self.alignment)
            self.label_items.append(item)
            
        self.layout().activate()
        self.updateGeometry()
    
    def setAlignment(self, alignment):
        """ Set alignment of text items in the widget
        """
        self.alignment = alignment
        layout = self.layout()
        for i in range(layout.count()):
            layout.setAlignment(layout.itemAt(i), alignment)
            
    def setVisible(self, bool):
        QGraphicsWidget.setVisible(self, bool)
        self.updateGeometry()
            
    def setFont(self, font):
        """ Set the font for the text.
        """
        QGraphicsWidget.setFont(self, font)
        for item in self.label_items:
            item.setFont(font)
        self.layout().invalidate()
        self.updateGeometry()
        
    def sizeHint(self, which, constraint=QRectF()):
        if not self.isVisible():
            return QSizeF(0, 0)
        else:
            return QGraphicsWidget.sizeHint(self, which, constraint)
            
    if DEBUG:
        def paint(self, painter, options, widget=0):
            rect =  self.geometry()
            rect.translate(-self.pos())
            painter.drawRect(rect)
        
class GraphicsLegendWidget(QGraphicsWidget):
    def __init__(self, heatmap_constructor, low, high, parent=None, scene=None):
        QGraphicsWidget.__init__(self, parent)
        layout = QGraphicsLinearLayout(Qt.Vertical)
        self.setLayout(layout)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(1)
        
        layout_labels = QGraphicsLinearLayout(Qt.Horizontal)
        layout.addItem(layout_labels)
        layout_labels.setContentsMargins(0, 0, 0, 0)
        label_lo = GtI("%.2f" % low, self)
        label_hi = GtI("%.2f" % high, self)
        self.item_low = GraphicsSimpleTextLayoutItem(label_lo, parent=self)
        self.item_high = GraphicsSimpleTextLayoutItem(label_hi, parent=self)
        
        layout_labels.addItem(self.item_low)
        layout.addStretch()
        layout_labels.addItem(self.item_high)
        
        self._pixmap = QPixmap(c_legendHeight, c_legendHeight)
        self.pixmap_item = QGraphicsPixmapItem(self._pixmap, self)
        self.pixmap_item = GraphicsPixmapLayoutItem(self.pixmap_item, parent=self)
        layout.addItem(self.pixmap_item)
        
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.set_legend(heatmap_constructor, low, high)
        
        if scene is not None:
            scene.addItem(self)
        
    def set_legend(self, heatmap_constructor, low, high):
        self.heatmap_constructor = heatmap_constructor
        self.low = low
        self.high = high
        self.color_table = None
        self._pixmap = None
        self._pixmap_args = None
        self.update()
        self.updateGeometry()
        
    def set_color_table(self, color_table):
        if qVersion() <= "4.5":
            self.color_table = signedPalette(color_table)
        else:
            self.color_table = color_table
        
        self._pixmap_args = None
        self.update()
        
    def update(self):
        crect = self.contentsRect()
        width = crect.width()
        height = c_legendHeight
        if not self.pixmap_item or self._pixmap_args != (width, height):
            bitmap = self.heatmap_constructor.getLegend(int(width), int(height), 1.0)
            image = QImage(bitmap, width, height, QImage.Format_Indexed8)
            color_table = self.color_table
            if color_table:
                image.setColorTable(color_table)
            self._pixmap = QPixmap.fromImage(image)
            
        self.pixmap_item.setPixmap(self._pixmap)
        self.item_low.setText("%.2f" % self.low)
        self.item_high.setText("%.2f" % self.high)
        self.layout().activate()
        QGraphicsWidget.update(self)
            
    if DEBUG:
        def paint(self, painter, options, widget=0):
            rect =  self.geometry()
            rect.translate(-self.pos())
            painter.drawRect(rect)

        
class HeatmapSelectionManager(QObject):
    """ Selection manager for heatmap rows
    """
    def __init__(self, parent=None):
        QObject.__init__(self, parent)
        self.selections = []
        self.selection_ranges = []
        self.heatmap_widgets = []
        self.selection_rects = []
        self.heatmaps = []
        self._heatmap_ranges = {}
        self._start_row = 0
        
    def set_heatmap_widgets(self, widgets):
        self.remove_rows(self.selections)
        self.heatmaps = widgets
        
        # Compute row ranges for all heatmaps
        self._heatmap_ranges = {}
        for group in widgets:
            start = end = 0
            for heatmap in group:
                end += heatmap.heatmap.height
                self._heatmap_ranges[heatmap] = (start, end)
                start = end
        
    def select_rows(self, rows, heatmap=None, clear=True):
        """ Add `rows` to selection. If `heatmap` is provided the rows
        are mapped from the local indices to global heatmap indics. If `clear`
        then remove previous rows.
        """
        if heatmap is not None:
            start, end = self._heatmap_ranges[heatmap]
            rows = [start + r for r in rows]
            
        old_selection = list(self.selections)
        if clear:
            self.selections = rows
        else:
            self.selections = sorted(set(self.selections + rows))

        if self.selections != old_selection:
            self.update_selection_rects()
            self.emit(SIGNAL("selection_changed()"))
            self.emit(SIGNAL("selection_rects_changed"), self.selection_rects)
            
    def remove_rows(self, rows):
        """ Remove `rows` from the selection.
        """
        old_selection = list(self.selections)
        self.selections = sorted(set(self.selections) - set(rows))
        if old_selection != self.selections:
            self.update_selection_rects()
            self.emit(SIGNAL("selection_changed()"))
            self.emit(SIGNAL("selection_rects_changed"), self.selection_rects)
            
    def combined_ranges(self, ranges):
        combined_ranges = set()
        for start, end in ranges:
            if start <= end:
                rng = range(start, end + 1)
            else:
                rng = range(start, end - 1, -1)
            combined_ranges.update(rng)
        return sorted(combined_ranges)
        
    def selection_start(self, heatmap_widget, event):
        """ Selection  started by `heatmap_widget` due to `event`.
        """
        pos = heatmap_widget.mapFromScene(event.scenePos())
        row, column = heatmap_widget.cell_at(pos)
        start, _ = self._heatmap_ranges[heatmap_widget]
        row = start + row
        self._start_row = row
        range = (row, row)
        if event.modifiers() & Qt.ControlModifier:
            self.selection_ranges.append(range)
        else:
            self.selection_ranges = [range]
        self.select_rows(self.combined_ranges(self.selection_ranges))
        
    def selection_update(self, heatmap_widget, event):
        """ Selection updated by `heatmap_widget due to `event` (mouse drag).
        """
        pos = heatmap_widget.mapFromScene(event.scenePos())
        row, column = heatmap_widget.cell_at(pos)
        start, _ = self._heatmap_ranges[heatmap_widget]
        row = start + row
        if self.selection_ranges:
            self.selection_ranges[-1] = (self._start_row, row)
        else:
            self.selection_ranges = [(row, row)]
            
        self.select_rows(self.combined_ranges(self.selection_ranges))
        
    def selection_finish(self, heatmap_widget, event):
        """ Selection finished by `heatmap_widget due to `event`.
        """
        pos = heatmap_widget.mapFromScene(event.scenePos())
        row, column = heatmap_widget.cell_at(pos)
        start, _ = self._heatmap_ranges[heatmap_widget]
        row = start + row
        range = (self._start_row, row)
        self.selection_ranges[-1] = range
        self.select_rows(self.combined_ranges(self.selection_ranges),
                         clear=not event.modifiers() & Qt.ControlModifier)
        self.emit(SIGNAL("selection_finished()"))
        
    def selection_add(self, start, end, heatmap=None, clear=True):
        """ Add a selection range from `start` to `end`.
        """ 
        if heatmap is not None:
            _start, _ = self._heatmap_ranges[heatmap]
            start = _start + start
            end = _start + end
        
        if clear:
            self.selection_ranges = []
        self.selection_ranges.append((start, end))
        self.select_rows(self.combined_ranges(self.selection_ranges))
        self.emit(SIGNAL("selection_finished()"))
        
    def update_selection_rects(self):
        """ Update the selection rects.
        """
        def continuous_ranges(selections):
            """ Group continuous ranges
            """
            selections = iter(selections)
            start = end = selections.next()
            try:
                while True:
                    new_end = selections.next()
                    if new_end > end + 1:
                        yield start, end
                        start = end = new_end
                    else:
                        end = new_end
            except StopIteration:
                yield start, end
                
        def group_selections(selections):
            """ Group selections along with heatmaps.
            """
            rows2hm = self.rows_to_heatmaps()
            selections = iter(selections)
            start = end = selections.next()
            end_heatmaps = rows2hm[end]
            try:
                while True:
                    new_end = selections.next()
                    new_end_heatmaps = rows2hm[new_end]
                    if new_end > end + 1 or new_end_heatmaps != end_heatmaps:
                        yield start, end, end_heatmaps
                        start = end = new_end
                        end_heatmaps = new_end_heatmaps
                    else:
                        end = new_end
                        
            except StopIteration:
                yield start, end, end_heatmaps
                
        def selection_rect(start, end, heatmaps):
            rect = QRectF()
            for heatmap in heatmaps:
                h_start, _ = self._heatmap_ranges[heatmap] 
                rect |= heatmap.mapToScene(heatmap.row_rect(start - h_start)).boundingRect()
                rect |= heatmap.mapToScene(heatmap.row_rect(end - h_start)).boundingRect()
            return rect
             
        self.selection_rects = []
        for start, end, heatmaps in group_selections(self.selections):
            rect = selection_rect(start, end, heatmaps)
            self.selection_rects.append(rect)
            
    def rows_to_heatmaps(self):
        heatmap_groups = zip(*self.heatmaps)
        rows2hm = {}
        for heatmaps in heatmap_groups:
            hm = heatmaps[0]
            start, end = self._heatmap_ranges[hm]
            rows2hm.update(dict.fromkeys(range(start, end), heatmaps))
        return rows2hm
    
         
##################################################################################################
# test script

if __name__=="__main__":
    a=QApplication(sys.argv)
    ow = OWHeatMap()
    
    ow.set_dataset(orange.ExampleTable("brown-selected"), 0)
    ow.handleNewSignals()
    ow.show()
    
    a.exec_()
    ow.saveSettings()

