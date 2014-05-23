"""
<name>Vulcano Plot</name>
<description>Plots fold change vs. p-value.)</description>
<priority>1020</priority>
<contact>Ales Erjavec (ales.erjavec@fri.uni-lj.si)</contact>
<icon>icons/VulcanoPlot.svg</icon>
"""

from __future__ import absolute_import

from operator import add
from collections import defaultdict

import numpy

import orange

from Orange.orng.orngDataCaching import data_hints
from Orange.OrangeWidgets.OWToolbars import ZoomSelectToolbar
from Orange.OrangeWidgets.OWWidget import *
from Orange.OrangeWidgets.OWGraph import *
from Orange.OrangeWidgets import OWGUI

from .. import obiExpression

NAME = "Vulcano Plot"
DESCRIPTION = "Plots fold change vs. p-value.)"
ICON = "icons/VulcanoPlot.svg"
PRIORITY = 1020

INPUTS = [("Examples", Orange.data.Table, "setData")]
OUTPUTS = [("Examples with selected attributes", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWVulcanoPlot.OWVulcanoPlot"]


class GraphSelections(QObject):
    """ Selection manager using a union of rectangle areas 
    """
    def __init__(self, parent):
        QObject.__init__(self, parent)
        self.selection = []
        
    def getPos(self, event):
        graph = self.parent()
        pos = graph.canvas().mapFrom(graph, event.pos())
        x = graph.invTransform(QwtPlot.xBottom, pos.x())
        y = graph.invTransform(QwtPlot.yLeft, pos.y())
        return QPointF(x, y)
        
    def start(self, event):
        pos = self.getPos(event)
        if event.modifiers() & Qt.ControlModifier:
            self.selection.append((pos, pos))
        else:
            self.selection = [(pos, pos)]
        self.emit(SIGNAL("selectionGeometryChanged()"))
    
    def update(self, event):
        pos = self.getPos(event)
        self.selection[-1] = self.selection[-1][:-1] + (pos,)
        self.emit(SIGNAL("selectionGeometryChanged()"))
    
    def end(self, event):
        self.update(event)
        
    def testSelection(self, data):
        if len(data) == 0:
            return []
        data = numpy.asarray(data)
        region = QPainterPath()
        for p1, p2 in self.selection:
            region.addRect(QRectF(p1, p2).normalized())
        def test(point):
            return region.contains(QPointF(point[0], point[1]))
        test = numpy.apply_along_axis(test, 1, data)
        return test
        
class SymetricSelections(GraphSelections):
    """ Selection manager using two symmetric areas extending to 'infinity'  
    """
    def __init__(self, parent, x=3, y=3):
        GraphSelections.__init__(self, parent)
        max = 100000
        self.selection = [(QPointF(-max, max), QPointF(-x, y)), (QPointF(max, max), QPointF(x, y))]
        self.updateAxes = None
        
    def updateSelection(self, axes, pos):
        if axes == QwtPlot.xBottom or axes == -1:
            self.selection[0][1].setX(-abs(pos.x()))
            self.selection[1][1].setX(abs(pos.x()))
        if axes == QwtPlot.yLeft or axes == -1:
            self.selection[0][1].setY(pos.y())
            self.selection[1][1].setY(pos.y())
            
        self.emit(SIGNAL("selectionGeometryChanged()"))
        
    def getAxesAndPos(self, event):
        graph = self.parent()
        pos = graph.canvas().mapFrom(graph, event.pos())
        x = graph.invTransform(QwtPlot.xBottom, pos.x())
        y = graph.invTransform(QwtPlot.yLeft, pos.y())
        
        offset = 3
        dx = abs(graph.invTransform(QwtPlot.xBottom, pos.x() + offset) - x)
        dy = abs(graph.invTransform(QwtPlot.yLeft, pos.y() + offset) - y)
        
        x = abs(x)
        
        cx = self.selection[1][1].x()
        cy = self.selection[1][1].y()

        bottom = QRectF(QPointF(cx, cy), QPointF(graph.maxX, cy)).adjusted(-dx, dy, dx, -dy).normalized()
        left = QRectF(QPointF(cx, graph.maxY), QPointF(cx, cy)).adjusted(-dx, dy, dx, -dy).normalized()
        
        if bottom.contains(QPointF(x, y)) or bottom.contains(QPointF(-x, y)):
            axes = QwtPlot.yLeft
        elif left.contains(QPointF(x, y)) or left.contains(QPointF(-x, y)):
            axes = QwtPlot.xBottom
        else:
            axes = -1
        return axes, QPointF(x, y)
        
    def start(self, event):
        axes, pos = self.getAxesAndPos(event)
        self.updateAxes = axes
        self.updateSelection(axes, pos)
        
    def update(self, event):
        _, pos = self.getAxesAndPos(event)
        self.updateSelection(self.updateAxes, pos)
    
    def end(self, event):
        self.update(event)
        self.updateAxes = None
        
    def testSelection(self, data):
        if len(data) == 0:
            return []
        data = numpy.asarray(data)
        cutoffX = self.selection[1][1].x()
        cutoffY = self.selection[1][1].y()
        return (numpy.abs(data[:, 0]) >= cutoffX) & (data[:, 1] >= cutoffY)

    def cuttoff(self):
        """
        Return the absolute x, y cutoff values.
        """
        return (self.selection[1][1].x(), self.selection[1][1].y())

from Orange.OrangeWidgets.OWItemModels import PyListModel


def item_selection(indices, model, selection=None, column=0):
    """ Create an QItemSelection for indices in model.
    """
    if selection is None:
        selection = QItemSelection()
        
    for i in indices:
        selection.select(model.index(i, column))
    return selection


class LabelSelectionWidget(QWidget):
    """ A widget for selection of label values.
    """
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self._values_model = PyListModel([], parent=self)
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        def group_box(title):
            box = QGroupBox(title)
            box.setFlat(True)
            lay = QVBoxLayout()
            lay.setContentsMargins(0, 0, 0, 0)
            box.setLayout(lay)
            return box
        
        self.labels_combo = QComboBox()
        self.values_view = QListView()
        self.values_view.setSelectionMode(QListView.ExtendedSelection)
        self.values_view.setModel(self._values_model)
        
        
        self.connect(self.labels_combo, SIGNAL("activated(int)"),
                     self.on_label_activated)
        
        self.connect(self.values_view.selectionModel(),
                     SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),
                     self.on_values_selection)
        
        l_box = group_box("Label")
        v_box = group_box("Values")
        
        l_box.layout().addWidget(self.labels_combo)
        v_box.layout().addWidget(self.values_view)
        
        layout.addWidget(l_box)
        layout.addWidget(v_box)
        
        self.setLayout(layout)
        
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        self._block_selection_emit = False
        
        self.labels = []
        self.set_labels([])
        
    def clear(self):
        self.labels_combo.clear()
        self._values_model[:] = []
        self.labels = []
        
    def set_labels(self, labels):
        """ Set the labels to display.
        """
        self.clear()
        if isinstance(labels, dict):
            labels = labels.items()
            
        self.labels = labels
        for label, values in labels:
            self.labels_combo.addItem(label)
            
        if labels:
            self.set_current_label(0)
#            self.labels_combo.setCurrentIndex(0)
            
    def set_selection(self, label, values):
        """ Set the selection to label and values
        """
        if isinstance(label, basestring):
            labels = [l for l, _ in self.labels]
            index = labels.index(label) if label in labels else -1
        else:
            index = label
            
        if index >= 0:
            if index != self.labels_combo.currentIndex():
                self.set_current_label(index)
                
            all_values = list(self._values_model)
            values = [v for v in values if v in all_values]
            selection = QItemSelection()
            for i, v in enumerate(self._values_model):
                if v in values:
                    index = self._values_model.index(i, 0)
                    selection.select(index, index)
            self.values_view.selectionModel().select(selection,  QItemSelectionModel.ClearAndSelect)
        else:
            self.values_view.selectionModel().clear()
            
    def set_current_label(self, index):
        """ Set the current label
        """
        self.labels_combo.setCurrentIndex(index)
        label, values = self.labels[index]
        # Block selection changed
        with self._blocked_signals():
            self._values_model[:] = values
        
    def on_label_activated(self, index):
        label, values = self.labels[index]
        with self._blocked_signals():
            self._values_model[:] = values
        self.emit(SIGNAL("label_activated()"))
        self.emit(SIGNAL("label_activated(int)"), index)
    
    def on_values_selection(self, selected, deselected):
        label, values = self.current_selection()
        self.emit(SIGNAL("selection_changed()"))
        self.emit(SIGNAL("selection_changed(PyQt_PyObject, PyQt_PyObject)"),
                  label, values)
        
    def selection_indexes(self):
        """ Return the values selection indices.
        """
        selection = self.values_view.selectionModel().selection()
        indexes = selection.indexes()
        return sorted(set([i.row() for i in indexes]))        
        
    def current_selection(self):
        """ Return the current label and selected values.
        """
        i = self.labels_combo.currentIndex()

        if i == -1:
            # When clearing the labels model / combobox
            return None, None

        label, all_values = self.labels[i]
        values = [all_values[i] for i in self.selection_indexes()]
        return label, values
    
    def _blocked_signals(self):
        """ Return a context handler blocking all emited signals from this
        object.
         
        """
        class block(object):
            def __enter__(blocker):
                self.blockSignals(True)
            def __exit__(blocker, *args):
                self.blockSignals(False)
                return False
        return block()
                
    def sizeHint(self):
        return QSize(100, 200)


class VulcanoGraph(OWGraph):
    def __init__(self, *args, **kwargs):
        OWGraph.__init__(self, *args, **kwargs)
        # Absolute cutoff values for symmetric selection mode.
        self.cutoffX = 2.0
        self.cutoffY = 2.0
        # maximum absolute x, y values.
        self.maxX, self.maxY = 10, 10
        self.symbolSize = 5
        self.symetricSelections = True

        self.selectedCurve = self.addCurve("", brushColor=Qt.red)
        self.unselectedCurve = self.addCurve("", brushColor=Qt.blue)

        self.plotValues = {}

        self.setAxisAutoScale(QwtPlot.xBottom)
        self.setAxisAutoScale(QwtPlot.yLeft)

        self.setSelection(
            SymetricSelections(
                self,
                x=min(self.maxX, self.cutoffX),
                y=min(self.maxY, self.cutoffY))
        )

    def setSelection(self, selection):
        self.selection = selection
        self.connect(self.selection, SIGNAL("selectionGeometryChanged()"), self.onSelectionChanged)
        if self.plotValues:
            self.updateSelectionArea()

    def onSelectionChanged(self):
        if self.symetricSelections:
            self.cutoffX, self.cutoffY = self.selection.cuttoff()

        self.replot_()
        self.emit(SIGNAL("selectionChanged()"))

    def splitSelected(self):
        test =  self.selection.testSelection(self.plotData)
        return (self.plotData[numpy.nonzero(test)], self.plotData[numpy.nonzero(~test)])

    def setPlotValues(self, values):
        self.plotValues = values
        self.plotData = numpy.array(values.values())
        if self.plotData.size:
            self.maxX = numpy.max(numpy.abs(self.plotData[:, 0]))
            self.maxY = numpy.max(self.plotData[:, 1])
        else:
            self.maxX, self.maxY = 10, 10

        self.setAxisScale(QwtPlot.xBottom, -self.maxX, self.maxX)
        self.setAxisScale(QwtPlot.yLeft, 0.0, self.maxY)
        self.replot_()
        self.emit(SIGNAL("selectionChanged()"))

    def createSelectionRectCurve(self, p1, p2):
        curve = self.addCurve("selection", style=QwtPlotCurve.Lines, penColor=Qt.red, symbol=QwtSymbol.NoSymbol)
        curve.setData([p1.x(), p2.x(), p2.x(), p1.x(), p1.x()], [p1.y(), p1.y(), p2.y(), p2.y(), p1.y()])
        
    def items(self, title=None):
        for item in self.itemList():
            if str(item.title().text()) == title:
                yield item
        
    def updateSelectionArea(self):
        for c in self.items(title="selection"):
            c.detach()
        for p1, p2 in self.selection.selection:
            self.createSelectionRectCurve(p1, p2)

    def replot_(self):
        if self.plotValues:
            selected, unselected = self.splitSelected()
            self.selectedCurve.setData(selected[:, 0], selected[:, 1])
            self.unselectedCurve.setData(unselected[:, 0], unselected[:, 1])
            self.updateSelectionArea()
        else:
            for curve in [self.selectedCurve, self.unselectedCurve]:
                curve.setData([], [])
        self.replot()

    def mousePressEvent(self, event):
        canvas = self.canvas()
        canvasPos = canvas.mapFrom(self, event.pos())
        if canvas.contentsRect().contains(canvasPos) and \
                self.state == SELECT and event.button() == Qt.LeftButton:
            self.selection.start(event)
        else:
            OWGraph.mousePressEvent(self, event)

    def mouseMoveEvent(self, event):
        if self.state == SELECT:
            if event.buttons() & Qt.LeftButton:
                self.selection.update(event)
            if isinstance(self.selection, SymetricSelections):
                axes, pos = self.selection.getAxesAndPos(event)
                cursors = {QwtPlot.xBottom: Qt.SizeHorCursor,
                           QwtPlot.yLeft: Qt.SizeVerCursor}
                self.canvas().setCursor(cursors.get(axes, self._cursor))
        else:
            OWGraph.mouseMoveEvent(self, event)

    def mouseReleaseEvent(self, event):
        if self.state == SELECT:
            if event.button() == Qt.LeftButton:
                self.selection.end(event)
        else:
            OWGraph.mouseReleaseEvent(self, event)

    def reselect(self, replot=True):
        if self.symetricSelections:
            self.setSelection(
                SymetricSelections(
                    self,
                    x=min(self.maxX, self.cutoffX),
                    y=min(self.maxY, self.cutoffY))
            )
        else:
            self.setSelection(GraphSelections(self))
            self.canvas().setCursor(self._cursor)
        if replot:
            self.replot_()

    def updateSymbolSize(self):
        def setSize(curve, size):
            symbol = curve.symbol()
            symbol.setSize(size)
            if QWT_VERSION_STR >= "5.2":
                curve.setSymbol(symbol)
        setSize(self.selectedCurve, self.symbolSize)
        setSize(self.unselectedCurve, self.symbolSize)
        self.replot()


from .OWGenotypeDistances import SetContextHandler
from .OWFeatureSelection import disable_controls


class OWVulcanoPlot(OWWidget):
    settingsList = ["graph.cutoffX", "graph.cutoffY", "graph.symbolSize",
                    "graph.symetricSelections", "showXTitle", "showYTitle"]

    contextHandlers = {
        "targets": SetContextHandler("targets")
    }

    def __init__(self, parent=None, signalManager=None, name="Vulcano Plot"):
        OWWidget.__init__(self, parent, signalManager, name, wantGraph=True)

        self.inputs = [("Examples", ExampleTable, self.setData)]
        self.outputs = [("Examples with selected attributes", ExampleTable)]

        self.genes_in_columns = False

        self.showXTitle = True
        self.showYTitle = True

        self.auto_commit = False
        self.selection_changed_flag = False
        self.target_group = None, []
        self.label_selections = []

        self.graph = VulcanoGraph(self)
        self.connect(self.graph, SIGNAL("selectionChanged()"),
                     self.on_selection_changed)
        self.mainArea.layout().addWidget(self.graph)

        self.loadSettings()

        ## GUI
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.infoLabel = OWGUI.label(box, self, "")
        self.infoLabel.setText("No data on input")
        self.infoLabel2 = OWGUI.label(box, self, "")
        self.infoLabel2.setText("0 selected genes")
        
        box = OWGUI.widgetBox(self.controlArea, "Target Labels")
        
        self.target_widget = LabelSelectionWidget(self)
        self.connect(self.target_widget,
                     SIGNAL("selection_changed(PyQt_PyObject, PyQt_PyObject)"),
                     self.on_target_changed)
        self.connect(self.target_widget,
                     SIGNAL("label_activated(int)"),
                     self.on_label_activated)
        
        box.layout().addWidget(self.target_widget)
        
        self.genesInColumnsCheck = OWGUI.checkBox(box, self, "genes_in_columns",
                                    "Genes in columns", 
                                    callback=[self.init_from_data, self.plot])

        box = OWGUI.widgetBox(self.controlArea, "Settings")
        OWGUI.hSlider(box, self, "graph.symbolSize", label="Symbol size:   ", minValue=2, maxValue=20, step=1, callback = self.graph.updateSymbolSize)
        OWGUI.checkBox(box, self, "showXTitle", "X axis title", callback=self.setAxesTitles)
        OWGUI.checkBox(box, self, "showYTitle", "Y axis title", callback=self.setAxesTitles)
        
        toolbar = ZoomSelectToolbar(self, self.controlArea, self.graph, buttons=[ZoomSelectToolbar.IconSelect, ZoomSelectToolbar.IconZoom, ZoomSelectToolbar.IconPan])
        
        top_layout = toolbar.layout()
        top_layout.setDirection(QBoxLayout.TopToBottom)
        button_layotu = QHBoxLayout()
        top_layout.insertLayout(0, button_layotu)
        
        for i in range(1, top_layout.count()):
            item = top_layout.itemAt(1)
            top_layout.removeItem(item)
            button_layotu.addItem(item)
        
        OWGUI.checkBox(toolbar, self, "graph.symetricSelections", "Symetric selection", callback=self.graph.reselect)

        box = OWGUI.widgetBox(self.controlArea, "Commit")
        b = OWGUI.button(box, self, "Commit", callback=self.commit, default=True)
        cb = OWGUI.checkBox(box, self, "auto_commit", "Commit automatically")
        OWGUI.setStopper(self, b, cb, "selection_changed_flag", self.commit_if)

        self.connect(self.graphButton, SIGNAL("clicked()"), self.graph.saveToFile)
        
        OWGUI.rubber(self.controlArea)

        self.data = None
        self.target_group = None, []
        self.current_selection = []

        self.graph.reselect(True)
        self.resize(800, 600)

    def clear(self):
        self.target_widget.set_labels([])
        self.targets = []
        self.label_selections = []
        self.target_group = None, []
        self.clear_graph()
        
    def clear_graph(self):
        self.values = {}
        self.graph.setPlotValues({})
        self.updateTooltips()

    def setData(self, data=None):
#         self.closeContext("")
        self.closeContext("targets")
        self.clear()
        self.data = data
        self.error(0)
        self.warning([0, 1])
        if data:
            self.genes_in_columns = not bool(data.domain.classVar)
            self.genesInColumnsCheck.setDisabled(not bool(data.domain.classVar))
            if self.genes_in_columns:
                self.genes_in_columns = not data_hints.get_hint(data, "genesinrows", not self.genes_in_columns)
#             self.openContext("", data)
        else:
            self.infoLabel.setText("No data on input.")
        self.init_from_data()

    def init_from_data(self):
        """ Init widget state from the data.
        """
        self.update_target_labels()
        self.error(0)
        if self.data:
            if not self.targets:
                if self.genes_in_columns:
                    self.error(0, "Data set with no column labels (attribute tags)")
                else:
                    self.error(0, "Data has no class.")
        
        self.openContext("targets", [(label, v) for label, vals in self.targets \
                                                for v in vals])
        
        if len(self.label_selections) != len(self.targets): # Some times this happens.
            self.label_selections = [[] for t in self.targets]
            
        if self.target_group == (None, []) and self.targets:
            label, values = self.targets[0]
            self.target_group = (label, values[:1])
            
        if self.target_group != (None, []):
            self.target_widget.set_selection(*self.target_group)
        else:
            self.clear_graph()

    def update_target_labels(self):
        if self.data:
            if self.genes_in_columns:
                items = [a.attributes.items() for a in self.data.domain.attributes]
                items = reduce(add, items, [])
                
                targets = defaultdict(set)
                for label, value in items:
                    targets[label].add(value)
                    
                targets = [(key, list(sorted(vals))) for key, vals in targets.items() \
                           if len(vals) >= 2]
                self.targets = targets
                
            else:
                var = self.data.domain.classVar
                values = list(var.values)
                if len(values) >= 2:
                    self.targets = [(var.name, values)]
                else:
                    self.targets = []
        else:
            self.targets = []
            
        if self.targets:
            label, values = self.targets[0]
            self.target_group = (label, values[:1])
        else:
            self.target_group = None, []
            
        self.label_selections = [[] for t in self.targets]
        self.target_widget.set_labels(self.targets)
                
        
    def on_label_activated(self, index):
        """ Try to restore a saved selection.
        """
        selected = self.label_selections[index]
        if not selected:
            selected = self.targets[index][1][:1]
            
        self.target_widget.set_selection(index, selected)
        
    def on_target_changed(self, label, values):
        self.target_group = label, values
        # Save the selection
        labels = [l for l, _ in self.targets]
        if label in labels:
            index = labels.index(label)
            self.label_selections[index] = values
            
        # replot
        if label and values:
            self.plot()
        else:
            self.clear_graph()
    
    @disable_controls
    def plot(self):
        self.values = {}
        self.current_selection = []
        target_label, target_values = self.target_group
        self.warning([0, 1])
        self.error(1)
        if self.data and target_values:
            target_label, target_values = self.target_group
            if self.genes_in_columns:
                target = set([(target_label, value) for value in target_values])
            else:
                target = set(target_values)
            
            ttest = obiExpression.ExpressionSignificance_TTest(self.data, useAttributeLabels=self.genes_in_columns)
            ind1, ind2 = ttest.test_indices(target)
            
            if not len(ind1) or not len(ind2):
                self.error(1, "Target labels most exclude/include at least one value.")
                
            if len(ind1) < 2 and len(ind2) < 2:
                self.warning(0, "Insufficient data to compute statistics. More than one measurement per class should be provided")
            
            self.progressBarInit()
            try:
                tt = ttest(target)
                self.progressBarSet(25)
                fold = obiExpression.ExpressionSignificance_FoldChange(self.data, useAttributeLabels=self.genes_in_columns)(target)
                self.progressBarSet(50)
            except ZeroDivisionError, ex:
                tt, fold = [], []
            self.infoLabel.setText("%i genes on input" % len(fold))
            
            invalid = set([key for (key, (t, p)), (_, f) in zip(tt, fold) if any(v is numpy.ma.masked for v in [t, p, f]) or f==0.0])
            tt = [t for t in tt if t[0] not in invalid]
            fold = [f for f in fold if f[0] not in invalid]
            self.progressBarSet(75)
            logratio = numpy.log2(numpy.abs([v for k, v in fold]))
            logpval = -numpy.log10([p for k, (t, p) in tt])
            self.values = dict(zip([k for k, v in tt], zip(logratio, logpval)))
            if not self.values:
                self.warning(1, "Could not compute statistics for any genes!")
            self.progressBarFinished()
        self.graph.setPlotValues(self.values)
        self.setAxesTitles()
        self.updateTooltips()

    def setAxesTitles(self):
        self.graph.setAxisTitle(QwtPlot.xBottom, "log<sub>2</sub> (ratio)" if self.showXTitle else "")
        self.graph.setAxisTitle(QwtPlot.yLeft, "-log<sub>10</sub> (p_value)" if self.showYTitle else "")

    def updateTooltips(self):
        self.graph.tips.removeAll()
        for key, (logratio, logpval) in self.values.items():
            self.graph.tips.addToolTip(logratio, logpval, "<b>%s</b><hr>log<sub>2</sub>(ratio): %.5f<br>p-value: %.5f" \
                                       %(str(key) if self.genes_in_columns else key.name, logratio, math.pow(10, -logpval)))

    def selection(self, items=None):
        """
        Return the current selection.
        """
        if items is None:
            items = sorted(self.values.items())
        values = [val for key, val in items]
        test = self.graph.selection.testSelection(values)
        return test

    def on_selection_changed(self):
        # Selection area on the plot has changed.
        nselected = self.graph.selectedCurve.data().size()
        self.infoLabel2.setText("%i selected genes" % nselected)

        if self.auto_commit:
            selection = list(self.selection())
            # Did the selection actually change
            if selection != self.current_selection:
                self.current_selection = selection
                self.commit()
        else:
            self.selection_changed_flag = True
    
    def commit(self):
        if self.data and self.genes_in_columns:
            items = sorted(self.values.items())
            test = self.selection(items)
            selected = [self.data[i] for t, (i, value) in zip(test, items) if t]
            if selected:
                data = orange.ExampleTable(self.data.domain, selected)
            else:
                data = None
            self.current_selection = list(test) # For testing in on_selection_changed
        elif self.data:
            attrs = [(attr, self.values[attr])  for attr in self.data.domain.attributes if attr in self.values]
            test = self.selection(attrs)
#            test = self.graph.selection.testSelection([val for attr, val in attrs])
            selected = [attr for t, (attr, val) in zip(test, attrs) if t]
            newdomain = orange.Domain(selected + [self.data.domain.classVar])
            newdomain.addmetas(self.data.domain.getmetas())
            data = orange.ExampleTable(newdomain, self.data)
            self.current_selection = list(test) # For testing in on_selection_changed
        else:
            data = None
        
        self.send("Examples with selected attributes", data)
        self.selection_changed_flag = False

    def commit_if(self):
        if self.auto_commit:
            self.commit()
        else:
            self.selection_changed_flag = True
            
    def settingsToWidgetCallbacktargets(self, handler, context):
        self.label_selections = list(getattr(context, "label_selections", self.label_selections))
        self.target_group = getattr(context, "target_group", self.target_group)
        
    def settingsFromWidgetCallbacktargets(self, handler, context):
        context.label_selections = list(self.label_selections)
        context.target_group = self.target_group


if __name__ == "__main__":
    ap = QApplication(sys.argv)
    w = OWVulcanoPlot()
    d = orange.ExampleTable("brown-selected.tab")
    w.setData(d)
    w.show()
    ap.exec_()
    w.saveSettings()
