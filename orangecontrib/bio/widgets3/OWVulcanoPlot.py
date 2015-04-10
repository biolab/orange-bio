"""
Volcano plot
------------

"""
import sys

from operator import add
from collections import defaultdict
from contextlib import contextmanager
from functools import reduce

import numpy
import scipy.stats

from PyQt4 import QtGui
from PyQt4.QtGui import (
    QWidget, QVBoxLayout, QPainterPath, QPainter, QItemSelection
)
from PyQt4.QtCore import (
    Qt, QEvent, QObject, QPointF, QRectF, QSize, QPoint, QRect
)
from PyQt4.QtCore import pyqtSignal as Signal

import pyqtgraph as pg

import Orange.data

from Orange.widgets.utils.datacaching import data_hints
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils.itemmodels import PyListModel

NAME = "Volcano Plot"
DESCRIPTION = "Plots fold change vs. p-value.)"
ICON = "icons/VolcanoPlot.svg"
PRIORITY = 1020

INPUTS = [("Data", Orange.data.Table, "set_data")]
OUTPUTS = [("Data subset", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWVulcanoPlot.OWVulcanoPlot"]


@contextmanager
def blocked_signals(qobj):
    """
    Context manager blocking the `qobj` signals.
    """
    state = qobj.signalsBlocked()
    qobj.blockSignals(True)
    try:
        yield qobj
    finally:
        qobj.blockSignals(state)


class LabelSelectionWidget(QWidget):
    """
    A widget for selection of label values.
    """
    label_activated = Signal(int)
    selection_changed = Signal(object, object)

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self._values_model = PyListModel([], parent=self)
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        def group_box(title):
            box = QtGui.QGroupBox(title)
            box.setFlat(True)
            lay = QVBoxLayout()
            lay.setContentsMargins(0, 0, 0, 0)
            box.setLayout(lay)
            return box

        self.labels_combo = QtGui.QComboBox()
        self.values_view = QtGui.QListView(
            selectionMode=QtGui.QListView.ExtendedSelection
        )
        self.values_view.setModel(self._values_model)
        self.labels_combo.activated[int].connect(self.on_label_activated)
        self.values_view.selectionModel().selectionChanged.connect(
            self.on_values_selection
        )

        l_box = group_box("Label")
        v_box = group_box("Values")

        l_box.layout().addWidget(self.labels_combo)
        v_box.layout().addWidget(self.values_view)

        layout.addWidget(l_box)
        layout.addWidget(v_box)

        self.setLayout(layout)

        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                           QtGui.QSizePolicy.Expanding)

        self._block_selection_emit = False

        self.labels = []
        self.set_labels([])

    def clear(self):
        self.labels_combo.clear()
        self._values_model[:] = []
        self.labels = []

    def set_labels(self, labels):
        """Set the labels to display.
        """
        self.clear()
        if isinstance(labels, dict):
            labels = labels.items()

        self.labels = labels
        for label, values in labels:
            self.labels_combo.addItem(label)

        if labels:
            self.set_current_label(0)

    def set_selection(self, label, values):
        """Set the selection to label and values
        """
        if isinstance(label, str):
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
            self.values_view.selectionModel().select(
                selection,  QtGui.QItemSelectionModel.ClearAndSelect)
        else:
            self.values_view.selectionModel().clear()

    def set_current_label(self, index):
        """Set the current label
        """
        self.labels_combo.setCurrentIndex(index)
        label, values = self.labels[index]
        # Block selection changed
        with blocked_signals(self):
            self._values_model[:] = values

    def on_label_activated(self, index):
        label, values = self.labels[index]
        with blocked_signals(self):
            self._values_model[:] = values

        self.label_activated.emit(index)

    def on_values_selection(self, selected, deselected):
        label, values = self.current_selection()
        self.selection_changed.emit(label, values)

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

    def sizeHint(self):
        return QSize(100, 200)


class GraphSelections(QObject):
    """
    Selection manager using a union of rectangle areas.
    """
    selectionStarted = Signal()
    selectionGeometryChanged = Signal()
    selectionFinished = Signal()

    def __init__(self, parent):
        QObject.__init__(self, parent)
        self.selection = []

    def mapToView(self, pos):
        """
        Map the `pos` in viewbox local into the view (plot) coordinates.
        """
        viewbox = self.parent()
        return viewbox.mapToView(pos)

    def start(self, event):
        pos = self.mapToView(event.pos())
        if event.modifiers() & Qt.ControlModifier:
            self.selection.append((pos, pos))
        else:
            self.selection = [(pos, pos)]
        self.selectionStarted.emit()
        self.selectionGeometryChanged.emit()

    def update(self, event):
        pos = self.mapToView(event.pos())
        self.selection[-1] = self.selection[-1][:-1] + (pos,)
        self.selectionGeometryChanged.emit()

    def end(self, event):
        self.update(event)
        self.selectionFinished.emit()

    def testSelection(self, data):
        """
        Given a

        Parameters
        ----------
        data : (N, 2) array
            Point coordinates

        """
        if len(data) == 0:
            return numpy.zeros(0, dtype=bool)

        data = numpy.asarray(data)
        region = QPainterPath()
        for p1, p2 in self.selection:
            region.addRect(QRectF(p1, p2).normalized())

        def test(point):
            return region.contains(QPointF(point[0], point[1]))

        test = numpy.apply_along_axis(test, 1, data)
        return test

    def mousePressEvent(self, event):
        """
        Handle the intercepted mouse event.

        Return True if the event has been handled, False otherwise
        (let the viewbox handle the event).
        """
        if event.button() == Qt.LeftButton:
            self.start(event)
            event.accept()
            return True
        else:
            return False

    def mouseMoveEvent(self, event):
        """
        Handle the intercepted mouse event.

        Return True if the event has been handled, False otherwise
        (let the viewbox handle the event).
        """
        if event.buttons() & Qt.LeftButton:
            self.update(event)
            event.accept()
            return True
        else:
            return False

    def mouseReleaseEvent(self, event):
        """
        Handle the intercepted mouse event.

        Return True if the event has been handled, False otherwise
        (let the viewbox handle the event).
        """
        if event.button() == Qt.LeftButton:
            self.update(event)
            self.end(event)
            event.accept()
            return True
        else:
            return False

    def eventFilter(self, obj, event):
        if obj is self.parent():
            if event.type() == QEvent.GraphicsSceneMousePress:
                return self.mousePressEvent(event)
            elif event.type() == QEvent.GraphicsSceneMouseMove:
                return self.mouseMoveEvent(event)
            elif event.type() == QEvent.GraphicsSceneMouseRelease:
                return self.mouseReleaseEvent(event)

        return QObject.eventFilter(self, obj, event)


class SymetricSelections(GraphSelections):
    """
    Selection manager using two symmetric areas extending to 'infinity'.
    """
    def __init__(self, parent, x=3, y=3):
        GraphSelections.__init__(self, parent)
        # fmax = numpy.finfo(float).max
        fmax = 2 ** 20
        self.selection = [(QPointF(-fmax, fmax), QPointF(-x, y)),
                          (QPointF(fmax, fmax), QPointF(x, y))]
        self.__constraints = Qt.Vertical | Qt.Horizontal

    def __setpos(self, constraint, pos):
        """
        Set the position of the selection rectangles 'free' corner.
        """
        if constraint & Qt.Horizontal:
            self.selection[0][1].setX(-abs(pos.x()))
            self.selection[1][1].setX(abs(pos.x()))

        if constraint & Qt.Vertical:
            self.selection[0][1].setY(pos.y())
            self.selection[1][1].setY(pos.y())
        self.selectionGeometryChanged.emit()

    def __testaxis(self, event):
        viewbox = self.parent()
        widget = event.widget()
        assert widget is not None
        view = widget.parent()
        assert isinstance(view, QtGui.QGraphicsView)
        pos = view.mapFromGlobal(event.screenPos())
        hitarea = view.mapToScene(QRect(pos - QPoint(2, 2), QSize(4, 4)))
        hitarea = viewbox.mapFromScene(hitarea)
        hitarea = viewbox.mapToView(hitarea).boundingRect()
        center = hitarea.center()

        if center.x() < 0:
            hitarea.moveCenter(QPointF(-center.x(), center.y()))

        cx = self.selection[1][1].x()
        cy = self.selection[1][1].y()
        if hitarea.contains(QPointF(cx, cy)):
            axes = Qt.Horizontal | Qt.Vertical
        elif hitarea.bottom() > cy > hitarea.top() and hitarea.left() > cx:
            axes = Qt.Vertical
        elif hitarea.left() < cx < hitarea.right() and hitarea.bottom() > cy:
            axes = Qt.Horizontal
        else:
            axes = 0
        return axes

    def start(self, event):
        constraints = self.__testaxis(event)
        pos = self.mapToView(event.pos())
        self.__constraints = (constraints or (Qt.Vertical | Qt.Horizontal))
        self.selectionStarted.emit()
        self.__setpos(self.__constraints, pos)

    def update(self, event):
        pos = self.mapToView(event.pos())
        self.__setpos(self.__constraints, pos)

    def end(self, event):
        self.update(event)
        self.__constraints = Qt.Vertical + Qt.Horizontal
        self.selectionFinished.emit()

    def testSelection(self, data):
        if len(data) == 0:
            return numpy.empty(0, dtype=bool)

        data = numpy.asarray(data)
        cutoffX = self.selection[1][1].x()
        cutoffY = self.selection[1][1].y()
        return (numpy.abs(data[:, 0]) >= cutoffX) & (data[:, 1] >= cutoffY)

    def cuttoff(self):
        """
        Return the absolute x, y cutoff values.
        """
        return (self.selection[1][1].x(), self.selection[1][1].y())

    def eventFilter(self, obj, event):
        if event.type() == QEvent.GraphicsSceneHoverMove:
            axes = self.__testaxis(event)
            if axes == Qt.Horizontal:
                self.parent().setCursor(Qt.SizeHorCursor)
            elif axes == Qt.Vertical:
                self.parent().setCursor(Qt.SizeVerCursor)
            else:
                self.parent().setCursor(Qt.ArrowCursor)
        return super().eventFilter(obj, event)


class ScatterPlotItem(pg.ScatterPlotItem):
    def paint(self, painter, option, widget):
        if self.opts["antialias"]:
            painter.setRenderHint(QPainter.Antialiasing, True)
        if self.opts["pxMode"]:
            painter.setRenderHint(QPainter.SmoothPixmapTransform, True)
        super().paint(painter, option, widget)


class VolcanoGraph(pg.PlotWidget):
    selectionChanged = Signal()

    #: Selection Modes
    NoSelection, SymetricSelection, RectSelection = 0, 1, 2

    def __init__(self, parent=None, symbolSize=5, **kwargs):
        pg.PlotWidget.__init__(self, parent, **kwargs)

        self.getAxis("bottom").setLabel("log<sub>2</sub> (ratio)")
        self.getAxis("left").setLabel("-log<sub>10</sub> (P_value)")

        # Absolute cutoff values for symmetric selection mode.
        self.cutoffX = 2.0
        self.cutoffY = 2.0
        # maximum absolute x, y values.
        self.maxX, self.maxY = 10, 10
        self.symbolSize = symbolSize

        self.__selectionMode = VolcanoGraph.NoSelection
        self.__selectiondelegate = None
        self._item = None
        self._selitem = None

        self.plotData = numpy.empty((0, 2))

    def setSelectionMode(self, mode):
        if mode != self.__selectionMode:
            viewbox = self.getViewBox()
            viewbox.setAcceptHoverEvents(True)
            if self.__selectiondelegate is not None:
                viewbox.removeEventFilter(self.__selectiondelegate)
                self.__selectiondelegate.deleteLater()
                self.__selectiondelegate = None

            self.__selectionMode = mode
            if mode == VolcanoGraph.SymetricSelection:
                self.__selectiondelegate = SymetricSelections(
                    viewbox,
                    x=min(self.maxX, self.cutoffX),
                    y=min(self.maxY, self.cutoffY)
                )
                viewbox.installEventFilter(self.__selectiondelegate)
            elif mode == VolcanoGraph.RectSelection:
                self.__selectiondelegate = GraphSelections(viewbox)
                viewbox.installEventFilter(self.__selectiondelegate)
            else:
                pass

            if self.__selectiondelegate is not None:
                self.__selectiondelegate.selectionGeometryChanged.connect(
                    self.__on_selectionChanged)

            if self.plotData.size:
                self.updateSelectionArea()

    def selectionMode(self):
        return self.__selectionMode

    def __on_selectionChanged(self):
        if self.__selectionMode == VolcanoGraph.SymetricSelection:
            self.cutoffX, self.cutoffY = self.__selectiondelegate.cuttoff()

        self.updateSelectionArea()

        self.selectionChanged.emit()

    def selectionMask(self):
        if self.__selectiondelegate is not None:
            return self.__selectiondelegate.testSelection(self.plotData)
        else:
            return numpy.zeros((len(self.plotData), ), dtype=bool)

    def setPlotData(self, data):
        self.plotData = numpy.asarray(data)
        if self.plotData.size:
            self.maxX = numpy.max(numpy.abs(self.plotData[:, 0]))
            self.maxY = numpy.max(self.plotData[:, 1])
        else:
            self.maxX, self.maxY = 10, 10

        self.replot_()
        self.selectionChanged.emit()

    def setSymbolSize(self, size):
        if size != self.symbolSize:
            self.symbolSize = size
            if self._item is not None:
                self._item.setSize(size)

    def updateSelectionArea(self):
        mask = self.selectionMask()
        brush = numpy.array([QtGui.QBrush(Qt.blue), QtGui.QBrush(Qt.red)])
        brush = brush[mask.astype(int)]
        self._item.setBrush(brush)

        if self._selitem is not None:
            self.removeItem(self._selitem)
            self._selitem = None

        if self.__selectiondelegate is not None:
            self._selitem = QtGui.QGraphicsItemGroup()
            self._selitem.dataBounds = \
                lambda axis, frac=1.0, orthoRange=None: None

            for p1, p2 in self.__selectiondelegate.selection:
                r = QRectF(p1, p2).normalized()
                ritem = QtGui.QGraphicsRectItem(r, self._selitem)
                ritem.setBrush(QtGui.QBrush(Qt.NoBrush))
                ritem.setPen(QtGui.QPen(Qt.red, 0))

            self.addItem(self._selitem)

    def replot_(self):
        if self.plotData.size:
            if self._item is not None:
                self.removeItem(self._item)

            mask = self.selectionMask()
            brush = numpy.array([QtGui.QBrush(Qt.blue), QtGui.QBrush(Qt.red)])
            brush = brush[mask.astype(int)]

            self._item = ScatterPlotItem(
                x=self.plotData[:, 0], y=self.plotData[:, 1],
                brush=brush, size=self.symbolSize, antialias=True,
                data=numpy.arange(len(self.plotData))
            )
            self.addItem(self._item)

            self.updateSelectionArea()
        else:
            self.removeItem(self._item)
            self._item = None

    def clear(self):
        self._item = None
        self._selitem = None
        self.plotData = numpy.empty((0, 2))
        super().clear()
        self.selectionChanged.emit()

    def sizeHint(self):
        return QSize(500, 500)

# from .OWGenotypeDistances import SetContextHandler


class OWVolcanoPlot(widget.OWWidget):
    name = NAME
    description = DESCRIPTION
    icon = "../widgets/icons/VulcanoPlot.svg"
    priority = PRIORITY

    inputs = INPUTS
    outputs = OUTPUTS

#     contextHandlers = {
#         "targets": SetContextHandler("targets")
#     }

    want_save_graph = False

    symbol_size = settings.Setting(5)
    symetric_selections = settings.Setting(True)
    auto_commit = settings.Setting(False)

    genes_in_columns = settings.Setting(False)
    label_selections = settings.Setting([])

    def __init__(self, parent=None):
        widget.OWWidget.__init__(self, parent)

        self.target_group = None, []
        self.label_selections = []

        self.data = None
        self.target_group = None, []
        self.current_selection = []

        self.graph = VolcanoGraph(symbolSize=self.symbol_size, background="w")
        self.graph.setSelectionMode(VolcanoGraph.SymetricSelection)

        self.graph.selectionChanged.connect(self.on_selection_changed)
        self.graph.scene().installEventFilter(self)
        self.graph.getViewBox().setMouseEnabled(False, False)
        self.graph.getViewBox().enableAutoRange(enable=True)
        self.mainArea.layout().addWidget(self.graph)

        ## GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.label(box, self, "")
        self.infoLabel.setText("No data on input")
        self.infoLabel2 = gui.label(box, self, "")
        self.infoLabel2.setText("0 selected genes")

        box = gui.widgetBox(self.controlArea, "Target Labels")

        self.target_widget = LabelSelectionWidget(self)
        self.target_widget.selection_changed.connect(self.on_target_changed)
        self.target_widget.label_activated.connect(self.on_label_activated)

        box.layout().addWidget(self.target_widget)

        self.genesInColumnsCheck = gui.checkBox(
            box, self, "genes_in_columns", "Genes in columns",
            callback=[self.init_from_data, self.plot]
        )

        box = gui.widgetBox(self.controlArea, "Settings")

        gui.hSlider(box, self, "symbol_size", label="Symbol size:   ",
                    minValue=2, maxValue=20, step=1,
                    callback=lambda:
                        self.graph.setSymbolSize(self.symbol_size))

        gui.checkBox(box, self,  "symetric_selections", "Symmetric selection",
                     callback=self.__on_selectionModeChanged)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit")
        gui.rubber(self.controlArea)

    def sizeHint(self):
        return QSize(800, 600)

    def clear(self):
        self.target_widget.set_labels([])
        self.targets = []
        self.label_selections = []
        self.target_group = None, []
        self.clear_graph()

    def clear_graph(self):
        self.graph.clear()

    def set_data(self, data=None):
#         self.closeContext("targets")
        self.clear()
        self.data = data
        self.error(0)
        self.warning([0, 1])
        if data is not None:
            self.genes_in_columns = not bool(data.domain.class_var)
            self.genesInColumnsCheck.setEnabled(bool(data.domain.class_var))
            if self.genes_in_columns:
                self.genes_in_columns = not data_hints.get_hint(data, "genesinrows", not self.genes_in_columns)
#             self.openContext("", data)
        else:
            self.infoLabel.setText("No data on input.")
        self.init_from_data()

    def init_from_data(self):
        """Init widget state from the data.
        """
        self.update_target_labels()
        self.error(0)
        if self.data:
            if not self.targets:
                if self.genes_in_columns:
                    self.error(0, "Data set with no column labels (attribute tags)")
                else:
                    self.error(0, "Data has no class.")

#         self.openContext("targets", [(label, v) for label, vals in self.targets \
#                                                 for v in vals])

        # Context handler might 'helpfully' restore some invalid selection
        if len(self.label_selections) != len(self.targets):
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
                items = reduce(add, map(list, items), [])

                targets = defaultdict(set)
                for label, value in items:
                    targets[label].add(value)

                targets = [(key, list(sorted(vals)))
                           for key, vals in targets.items() if len(vals) >= 2]
                self.targets = targets
            else:
                var = self.data.domain.class_var
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
        """Try to restore a saved selection.
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

    def __on_selectionModeChanged(self):
        if self.symetric_selections:
            self.graph.setSelectionMode(VolcanoGraph.SymetricSelection)
        else:
            self.graph.setSelectionMode(VolcanoGraph.RectSelection)

    def group_indices(self):
        """
        Return the selection masks for the current group selection.
        """
        key, values = self.target_group
        if self.genes_in_columns:
            target = set([(key, value) for value in values])
            I1 = [bool(set(var.attributes.items()).intersection(target))
                  for var in self.data.domain.attributes]
            I1 = numpy.array(I1)
            I2 = numpy.ones_like(I1, dtype=bool)
            I2[I1] = False
            return I1, I2
        else:
            target = set(values)
            class_var = self.data.domain.class_var
            target_ind = [class_var.values.index(t) for t in target]

            X, _ = self.data.get_column_view(self.data.domain.class_var)
            I1 = numpy.zeros_like(X, dtype=bool)
            for i in target_ind:
                I1 |= X == i
            I2 = numpy.ones_like(I1, dtype=bool)
            I2[I1] = False
            return I1, I2

    def plot(self):
        self.graph.clear()

        self.current_selection = []
        _, target_values = self.target_group
        self.warning([0, 1])
        self.error(1)

        if self.data and target_values:
            X = self.data.X
            I1, I2 = self.group_indices()
            if not self.genes_in_columns:
                X = X.T

            N1, N2 = numpy.count_nonzero(I1), numpy.count_nonzero(I2)

            if not N1 or not N2:
                self.error(
                    1, "Target labels most exclude/include at least one value."
                )

            if N1 < 2 and N2 < 2:
                self.warning(
                    0, "Insufficient data to compute statistics. "
                       "More than one measurement per class should be provided"
                )

            fold = numpy.log2(numpy.mean(X[:, I1], axis=1) /
                              numpy.mean(X[:, I2], axis=1))

            # TODO: handle missing values better (mstats)
            _, P = scipy.stats.ttest_ind(X[:, I1], X[:, I2], axis=1,
                                         equal_var=True)
            logP = numpy.log10(P)

            mask = numpy.isfinite(fold) & numpy.isfinite(logP)
            self.validindices = numpy.flatnonzero(mask)
            self.graph.setPlotData(numpy.array([fold[mask], -logP[mask]]).T)

            self.infoLabel.setText("%i genes on input" % len(fold))
            # ("{displayed} displayed, {undef} with undefined ratio "
            #  "or t-statistics.")

            if not len(numpy.flatnonzero(mask)):
                self.warning(1, "Could not compute statistics for any genes!")

    def selection(self):
        """
        Return the current item selection mask.
        """
        mask = self.graph.selectionMask()
        return self.validindices[mask]

    def on_selection_changed(self):
        # Selection area on the plot has changed.
        mask = self.graph.selectionMask()
        nselected = numpy.count_nonzero(mask)
        self.infoLabel2.setText("%i selected genes" % nselected)

        if self.auto_commit:
            selection = list(self.selection())
            # Did the selection actually change
            if selection != self.current_selection:
                self.current_selection = selection
                self.commit()

    def commit(self):
        """
        Commit (send) the selected items.
        """
        if self.data and self.genes_in_columns:
            selected = self.selection()
            if selected.size:
                data = self.data[selected]
            else:
                data = None
            self.current_selection = list(selected)  # For testing in on_selection_changed
        elif self.data:
            selected = self.selection()
            attrs = [self.data.domain[i] for i in selected]
            data = self.data[:, tuple(attrs) + self.data.domain.metas]

            self.current_selection = list(selected)  # For testing in on_selection_changed
        else:
            data = None

        self.send("Data subset", data)

    def _handleHelpEvent(self, event):
        # Handle a help event for the plot
        viewbox = self.graph.getViewBox()
        pos = viewbox.mapToView(viewbox.mapFromScene(event.scenePos()))
        logr = pos.x()
        logp = pos.y()
        points = self.graph._item.pointsAt(pos)
        if points:
            point = points[0]
            pos = point.pos()
            logr, logp = pos.x(), pos.y()
            index = point.data()

            text = ("<b>{index}</b><hr/>log<sub>2</sub>(ratio): {logr:.5f}<br/>"
                    "p-value: {p:.5f}").format(logr=logr, p=10 ** -logp,
                                               index=index)

            if len(points) > 1:
                text += "<br/>... (and {} not shown)".format(len(points) - 1)
        else:
            text = ""

        QtGui.QToolTip.showText(event.screenPos(), text, widget=self.graph)
        return True

    def eventFilter(self, obj, event):
        if obj is self.graph.scene():
            if event.type() == QEvent.GraphicsSceneHelp:
                return self._handleHelpEvent(event)
        return super().eventFilter(obj, event)


def test_main():
    ap = QtGui.QApplication(sys.argv)
    w = OWVolcanoPlot()
    w.setAttribute(Qt.WA_DeleteOnClose, True)

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "brown-selected"

    d = Orange.data.Table(filename)
    w.set_data(d)
    w.show()
    w.raise_()
    rval = ap.exec_()
    w.saveSettings()
    return rval


if __name__ == "__main__":
    sys.exit(test_main())
