"""
Volcano plot
------------

"""
import sys

from contextlib import contextmanager
from functools import partial
from typing import List

import numpy
import scipy.stats

from AnyQt.QtWidgets import (
    QGraphicsView, QGraphicsItemGroup, QGraphicsRectItem, QToolTip,
    QFormLayout, QCheckBox
)
from AnyQt.QtGui import QPainter, QBrush, QPen, QStandardItemModel

from AnyQt.QtCore import (
    Qt, QEvent, QObject, QPointF, QRectF, QSize, QPoint, QRect, Signal
)

import pyqtgraph as pg

import Orange.data

from Orange.widgets.utils.datacaching import data_hints
from Orange.widgets import widget, gui, settings

from orangecontrib.bio.widgets3.utils import gui as guiutils
from orangecontrib.bio.widgets3.utils import group as grouputils
from orangecontrib.bio.widgets3.utils.settings import SetContextHandler

NAME = "Volcano Plot"
DESCRIPTION = "Plots fold change vs. p-value."
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

        def contained(a, left, top, right, bottom):
            assert left <= right and bottom <= top
            x, y = a.T
            return (x >= left) & (x <= right) & (y <= top) & (y >= bottom)

        data = numpy.asarray(data)
        selected = numpy.zeros(len(data), dtype=bool)
        for p1, p2 in self.selection:
            r = QRectF(p1, p2).normalized()
            # Note the inverted top/bottom (Qt coordinate system)
            selected |= contained(data, r.left(), r.bottom(),
                                  r.right(), r.top())
        return selected

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
        assert isinstance(view, QGraphicsView)
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

    def cutoff(self):
        """
        Return the absolute x, y cutoff values.
        """
        return (self.selection[1][1].x(), self.selection[1][1].y())

    def setCutoff(self, x, y):
        if (x, y) != self.cutoff():
            self.__setpos(Qt.Vertical | Qt.Horizontal, QPointF(x, y))

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
        self.maxX, self.maxY = 10., 10.
        self.symbolSize = symbolSize

        self.__selectionMode = VolcanoGraph.NoSelection
        self.__selectiondelegate = None
        self._item = None
        self._selitem = None
        self.plotData = numpy.empty((0, 2))
        self._stylebrush = numpy.array([
            pg.mkBrush((0, 0, 255, 100)),  # normal points
            pg.mkBrush((255, 0, 0, 100))   # selected points
        ])

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

            self.selectionChanged.emit()

    def selectionMode(self):
        return self.__selectionMode

    def setSelectionCut(self, x, y):
        if (x, y) != (self.cutoffX, self.cutoffY):
            self.cutoffX, self.cutoffY = x, y
            if self.__selectionMode == VolcanoGraph.SymetricSelection:
                assert isinstance(self.__selectiondelegate, SymetricSelections)
                self.__selectiondelegate.setCutoff(x, y)

    def selectionCut(self):
        return self.cutoffX, self.cutoffY

    def __on_selectionChanged(self):
        if self.__selectionMode == VolcanoGraph.SymetricSelection:
            self.cutoffX, self.cutoffY = self.__selectiondelegate.cutoff()

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
            self.maxX, self.maxY = 10., 10.

        self.replot_()

        self.setRange(QRectF(-self.maxX, 0, 2 * self.maxX, self.maxY))
        self.selectionChanged.emit()

    def setSymbolSize(self, size):
        if size != self.symbolSize:
            self.symbolSize = size
            if self._item is not None:
                self._item.setSize(size)

    def updateSelectionArea(self):
        mask = self.selectionMask()
        brush = self._stylebrush[mask.astype(int)]
        if self._item is not None:
            self._item.setBrush(brush)

        if self._selitem is not None:
            self.removeItem(self._selitem)
            self._selitem = None

        if self.__selectiondelegate is not None:
            self._selitem = QGraphicsItemGroup()
            self._selitem.dataBounds = \
                lambda axis, frac=1.0, orthoRange=None: None

            for p1, p2 in self.__selectiondelegate.selection:
                r = QRectF(p1, p2).normalized()
                ritem = QGraphicsRectItem(r, self._selitem)
                ritem.setBrush(QBrush(Qt.NoBrush))
                ritem.setPen(QPen(Qt.red, 0))

            self.addItem(self._selitem)

    def replot_(self):
        if self.plotData.size:
            if self._item is not None:
                self.removeItem(self._item)

            mask = self.selectionMask()
            brush = self._stylebrush[mask.astype(int)]
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

    def clearPlot(self):
        self._item = None
        self._selitem = None
        self.plotData = numpy.empty((0, 2))
        self.clear()
        self.selectionChanged.emit()

    def sizeHint(self):
        return QSize(500, 500)


class OWVolcanoPlot(widget.OWWidget):
    name = NAME
    description = DESCRIPTION
    icon = "../widgets/icons/VulcanoPlot.svg"
    priority = PRIORITY

    inputs = INPUTS
    outputs = OUTPUTS

    settingsHandler = SetContextHandler()

    symbol_size = settings.Setting(5)
    symetric_selections = settings.Setting(True)
    auto_commit = settings.Setting(False)

    current_group_index = settings.ContextSetting(-1)
    #: stored selection indices (List[int]) for every split group.
    stored_selections = settings.ContextSetting([])  # type: List[int]

    #: The (abs) selection cut on the log(ratio) scale (X axis)
    cut_log_r = settings.Setting(2.0)  # type: float
    #: The selection cut on on log(p_value) scale (Y axis)
    cut_log_p = settings.Setting(2.0)  # type: float

    graph_name = "graph.plotItem"

    def __init__(self, parent=None):
        widget.OWWidget.__init__(self, parent)

        self.data = None
        self.target_group = None, []
        self.targets = []
        self.current_selection = []
        self.validindices = numpy.empty((0,), dtype=int)

        self.graph = VolcanoGraph(symbolSize=self.symbol_size, background="w")
        self.graph.setSelectionCut(self.cut_log_r, self.cut_log_p)
        self.graph.setSelectionMode(
            VolcanoGraph.SymetricSelection if self.symetric_selections else
            VolcanoGraph.RectSelection)

        self.graph.selectionChanged.connect(self.on_selection_changed)
        self.graph.scene().installEventFilter(self)
        self.graph.getViewBox().setMouseEnabled(False, False)
        self.graph.getViewBox().enableAutoRange(enable=True)
        self.graph.getViewBox().setMenuEnabled(False)
        self.graph.getPlotItem().hideButtons()
        self.graph.plotItem.showGrid(True, True, 0.3)
        self.mainArea.layout().addWidget(self.graph)

        ## GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.infoLabel = gui.label(box, self, "")
        self.infoLabel.setText("No data on input")
        self.infoLabel2 = gui.label(box, self, "")
        self.infoLabel2.setText("0 selected genes")

        box = gui.widgetBox(self.controlArea, "Target Labels")

        self.target_widget = guiutils.LabelSelectionWidget(
            self,
            groupChanged=self.on_label_activated,
            groupSelectionChanged=self.on_target_changed)

        box.layout().addWidget(self.target_widget)

        box = gui.widgetBox(self.controlArea, "Selection")
        cb = gui.checkBox(box, self, "symetric_selections", "Symmetric selection",
                          callback=self.__on_selectionModeChanged)  # type: QCheckBox
        form = QFormLayout(
            labelAlignment=Qt.AlignRight,
            formAlignment=Qt.AlignLeft,
            fieldGrowthPolicy=QFormLayout.AllNonFixedFieldsGrow,
        )

        form.addRow(
            "|log(ratio)|:",
            gui.spin(box, self, "cut_log_r", 0.0, 20.0,
                     step=0.01, decimals=4, spinType=float,
                     callback=self.on_cut_changed)
        )
        form.addRow(
            "-log<sub>10</sub>(P_value):",
            gui.spin(box, self, "cut_log_p", 0.0, 20.0,
                     step=0.01, decimals=4, spinType=float,
                     callback=self.on_cut_changed)
        )
        box.layout().addLayout(form)

        def set_enable_all(layout, enabled):
            for i in range(layout.count()):
                item = layout.itemAt(i)
                if item is not None and item.widget() is not None:
                    item.widget().setEnabled(enabled)
        set_enable_all(form, self.symetric_selections)
        cb.toggled.connect(partial(set_enable_all, form))
        box = gui.widgetBox(self.controlArea, "Settings")

        gui.hSlider(box, self, "symbol_size", label="Symbol size:   ",
                    minValue=2, maxValue=20, step=1,
                    callback=lambda:
                        self.graph.setSymbolSize(self.symbol_size))

        gui.rubber(self.controlArea)
        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit")

    def sizeHint(self):
        return QSize(800, 600)

    def clear(self):
        self.targets = []
        self.stored_selections = []
        self.target_widget.setModel(None)
        self.validindices = numpy.empty((0,), dtype=int)
        self.clear_graph()

    def clear_graph(self):
        self.graph.clearPlot()

    def set_data(self, data=None):
        self.closeContext()
        self.clear()
        self.data = data
        self.error(0)
        self.warning([0, 1])
        if data is not None:
            self.init_from_data()

            if not self.targets:
                self.error(0, "Data has no suitable defined groups/labels")

            genesinrows = data_hints.get_hint(data, "genesinrows", True)

            # Select the first RowGroup if genesinrows hint is False
            if not genesinrows:
                ind = [i for i, grp in enumerate(self.targets)
                       if isinstance(grp, grouputils.RowGroup)]
                if ind:
                    self.current_group_index = ind[0]
                else:
                    self.current_group_index = 0 if self.targets else -1

            else:
                self.current_group_index = 0

            # TODO: Why does the current_group_index not get restored
            self.openContext([(grp.name, v)
                              for grp in self.targets for v in grp.values])
            model = self.target_widget.model()

            # FIX: This assumes fixed target key/values order.
            indices = [model.index(i, 0, model.index(keyind, 0,))
                       for keyind, selected in enumerate(self.stored_selections)
                       for i in selected]

            selection = guiutils.itemselection(indices)
            self.target_widget.setCurrentGroupIndex(self.current_group_index)
            self.target_widget.setSelection(selection)

            self.plot()
        else:
            self.infoLabel.setText("No data on input.")

        self.unconditional_commit()

    def init_from_data(self):
        """Initialize widget state after receiving new data.
        """
        if self.data is not None:
            column_groups, row_groups = grouputils.group_candidates(self.data)
            self.targets = column_groups + row_groups
            self.stored_selections = [[0] for _ in self.targets]

            targetitems = [guiutils.standarditem_from(desc)
                           for desc in self.targets]
            model = QStandardItemModel()
            for item in targetitems:
                model.appendRow(item)

            with blocked_signals(self.target_widget):
                self.target_widget.setModel(model)

        else:
            self.targets = []
            self.stored_selections = []
            with blocked_signals(self.target_widget):
                self.target_widget.setModel(None)

    def selected_split(self):
        groupindex = self.target_widget.currentGroupIndex()
        if not (0 <= groupindex < len(self.targets)):
            return None, []
        group = self.targets[groupindex]
        selection = self.target_widget.currentGroupSelection()
        selection = [ind.row() for ind in selection.indexes()]
        return group, selection

    def on_label_activated(self, index):
        self.current_group_index = index
        self.plot()

    def on_target_changed(self):
        # Save the current selection to persistent settings
        rootind = self.target_widget.currentGroupIndex()
        selection = self.target_widget.currentGroupSelection()
        indices = [ind.row() for ind in selection.indexes()]

        if rootind != -1 and indices:
            self.stored_selections[rootind] = indices

        self._update_plot()

    def on_cut_changed(self):
        self.graph.setSelectionCut(self.cut_log_r, self.cut_log_p)

    def _update_plot(self):
        group, targetindices = self.selected_split()
        if group is not None and targetindices:
            self.plot()
        else:
            self.clear_graph()

    def __on_selectionModeChanged(self):
        if self.symetric_selections:
            self.graph.setSelectionMode(VolcanoGraph.SymetricSelection)
        else:
            self.graph.setSelectionMode(VolcanoGraph.RectSelection)

    def plot(self):
        self.graph.clearPlot()
        self.validindices = numpy.empty((0,), dtype=int)
        self.current_selection = []
        group, target_indices = self.selected_split()
        self.warning([0, 1])
        self.error(1)

        if self.data and group is not None and target_indices:
            X = self.data.X
            I1 = grouputils.group_selection_mask(
                self.data, group, target_indices)
            I2 = ~I1
            if isinstance(group, grouputils.RowGroup):
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

            X1, X2 = X[:, I1], X[:, I2]
            if numpy.any(X1 < 0.0) or numpy.any(X2 < 0):
                self.error(
                    "Negative values in the input. The inputs cannot be in "
                    "ratio scale."
                )
                X1 = numpy.full_like(X1, numpy.nan)
                X2 = numpy.full_like(X2, numpy.nan)

            with numpy.errstate(divide="ignore", invalid="ignore"):
                fold = numpy.log2(numpy.mean(X1, axis=1) /
                                  numpy.mean(X2, axis=1))
                # TODO: handle missing values better (mstats)
                _, P = scipy.stats.ttest_ind(X1, X2, axis=1, equal_var=True)
                logP = numpy.log10(P)
                if numpy.isscalar(logP):
                    # ttest_ind does not preserve output shape if either
                    # a or b is empty
                    logP = numpy.full(fold.shape, numpy.nan)

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
        if self.graph.selectionMode() == VolcanoGraph.SymetricSelection:
            self.cut_log_r, self.cut_log_p = self.graph.selectionCut()

        if self.data is not None:
            mask = self.graph.selectionMask()
            nselected = numpy.count_nonzero(mask)
            self.infoLabel2.setText("%i selected genes" % nselected)

            selection = list(self.selection())
            # Did the selection actually change
            if selection != self.current_selection:
                self.current_selection = selection
                self.commit()

    def commit(self):
        """
        Commit (send) the selected items.
        """
        if self.data is not None and \
                isinstance(self.selected_split()[0], grouputils.ColumnGroup):
            selected = self.selection()
            if selected.size:
                data = self.data[selected]
            else:
                data = None
            self.current_selection = list(selected)  # For testing in on_selection_changed
        elif self.data is not None:
            selected = self.selection()
            attrs = [self.data.domain[i] for i in selected]
            domain = Orange.data.Domain(
                attrs, self.data.domain.class_vars, self.data.domain.metas)
            data = self.data.from_table(domain, self.data)
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
        text = ""
        points = []
        if self.graph._item is not None:
            points = self.graph._item.pointsAt(pos)
        if points:
            point = points[0]
            pos = point.pos()
            logr, logp = pos.x(), pos.y()
            index = point.data()

            text = ("<b>{index}</b><hr/>"
                    "log<sub>2</sub>(ratio): {logr:.5f}<br/>"
                    "p-value: {p:.5f}").format(logr=logr, p=10 ** -logp,
                                               index=index)

            if len(points) > 1:
                text += "<br/>... (and {} not shown)".format(len(points) - 1)
        else:
            text = ""

        QToolTip.showText(event.screenPos(), text, widget=self.graph)
        return True

    def eventFilter(self, obj, event):
        if obj is self.graph.scene():
            if event.type() == QEvent.GraphicsSceneHelp:
                return self._handleHelpEvent(event)
        return super().eventFilter(obj, event)

    def send_report(self):
        group, selection = self.selected_split()
        self.report_plot()
        caption = []
        if group:
            target = group.name
            values = ", ".join(numpy.array(group.values)[selection])
            caption.append("{var} = {value}".format(var=target, value=values))
        caption.append(self.infoLabel2.text())
        self.report_caption(", ".join(caption))


def main(argv=None):
    from AnyQt.QtWidgets import QApplication
    app = QApplication(list(argv) if argv else [])
    argv = app.arguments()
    w = OWVolcanoPlot()

    if len(argv) > 1:
        filename = argv[1]
    else:
        filename = "geo-gds360"
    d = Orange.data.Table(filename)
    w.set_data(d)
    w.handleNewSignals()
    w.show()
    w.raise_()
    rval = app.exec_()
    w.set_data(None)
    w.handleNewSignals()
    w.saveSettings()
    w.onDeleteWidget()
    return rval


if __name__ == "__main__":
    sys.exit(main(sys.argv))
