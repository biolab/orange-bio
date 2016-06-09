import sys
import operator
from collections import defaultdict
from contextlib import contextmanager
from functools import reduce
from xml.sax.saxutils import escape

import numpy
import scipy.stats

from AnyQt.QtWidgets import (
    QListView, QGraphicsView, QGraphicsScene, QSizePolicy,
    QGraphicsGridLayout, QGraphicsLayoutItem, QGraphicsSimpleTextItem,
    QGraphicsWidget, QGraphicsLineItem
)
from AnyQt.QtGui import QPainter, QColor, QPen
from AnyQt.QtCore import Qt, QEvent, QSizeF, QItemSelectionModel

import Orange
from Orange.widgets import widget, gui, settings
from Orange.widgets.utils import itemmodels

from ..utils import group as exp

from .utils.settings import SetContextHandler


def take_columns(data, indices):
    """
    Return an array from the data's columns indexed by indices.

    None values in `indices` are interpreted as indexing a missing
    `column`, i.e. all the values are NaN.

    The resulting array is of shape (N * k,) where the data from
    column 0 is located in `res[:N]`, column2 in `res[N: 2 * N]`, ...
    """
    columns = []
    for index in indices:
        if index is not None:
            col_data, _ = data.get_column_view(index)
        else:
            col_data = numpy.full(len(data), numpy.nan)
        columns.append(col_data)
    return numpy.hstack(columns)


def dist_eucl(a1, a2):
    mask = ~(numpy.isnan(a1) | numpy.isnan(a2))
    a1, a2 = a1[mask], a2[mask]
    return numpy.sqrt(numpy.sum((a1 - a2) ** 2))


def dist_pcorr(a1, a2):
    mask = ~(numpy.isnan(a1) | numpy.isnan(a2))
    a1, a2 = a1[mask], a2[mask]
    P = numpy.corrcoef(a1, a2)[0, 1]
    return (1.0 - P) / 2.0


def dist_spearman(a1, a2):
    mask = ~(numpy.isnan(a1) | numpy.isnan(a2))
    a1, a2 = a1[mask], a2[mask]
    R, _ = scipy.stats.spearmanr(a1, a2)
    return (1.0 - R) / 2.0


@contextmanager
def widget_disable(widget):
    """A context to disable the widget (enabled property)
    """
    widget.setEnabled(False)
    try:
        yield
    finally:
        widget.setEnabled(True)


@contextmanager
def disable_updates(widget):
    """
    A context that sets '_disable_updates' member to True, and
    then restores it.

    """
    widget._disable_updates = True
    try:
        yield
    finally:
        widget._disable_updates = False


def group_label(splits, groups):
    """Return group label.
    """
    labels = ["%s=%s" % (split, group) \
              for split, group in zip(splits, groups)]
    return " | ".join(labels)


def sort_label(sort, attr):
    """Return within group sorted items label for attribute.
    """
    items = [(key, attr.attributes.get(key, "?")) \
             for key in sort]
    labels = ["%s=%s" % tuple(item) for item in items]
    return " | ".join(labels)


def float_if_posible(val):
    """Return val as float if possible otherwise return the value unchanged.
    """
    try:
        return float(val)
    except ValueError:
        return val


def experiment_description(feature):
    """Return experiment description from ``feature.attributes``.
    """
    text = ""
    if feature.attributes:
        items = feature.attributes.items()
        items = [(escape(key), escape(value)) for key, value in items]
        labels = map("%s = %s".__mod__, items)
        text += "<b>%s</b><br/>" % escape(feature.name)
        text += "<br/>".join(labels)
    return text


class OWQualityControl(widget.OWWidget):
    name = "Quality Control"
    description = "Experiment quality control"
    icon = "../widgets/icons/QualityControl.svg"
    priority = 5000

    inputs = [("Experiment Data", Orange.data.Table, "set_data")]
    outputs = []

    DISTANCE_FUNCTIONS = [("Distance from Pearson correlation",
                           dist_pcorr),
                          ("Euclidean distance",
                           dist_eucl),
                          ("Distance from Spearman correlation",
                           dist_spearman)]

    settingsHandler = SetContextHandler()

    split_by_labels = settings.ContextSetting({})
    sort_by_labels = settings.ContextSetting({})

    selected_distance_index = settings.Setting(0)

    def __init__(self, parent=None):
        super().__init__(parent)

        ## Attributes
        self.data = None
        self.distances = None
        self.groups = None
        self.unique_pos = None
        self.base_group_index = 0

        ## GUI
        box = gui.widgetBox(self.controlArea, "Info")
        self.info_box = gui.widgetLabel(box, "\n")

        ## Separate By box
        box = gui.widgetBox(self.controlArea, "Separate By")
        self.split_by_model = itemmodels.PyListModel(parent=self)
        self.split_by_view = QListView()
        self.split_by_view.setSelectionMode(QListView.ExtendedSelection)
        self.split_by_view.setModel(self.split_by_model)
        box.layout().addWidget(self.split_by_view)

        self.split_by_view.selectionModel().selectionChanged.connect(
            self.on_split_key_changed)

        ## Sort By box
        box = gui.widgetBox(self.controlArea, "Sort By")
        self.sort_by_model = itemmodels.PyListModel(parent=self)
        self.sort_by_view = QListView()
        self.sort_by_view.setSelectionMode(QListView.ExtendedSelection)
        self.sort_by_view.setModel(self.sort_by_model)
        box.layout().addWidget(self.sort_by_view)

        self.sort_by_view.selectionModel().selectionChanged.connect(
            self.on_sort_key_changed)

        ## Distance box
        box = gui.widgetBox(self.controlArea, "Distance Measure")
        gui.comboBox(box, self, "selected_distance_index",
                     items=[name for name, _ in self.DISTANCE_FUNCTIONS],
                     callback=self.on_distance_measure_changed)

        self.scene = QGraphicsScene()
        self.scene_view = QGraphicsView(self.scene)
        self.scene_view.setRenderHints(QPainter.Antialiasing)
        self.scene_view.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.mainArea.layout().addWidget(self.scene_view)

        self.scene_view.installEventFilter(self)

        self._disable_updates = False
        self._cached_distances = {}
        self._base_index_hints = {}
        self.main_widget = None

        self.resize(800, 600)

    def clear(self):
        """Clear the widget state."""
        self.data = None
        self.distances = None
        self.groups = None
        self.unique_pos = None

        with disable_updates(self):
            self.split_by_model[:] = []
            self.sort_by_model[:] = []

        self.main_widget = None
        self.scene.clear()
        self.info_box.setText("\n")
        self._cached_distances = {}

    def set_data(self, data=None):
        """Set input experiment data."""
        self.closeContext()
        self.clear()

        self.error(0)
        self.warning(0)

        if data is not None:
            keys = self.get_suitable_keys(data)
            if not keys:
                self.error(0, "Data has no suitable feature labels.")
                data = None

        self.data = data
        if data is not None:
            self.on_new_data()

    def update_label_candidates(self):
        """Update the label candidates selection GUI 
        (Group/Sort By views).

        """
        keys = self.get_suitable_keys(self.data)
        with disable_updates(self):
            self.split_by_model[:] = keys
            self.sort_by_model[:] = keys

    def get_suitable_keys(self, data):
        """ Return suitable attr label keys from the data where
        the key has at least two unique values in the data.

        """
        attrs = [attr.attributes.items() for attr in data.domain.attributes]
        attrs = reduce(operator.iadd, attrs, [])
        # in case someone put non string values in attributes dict
        attrs = [(str(key), str(value)) for key, value in attrs]
        attrs = set(attrs)
        values = defaultdict(set)
        for key, value in attrs:
            values[key].add(value)
        keys = [key for key in values if len(values[key]) > 1]
        return keys

    def selected_split_by_labels(self):
        """Return the current selected split labels.
        """
        sel_m = self.split_by_view.selectionModel()
        indices = [r.row() for r in sel_m.selectedRows()]
        return [self.sort_by_model[i] for i in indices]

    def selected_sort_by_labels(self):
        """Return the current selected sort labels
        """
        sel_m = self.sort_by_view.selectionModel()
        indices = [r.row() for r in sel_m.selectedRows()]
        return [self.sort_by_model[i] for i in indices]

    def selected_distance(self):
        """Return the selected distance function.
        """
        return self.DISTANCE_FUNCTIONS[self.selected_distance_index][1]

    def selected_base_group_index(self):
        """Return the selected base group index
        """
        return self.base_group_index

    def selected_base_indices(self, base_group_index=None):
        indices = []
        for g, ind in self.groups:
            if base_group_index is None:
                label = group_label(self.selected_split_by_labels(), g)
                ind = [i for i in ind if i is not None]
                i = self._base_index_hints.get(label, ind[0] if ind else None)
            else:
                i = ind[base_group_index]
            indices.append(i)
        return indices

    def on_new_data(self):
        """We have new data and need to recompute all.
        """
        self.closeContext()

        self.update_label_candidates()
        self.info_box.setText(
            "%s genes \n%s experiments" %
            (len(self.data),  len(self.data.domain.attributes))
        )

        self.base_group_index = 0

        keys = self.get_suitable_keys(self.data)
        self.openContext(keys)

        ## Restore saved context settings (split/sort selection)
        split_by_labels = self.split_by_labels
        sort_by_labels = self.sort_by_labels

        def select(model, selection_model, selected_items):
            """Select items in a Qt item model view
            """
            all_items = list(model)
            try:
                indices = [all_items.index(item) for item in selected_items]
            except:
                indices = []
            for ind in indices:
                selection_model.select(model.index(ind),
                                       QItemSelectionModel.Select)

        with disable_updates(self):
            select(self.split_by_view.model(),
                   self.split_by_view.selectionModel(),
                   split_by_labels)

            select(self.sort_by_view.model(),
                   self.sort_by_view.selectionModel(),
                   sort_by_labels)

        with widget_disable(self):
            self.split_and_update()

    def on_split_key_changed(self, *args):
        """Split key has changed
        """
        with widget_disable(self):
            if not self._disable_updates:
                self.base_group_index = 0
                self.split_by_labels = self.selected_split_by_labels()
                self.split_and_update()

    def on_sort_key_changed(self, *args):
        """Sort key has changed
        """
        with widget_disable(self):
            if not self._disable_updates:
                self.base_group_index = 0
                self.sort_by_labels = self.selected_sort_by_labels()
                self.split_and_update()

    def on_distance_measure_changed(self):
        """Distance measure has changed
        """
        if self.data is not None:
            with widget_disable(self):
                self.update_distances()
                self.replot_experiments()

    def on_view_resize(self, size):
        """The view with the quality plot has changed
        """
        if self.main_widget:
            current = self.main_widget.size()
            self.main_widget.resize(size.width() - 6,
                                    current.height())

            self.scene.setSceneRect(self.scene.itemsBoundingRect())

    def on_rug_item_clicked(self, item):
        """An ``item`` in the quality plot has been clicked.
        """
        update = False
        sort_by_labels = self.selected_sort_by_labels()
        if sort_by_labels and item.in_group:
            ## The item is part of the group
            if item.group_index != self.base_group_index:
                self.base_group_index = item.group_index
                update = True

        else:
            if sort_by_labels:
                # If the user clicked on an background item it
                # invalidates the sorted labels selection
                with disable_updates(self):
                    self.sort_by_view.selectionModel().clear()
                    update = True

            index = item.index
            group = item.group
            label = group_label(self.selected_split_by_labels(), group)

            if self._base_index_hints.get(label, 0) != index:
                self._base_index_hints[label] = index
                update = True

        if update:
            with widget_disable(self):
                self.split_and_update()

    def eventFilter(self, obj, event):
        if obj is self.scene_view and event.type() == QEvent.Resize:
            self.on_view_resize(event.size())
        return super().eventFilter(obj, event)

    def split_and_update(self):
        """
        Split the data based on the selected sort/split labels
        and update the quality plot.

        """
        split_labels = self.selected_split_by_labels()
        sort_labels = self.selected_sort_by_labels()

        self.warning(0)
        if not split_labels:
            self.warning(0, "No separate by label selected.")

        self.groups, self.unique_pos = \
                exp.separate_by(self.data, split_labels,
                                consider=sort_labels,
                                add_empty=True)

        self.groups = sorted(self.groups.items(),
                             key=lambda t: list(map(float_if_posible, t[0])))
        self.unique_pos = sorted(self.unique_pos.items(),
                                 key=lambda t: list(map(float_if_posible, t[0])))

        if self.groups:
            if sort_labels:
                group_base = self.selected_base_group_index()
                base_indices = self.selected_base_indices(group_base)
            else:
                base_indices = self.selected_base_indices()
            self.update_distances(base_indices)
            self.replot_experiments()

    def get_cached_distances(self, measure):
        if measure not in self._cached_distances:
            attrs = self.data.domain.attributes
            mat = numpy.zeros((len(attrs), len(attrs)))

            self._cached_distances[measure] = \
                (mat, set(zip(range(len(attrs)), range(len(attrs)))))

        return self._cached_distances[measure]

    def get_cached_distance(self, measure, i, j):
        matrix, computed = self.get_cached_distances(measure)
        key = (i, j) if i < j else (j, i)
        if key in computed:
            return matrix[i, j]
        else:
            return None

    def get_distance(self, measure, i, j):
        d = self.get_cached_distance(measure, i, j)
        if d is None:
            vec_i = take_columns(self.data, [i])
            vec_j = take_columns(self.data, [j])
            d = measure(vec_i, vec_j)

            mat, computed = self.get_cached_distances(measure)
            mat[i, j] = d
            key = key = (i, j) if i < j else (j, i)
            computed.add(key)
        return d

    def store_distance(self, measure, i, j, dist):
        matrix, computed = self.get_cached_distances(measure)
        key = (i, j) if i < j else (j, i)
        matrix[j, i] = matrix[i, j] = dist
        computed.add(key)

    def update_distances(self, base_indices=()):
        """Recompute the experiment distances.
        """
        distance = self.selected_distance()
        if base_indices == ():
            base_group_index = self.selected_base_group_index()
            base_indices = [ind[base_group_index] \
                            for _, ind in self.groups]

        assert(len(base_indices) == len(self.groups))

        base_distances = []
        attributes = self.data.domain.attributes
        pb = gui.ProgressBar(self, len(self.groups) * len(attributes))

        for (group, indices), base_index in zip(self.groups, base_indices):
            # Base column of the group
            if base_index is not None:
                base_vec = take_columns(self.data, [base_index])
                distances = []
                # Compute the distances between base column
                # and all the rest data columns.
                for i in range(len(attributes)):
                    if i == base_index:
                        distances.append(0.0)
                    elif self.get_cached_distance(distance, i, base_index) is not None:
                        distances.append(self.get_cached_distance(distance, i, base_index))
                    else:
                        vec_i = take_columns(self.data, [i])
                        dist = distance(base_vec, vec_i)
                        self.store_distance(distance, i, base_index, dist)
                        distances.append(dist)
                    pb.advance()

                base_distances.append(distances)
            else:
                base_distances.append(None)

        pb.finish()
        self.distances = base_distances

    def replot_experiments(self):
        """Replot the whole quality plot.
        """
        self.scene.clear()
        labels = []

        max_dist = numpy.nanmax(list(filter(None, self.distances)))
        rug_widgets = []

        group_pen = QPen(Qt.black)
        group_pen.setWidth(2)
        group_pen.setCapStyle(Qt.RoundCap)
        background_pen = QPen(QColor(0, 0, 250, 150))
        background_pen.setWidth(1)
        background_pen.setCapStyle(Qt.RoundCap)

        main_widget = QGraphicsWidget()
        layout = QGraphicsGridLayout()
        attributes = self.data.domain.attributes
        if self.data is not None:
            for (group, indices), dist_vec in zip(self.groups, self.distances):
                indices_set = set(indices)
                rug_items = []
                if dist_vec is not None:
                    for i, attr in enumerate(attributes):
                        # Is this a within group distance or background
                        in_group = i in indices_set
                        if in_group:
                            rug_item = ClickableRugItem(dist_vec[i] / max_dist,
                                           1.0, self.on_rug_item_clicked)
                            rug_item.setPen(group_pen)
                            tooltip = experiment_description(attr)
                            rug_item.setToolTip(tooltip)
                            rug_item.group_index = indices.index(i)
                            rug_item.setZValue(rug_item.zValue() + 1)
                        else:
                            rug_item = ClickableRugItem(dist_vec[i] / max_dist,
                                           0.85, self.on_rug_item_clicked)
                            rug_item.setPen(background_pen)
                            tooltip = experiment_description(attr)
                            rug_item.setToolTip(tooltip)

                        rug_item.group = group
                        rug_item.index = i
                        rug_item.in_group = in_group

                        rug_items.append(rug_item)

                rug_widget = RugGraphicsWidget(parent=main_widget)
                rug_widget.set_rug(rug_items)

                rug_widgets.append(rug_widget)

                label = group_label(self.selected_split_by_labels(), group)
                label_item = QGraphicsSimpleTextItem(label, main_widget)
                label_item = GraphicsSimpleTextLayoutItem(label_item, parent=layout)
                label_item.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
                labels.append(label_item)

        for i, (label, rug_w) in enumerate(zip(labels, rug_widgets)):
            layout.addItem(label, i, 0, Qt.AlignVCenter)
            layout.addItem(rug_w, i, 1)
            layout.setRowMaximumHeight(i, 30)

        main_widget.setLayout(layout)
        self.scene.addItem(main_widget)
        self.main_widget = main_widget
        self.rug_widgets = rug_widgets
        self.labels = labels
        self.on_view_resize(self.scene_view.size())


class GraphicsSimpleTextLayoutItem(QGraphicsLayoutItem):
    """ A Graphics layout item wrapping a QGraphicsSimpleTextItem alowing it
    to be managed by a layout.

    """
    def __init__(self, text_item, orientation=Qt.Horizontal, parent=None):
        super().__init__(parent)
        self.orientation = orientation
        self.text_item = text_item
        if orientation == Qt.Vertical:
            self.text_item.rotate(-90)
            self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        else:
            self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)

    def setGeometry(self, rect):
        super().setGeometry(rect)
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

    def updateGeometry(self):
        super().updateGeometry()
        parent = self.parentLayoutItem()
        if parent.isLayout():
            parent.updateGeometry()

    def setFont(self, font):
        self.text_item.setFont(font)
        self.updateGeometry()

    def setText(self, text):
        self.text_item.setText(text)
        self.updateGeometry()


class RugGraphicsWidget(QGraphicsWidget):
    def __init__(self, parent=None, rug=None):
        QGraphicsWidget.__init__(self, parent)
        self.rug_items = []
        self.set_rug(rug)
        self.setMaximumHeight(30)
        self.setMinimumHeight(30)

    def clear(self):
        """
        Clear all rug items from this widget and remove them 
        from the scene.

        """
        for item in self.rug_items:
            item.setParent(None)
            if self.scene() is not None:
                self.scene().removeItem(item)

    def set_rug(self, rug):
        """
        Set the rug items.

        ``rug`` must be a list of floats or already initialized
        instances of RugItem. The widget takes ownership of all
        items.

        """
        rug = rug if rug is not None else []
        self.clear()
        self.add_rug(rug)

    def add_rug(self, rug):
        """
        Add rug items.

        See :obj:`set_rug`

        """
        items = []
        for item in rug:
            if isinstance(item, float):
                item = RugItem(value=item)
                items.append(item)
            elif isinstance(item, RugItem):
                items.append(item)

        for item in items:
            item.setParentItem(self)

        self.rug_items.extend(items)

        self.update_rug_geometry()

    def update_rug_geometry(self):
        """Recompute the rug items positions within this widget.
        """
        size = self.size()
        height = size.height()
        width = size.width()

        for item in self.rug_items:
            offset = (1.0 - item.height) * height / 2.0
            item.setPos(width * item.value, 0)
            item.setLine(0., offset, 0., height - offset)

    def resizeEvent(self, event):
        """Reimplemented from QGraphicsWidget
        """
        QGraphicsWidget.resizeEvent(self, event)
        self.update_rug_geometry()


class RugItem(QGraphicsLineItem):
    def __init__(self, value, height):
        QGraphicsLineItem.__init__(self)
        self.value = value
        self.height = height

    def set_height(self, height):
        """Set the height of this item (in ratio of the rug height)
        """
        self.height = height


class ClickableRugItem(RugItem):
    def __init__(self, value, height, on_pressed):
        RugItem.__init__(self, value, height)
        self.on_pressed = on_pressed
        self.setAcceptedMouseButtons(Qt.LeftButton)
        self.setAcceptHoverEvents(True)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton and self.on_pressed:
            self.on_pressed(self)

    def hoverEnterEvent(self, event):
        pen = QPen(self.pen())
        pen.setWidthF(3)
        self.setPen(pen)
        return RugItem.hoverEnterEvent(self, event)

    def hoverLeaveEvent(self, event):
        pen = QPen(self.pen())
        pen.setWidth(2)
        self.setPen(pen)
        return RugItem.hoverLeaveEvent(self, event)


def test_main(argv=sys.argv):
    from AnyQt.QtWidgets import QApplication
    if len(argv) > 1:
        filename = argv[1]
    else:
        filename = "GDS624.tab"
    data = Orange.data.Table(filename)

    app = QApplication(argv)
    w = OWQualityControl()
    w.set_data(data)
    w.handleNewSignals()
    w.show()
    w.raise_()
    r = app.exec_()
    w.set_data(None)
    w.handleNewSignals()
    w.saveSettings()
    return r

if __name__ == "__main__":
    sys.exit(test_main())
