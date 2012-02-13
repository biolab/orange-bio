"""
<name>Quality Control</name>

"""

import sys

import Orange

from OWWidget import *
from OWItemModels import PyListModel
from OWGraphics import GraphicsSimpleTextLayoutItem, GtI
import OWGUI

from OWGenotypeDistances import SetContextHandler

import obiExperiments as exp
import numpy

from collections import defaultdict
from contextlib import contextmanager
from pprint import pprint

DEBUG = False

@contextmanager
def control_disable(widget):
    widget.setEnabled(False)
    yield
    widget.setEnabled(True)
    
@contextmanager
def disable_updates(widget):
    widget._disable_updates = True
    yield
    widget._disable_updates = False

def group_label(splits, groups):
    """
    Return group label
    """
    labels = ["{}={}".format(split, group) \
              for split, group in zip(splits, groups)]
    return " | ".join(labels)


def sort_label(sort, attr):
    """
    Return within group sorted items label for attribute.
    """
    items = [(key, attr.attributes.get(key, "?")) \
             for key in sort]
    labels = ["{}={}".format(*item) for item in items]
    return " | ".join(labels)

def float_if_posible(val):
    try:
        return float(val)
    except ValueError:
        return val

class OWQualityControl(OWWidget):
    contextHandlers = {"": SetContextHandler("")}
    settingsList = []
    
    DISTANCE_FUNCTIONS = [("Distance from Pearson correlation", exp.dist_pcorr),
                          ("Euclidean distance", exp.dist_eucl),
                          ("Distance from Spearman correlation", exp.dist_spearman)]

    def __init__(self, parent=None, signalManager=None,
                 title="Quality Control"):
        OWWidget.__init__(self, parent, signalManager, title,
                          wantGraph=True)

        self.inputs = [("Experiment Data", Orange.data.Table, self.set_data)]

        ## Settings
        self.selected_distance_index = 0
        self.replicate_id = "replicate"

        ## Attributes
        self.data = None
        self.distances = None
        self.groups = None
        self.unique_pos = None
        self.base_index = 0

        ## GUI
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info_box = OWGUI.widgetLabel(box, "\n")

        ## Split By box
        box = OWGUI.widgetBox(self.controlArea, "Split By")
        self.split_by_model = PyListModel()
        self.split_by_view = QListView()
        self.split_by_view.setSelectionMode(QListView.ExtendedSelection)
        self.split_by_view.setModel(self.split_by_model)
        box.layout().addWidget(self.split_by_view)

        self.connect(self.split_by_view.selectionModel(),
                     SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),
                     self.on_split_key_changed)

        ## Sort By box
        box = OWGUI.widgetBox(self.controlArea, "Sort By")
        self.sort_by_model = PyListModel()
        self.sort_by_view = QListView()
        self.sort_by_view.setSelectionMode(QListView.ExtendedSelection)
        self.sort_by_view.setModel(self.sort_by_model)
        box.layout().addWidget(self.sort_by_view)
        
        self.connect(self.sort_by_view.selectionModel(),
                     SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),
                     self.on_sort_key_changed)

        ## Distance box
        box = OWGUI.widgetBox(self.controlArea, "Distance Measure")
        OWGUI.comboBox(box, self, "selected_distance_index",
                       items=[t[0] for t in self.DISTANCE_FUNCTIONS],
                       callback=self.on_distance_measure_changed)

        self.connect(self.graphButton,
                     SIGNAL("clicked()"),
                     self.save_graph)
        
        self.scene = QGraphicsScene()
        self.scene_view = QualityGraphicsView(self.scene)
        self.scene_view.setRenderHints(QPainter.Antialiasing)
        self.mainArea.layout().addWidget(self.scene_view)
        
        self.connect(self.scene_view,
                     SIGNAL("view_size_changed(QSize)"),
                     self.on_view_resize)

        self._disable_updates = False
        self.main_widget = None
        
        self.resize(800, 600)

    def clear(self):
        """Clear the widget state.
        """
        self.data = None
        self.distances = None
        self.groups = None
        self.unique_pos = None
        
        with disable_updates(self):
            self.split_by_model[:] = []
            self.sort_by_model[:] = []

        self.scene.clear()
        self.info_box.setText("\n")

    def set_data(self, data=None):
        """Set input experiment data.
        """
        self.data = data

    def handleNewSignals(self):
        """Called after all signals have been set. 
        """
        if self.data:
            self.on_new_data()
        else:
            self.closeContext("")
            self.clear()

    def update_label_candidates(self):
        """Update the label candidates selection GUI 
        (Group/Sort By views).
        
        """
        keys = self.get_suitable_keys(self.data)
        with disable_updates(self):
            self.split_by_model[:] = keys
            self.sort_by_model[:] = keys
        
    def get_suitable_keys(self, data):
        """ Return suitable attr label keys from the data where the key has at least
        two unique values in the data.
        
        """
        attrs = [attr.attributes.items() for attr in data.domain.attributes]
        attrs  = reduce(list.__add__, attrs, [])
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
    
    def selected_base_index(self):
        return self.base_index

    def on_new_data(self):
        """We have new data and need to recompute all.
        """
        self.closeContext("")
        
        self.update_label_candidates()
        self.info_box.setText("{} genes \n{} experiments".format(
                                len(self.data), 
                                len(self.data.domain.attributes)
                                )
                              )
        
        keys = self.get_suitable_keys(self.data)
        self.openContext("", keys)
        
        ## Restore saved context
        context = self.currentContexts[""]
        split_by_labels= getattr(context, "split_by_labels", set())
        sort_by_labels = getattr(context, "sort_by_labels", set())
        
        def select(model, selection_model, selected_items):
            all_items = list(model)
            try:
                indices = [all_items.index(item) for item in selected_items]
            except:
                indices = []
            for ind in indices:
                selection_model.select(model.index(ind), QItemSelectionModel.Select)
                
#        self._disable_updates = True
#        try:
        with disable_updates(self):
            select(self.split_by_view.model(),
                   self.split_by_view.selectionModel(),
                   split_by_labels)
            
            select(self.sort_by_view.model(),
                   self.sort_by_view.selectionModel(),
                   sort_by_labels)
#        finally:
#            self._disable_updates = False
            
        self.split_and_update()
        
    def on_split_key_changed(self, *args):
        with control_disable(self):
            if not self._disable_updates:
                context = self.currentContexts[""]
                context.split_by_labels = self.selected_split_by_labels()
                self.split_and_update()
    
    def on_sort_key_changed(self, *args):
        with control_disable(self):
            if not self._disable_updates:
                context = self.currentContexts[""]
                context.sort_by_labels = self.selected_sort_by_labels()
                self.split_and_update()
        
    def on_distance_measure_changed(self):
        self.update_distances()
        self.replot_experiments()
        
    def on_view_resize(self, size):
        if self.main_widget:
            current = self.main_widget.size()
            self.main_widget.resize(size.width() - 2, 
                                    current.height())
            
            self.scene.setSceneRect(self.scene.itemsBoundingRect())
        
    def on_rug_item_clicked(self, item):
        base = item.index
        if base != self.base_index:
            self.base_index = base
            self.split_and_update()
        
    def split_and_update(self):
        """
        Split the data based on the selected sort/split labels
        and update the quality plot.
        
        """
        self.groups, self.unique_pos = \
                exp.separate_by(self.data,
                                self.selected_split_by_labels(),
                                consider=self.selected_sort_by_labels(),
                                add_empty=True)
        
        
        self.groups = sorted(self.groups.items(),
                             key=lambda t: map(float_if_posible, t[0]))
        self.unique_pos = sorted(self.unique_pos.items(),
                                 key=lambda t: map(float_if_posible, t[0]))
        
        pprint(self.groups)
        pprint(self.unique_pos)
        if self.groups:
            # TODO: Check if the groups of base experiment have changed
            self.update_distances()
            self.replot_experiments()

    def update_distances(self):
        """Recompute the experiment distances.
        """
        distance = self.selected_distance()
        base_index = self.base_index
        base_distances = []
        pb = OWGUI.ProgressBar(self, len(self.groups) * \
                               len(self.data.domain.attributes))
        for group, indices in self.groups:
            # Base column of the group
            if indices[base_index] is not None:
                base_vec = exp.linearize(self.data, [indices[base_index]])
                distances = []
                # Compute the distances between base replicate 
                # and all the rest data columns.
                for i in range(len(self.data.domain.attributes)):
                    if i == indices[base_index]:
                        distances.append(0.0)
                    else:
                        vec_i = exp.linearize(self.data, [i])
                        distances.append(distance(base_vec, vec_i))
                    pb.advance()
                    
                base_distances.append(distances)
            else:
                base_distances.append(None)
        pb.finish()
        self.distances = base_distances

    def replot_experiments(self):
        """Replot the whole quality plot
        """
        self.scene.clear()
        labels = []
        ## Base replicate=1 TODO: the index should be set
        base_index = self.base_index
        max_dist = numpy.max(filter(None, self.distances))
        print max_dist
        rug_widgets = []
        
        group_pen = QPen(QColor(0, 0, 0))
        group_pen.setWidth(2)
        group_pen.setCapStyle(Qt.RoundCap)
        background_pen = QPen(QColor(0, 0, 250, 150))
        background_pen.setWidth(1)
        background_pen.setCapStyle(Qt.RoundCap)
        
        main_widget = QualityControlWidget()
        layout = QGraphicsGridLayout()
        split_by = self.selected_split_by_labels()
        sort_by = self.selected_sort_by_labels()
        attributes = self.data.domain.attributes
        if self.data:
            for (group, indices), dist_vec in zip(self.groups, self.distances):
                base = indices[base_index]
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
                            tooltip = sort_label(split_by, attr)
                            tooltip += "\n" + sort_label(sort_by, attr)
                            rug_item.setToolTip(tooltip)
                            rug_item.setToolTip(sort_label(sort_by, attr))
                            rug_item.index = indices.index(i)
                        else:
#                            rug_item = RugItem(dist_vec[i] / max_dist, 0.85)
                            rug_item = ClickableRugItem(dist_vec[i] / max_dist,
                                           0.85, None)#self.on_rug_item_clicked)
                            tooltip = sort_label(split_by, attr)
                            tooltip += "\n" + sort_label(sort_by, attr)
                            rug_item.setToolTip(tooltip)
                            rug_item.setPen(background_pen)
                        rug_items.append(rug_item)
                    
                rug_widget = RugGraphicsWidget()
                rug_widget.set_rug(rug_items)
                
                rug_widgets.append(rug_widget)
                
                label = group_label(self.selected_split_by_labels(), group)
#                label_item = QGraphicsSimpleTextItem(label)
                label_item = GtI(label, main_widget)
                label_item = GraphicsSimpleTextLayoutItem(label_item)
                label_item.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
                labels.append(label_item)
        
        for i, (label, w) in enumerate(zip(labels, rug_widgets)):
            layout.addItem(label, i, 0, Qt.AlignVCenter)
            layout.addItem(w, i, 1)
            layout.setRowMaximumHeight(i, 30)
            
        main_widget.setLayout(layout)
        self.scene.addItem(main_widget)
        main_widget.show()
        self.main_widget = main_widget
        self.rug_widgets = rug_widgets
        self.labels = labels
        self.on_view_resize(self.scene_view.size())
        
    def save_graph(self):
        from OWDlgs import OWChooseImageSizeDlg
        dlg = OWChooseImageSizeDlg(self.scene, parent=self)
        dlg.exec_()


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
        size = self.size()
        height = size.height()
        width = size.width()
        
        for item in self.rug_items:
            offset = (1.0 - item.height) * height / 2.0
            item.setPos(width * item.value, 0)
            item.setLine(0., offset, 0., height - offset)

    def resizeEvent(self, event):
        QGraphicsWidget.resizeEvent(self, event)
        self.update_rug_geometry()

    def setGeometry(self, geom):
        QGraphicsWidget.setGeometry(self, geom)


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


class QualityGraphicsView(QGraphicsView):
    def resizeEvent(self, event):
        QGraphicsView.resizeEvent(self, event)
        self.emit(SIGNAL("view_size_changed(QSize)"),
                  event.size())


class QualityControlWidget(QGraphicsWidget):
    if DEBUG:
        def paint(self, painter, options, widget=0):
            rect =  self.geometry()
            rect.translate(-self.pos())
            painter.drawRect(rect)
            
if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWQualityControl()
#    data = Orange.data.Table("doc:dicty-abc-sample.tab")
    data = Orange.data.Table("doc:pipa.tab")

    w.set_data(data)
    w.show()
    w.handleNewSignals()
    app.exec_()
    w.set_data(None)
    w.handleNewSignals()
    w.saveSettings()
