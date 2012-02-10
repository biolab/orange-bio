"""
<name>Quality control</name>

"""

import sys

import Orange

from OWWidget import *
from OWItemModels import PyListModel
import OWGUI

import obiExperiments as exp
import numpy

from pprint import pprint

def compute_within_group_distances(data, groups, distance=exp.dist_pcorr):
    matrices = []
    for _, indices in groups:
        matrix = Orange.misc.SymMatrix(len(group_i))
        for i, index_i in enumerate(indices):
            for j, index_j in enumerate(indics[i + 1:], i + 1):
                vec_i = exp.linearize(data, [index_i])
                vec_j = exp.linearize(data, [index_j])
                matrix[i, j] = distance(vec_i, vec_j)
        metrices.append(matrix)
    return matrices

def group_label(splits, group):
    """
    Return group label
    """
    labels = ["{}={}".format(split, group) \
              for split, group in zip(splits, groups)]
    return " | ".join(labels)



class OWQualityControl(OWWidget):
    settingsList = []
    DISTANCE_FUNCTIONS = [("Distance from Pearson correlation", exp.dist_pcorr),
                          ("Euclidean distance", exp.dist_eucl),
                          ("Distance from Spearman correlation", exp.dist_spearman)]

    def __init__(self, parent=None, signalManager=None,
                 title="Quality Control"):
        OWWidget.__init__(self, parent, signalManager, title,
                          wantGraph=True)

        self.inputs = [("Experiment Data", Orange.data.Table, self.set_data)]

        ## Setings
        self.selected_distance_index = 0
        self.replicate_id = "replicate"

        ## Attributes
        self.data = None
        self.distances = None
        self.groups = None
        self.group_indices = None
        self.separate_by = None
        self.sort_by = None

        ## GUI
        box = OWGUI.widgetBox(self.controlArea, "Info")
        self.info_box = OWGUI.widgetLabel(box, "\n")

        ## Split By box
        box = OWGUI.widgetBox(self.controlArea, "Split By")
        self.split_by_model = PyListModel()
        self.split_by_view = QListView()
        self.split_by_view.setModel(self.split_by_model)
        box.layout().addWidget(self.split_by_view)

        ## Sort By box
        box = OWGUI.widgetBox(self.controlArea, "Sort By")
        self.sort_by_model = PyListModel()
        self.sort_by_view = QListView()
        self.sort_by_view.setModel(self.sort_by_model)
        box.layout().addWidget(self.sort_by_view)

        ## Distance box
        box = OWGUI.widgetBox(self.controlArea, "Distance Measure")
        OWGUI.comboBox(box, self, "selected_distance_index",
                       items=[t[0] for t in self.DISTANCE_FUNCTIONS])

        self.scene = QGraphicsScene()
        self.scene_view = QGraphicsView(self.scene)
        self.mainArea.layout().addWidget(self.scene_view)

        self.resize(800, 600)

    def clear(self):
        """Clear the widget state.
        """
        self.data = None
        self.distances = None
        self.groups = None
        self.group_indices = None
        self.separate_by = None
        self.sort_by = None

        self.scene.clear()

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
            self.clear()

    def update_label_candidates(self):
        """Update the label candidates selection GUI 
        (Group/Sort By views).
        
        """
        pass

    def selected_split_by_labels(self):
        """Return the current selected split labels.
        """
        return ["tp"]

    def selected_sort_by_labels(self):
        """Return the current selected sort labels
        """
        return ["replicate"]

    def selected_distance(self):
        """Return the selected distance function.
        """
        return self.DISTANCE_FUNCTIONS[self.selected_distance_index][1]

    def on_new_data(self):
        """We have new data and need to recompute all.
        """
        self.update_label_candidates()

        self.separate_by = self.selected_split_by_labels()
        self.sort_by = self.selected_sort_by_labels()

        self.split_and_update()

    def split_and_update(self):
        """
        Split the data based on the selected sort/split labels
        and update the quality plot.
        
        """
        self.groups, self.group_indices = \
                exp.separate_by(self.data, self.separate_by,
                                consider=self.sort_by,
                                add_empty=False)

        self.groups = sorted(self.groups.items())
        self.group_indices = sorted(self.group_indices.items())
        
        pprint(self.groups)
        pprint(self.group_indices)

        self.update_distances()
        self.replot_experiments()

    def update_distances(self):
        """Recompute the experiment distances.
        """
        distance = self.selected_distance()
        base_index = 0
        base_distances = []
        for group, indices in self.groups:
            # Base column of the group
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
            base_distances.append(distances)
            
        self.distances = base_distances
#        distance = self.selected_distance()
#        self.distances = \
#            compute_within_group_distances(self.data, self.groups, distance)
#        self.backgroud_distances = \
#            compute_within_group_distances( \
#                            self.data, 
#                            [(None,  
#                             range(len(self.data.attributes)))],
#                            distance
#                            )
#        print self.distances
#        print self.backgroud_distances

    def replot_experiments(self):
        """Replot the whole quality plot
        """
        self.scene.clear()
        labels = []
        ## Base replicate=1 TODO: the index should be set
        base_index = 0
        max_dist = numpy.max(self.distances)
        rug_widgets = []
        if self.data:
            for (group, indices), dist_vec in zip(self.groups, self.distances):
                base = indices[base_index]
                indices_set = set(indices)
                rug_items = []
                for i in range(len(self.data.domain.attributes)):
                    # Is this a within group distance or background
                    in_group = i in indices_set 
                    rug_item = RugItem(dist_vec[i] / max_dist)
                    rug_items.append(rug_item)
                    
                rug_widget = RugGraphicsWidget()
                rug_widget.set_rug(rug_items)
                rug_widgets.append(rug_widget)
                
        layout = QGraphicsLinearLayout(Qt.Vertical)
        for w in rug_widgets:
            layout.addItem(w)
            
        main_widget = QGraphicsWidget()
        main_widget.setLayout(layout)
        main_widget.resize(600, 600)
        self.scene.addItem(main_widget)
        
                    


class RugGraphicsWidget(QGraphicsWidget):
    def __init__(self, parent=None, rug=None):
        QGraphicsWidget.__init__(self, parent)
        self.rug_items = []
        self.set_rug(rug)

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
            item.setPos(width * item.value, 0)
            item.setLine(0., 0., 0., height)
            

    def resizeEvent(self, event):
        QGraphicsWidget.resizeEvent(self, event)
        self.update_rug_geometry()

    def setGeometry(self, geom):
        QGraphicsWidget.setGeometry(self, geom)


class RugItem(QGraphicsLineItem):
    def __init__(self, value):
        QGraphicsLineItem.__init__(self)
        self.value = value

    def set_height(self, height):
        """Set the height of this item.
        """
        self.height = height



if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWQualityControl()
    data = Orange.data.Table("doc:dicty-abc-sample.tab")

    w.set_data(data)
    w.show()
    w.handleNewSignals()
    app.exec_()

