"""<name>Genotype Distances</name>
<description>Compute distances between expression profiles of different experimental factors.</description>
<icon>icons/GenotypeDistances.svg</icon>
<priority>1050</priority>
<contact>Ales Erjavec (ales.erjavec(@at@).fri.uni-lj.si)</contact>
"""

from __future__ import absolute_import

from collections import defaultdict
from operator import add
import math

import numpy

import Orange
from Orange.OrangeWidgets import OWGUI
from Orange.OrangeWidgets.OWItemModels import PyListModel
from Orange.OrangeWidgets.OWWidget import *

from ..obiExperiments import separate_by, data_type, linearize, dist_pcorr, dist_eucl, dist_spearman

NAME = "Expression Profile Distances"
DESCRIPTION = "Compute distances between expression profiles of different experimental factors."
ICON = "icons/GenotypeDistances.svg"
PRIORITY = 1050

INPUTS = [("Example Table", Orange.data.Table, "set_data")]
OUTPUTS = [("Distances", Orange.core.SymMatrix),
           ("Sorted Example Table", Orange.data.Table)]

REPLACES = ["_bioinformatics.widgets.OWGenotypeDistances.OWGenotypeDistances"]


def clone_attr(attr):
    newattr = attr.clone()
    def get_value_from(ex, w=None):
        if attr in ex.domain:
            v = str(ex[attr])
        else:
            v = "?"
        return newattr(v)
        
    newattr.get_value_from = get_value_from
    newattr.source_variable = attr
    newattr.attributes.update(attr.attributes)
    return newattr
    
def name_gen(a):
    i = 1
    while True:
        yield a + " {0}".format(i)
        i += 1
        
missing_name_gen = name_gen("!!missing")
inactive_name_gen = name_gen("!!inactive")


class MyHeaderView(QHeaderView):
    def __init__(self, *args):
        QHeaderView.__init__(self, *args)
        
    def mouseMoveEvent(self, event):
        event.ignore()
        
    def wheelEvent(self, event):
        event.ignore()
        
        
class SetContextHandler(ContextHandler):
    def match(self, context, imperfect, items):
        items = set(items)
        saved_items = set(getattr(context, "items", []))
        if imperfect:
            return float(len(items.intersection(saved_items)))/(len(items.union(saved_items)) or 1)
        else:
            return items == saved_items
        
    def findOrCreateContext(self, widget, items):        
        index, context, score = self.findMatch(widget, self.findImperfect, items)
        if context:
            if index < 0:
                self.addContext(widget, context)
            else:
                self.moveContextUp(widget, index)
            return context, False
        else:
            context = self.newContext()
            context.items = items
            self.addContext(widget, context)
            return context, True
        
class OWGenotypeDistances(OWWidget):
    contextHandlers = {"": SetContextHandler("")}
    settingsList = ["auto_commit"]
    
    DISTANCE_FUNCTIONS = [("Distance from Pearson correlation", dist_pcorr),
                          ("Euclidean distance", dist_eucl),
                          ("Distance from Spearman correlation", dist_spearman)]
    
    def __init__(self, parent=None, signalManager=None, title="Expression Profile Distances"):
        OWWidget.__init__(self, parent, signalManager, title)
        
        self.inputs = [("Example Table", ExampleTable, self.set_data)]
        self.outputs = [("Distances", Orange.core.SymMatrix), ("Sorted Example Table", ExampleTable)]
        
        self.distance_measure = 0
        self.auto_commit = False
        self.changed_flag = False
        
        self.loadSettings()
        
        ########
        # GUI
        ########
        
        self.info_box = OWGUI.widgetLabel(OWGUI.widgetBox(self.controlArea, "Input",
                                                         addSpace=True),
                                         "No data on input\n")
        
        box = OWGUI.widgetBox(self.controlArea, "Separate By",
                              addSpace=True)
        self.separate_view = QListView()
        self.separate_view.setSelectionMode(QListView.MultiSelection)
        box.layout().addWidget(self.separate_view)
        
        box = OWGUI.widgetBox(self.controlArea, "Sort By",
                              addSpace=True)
        self.relevant_view = QListView()
        self.relevant_view.setSelectionMode (QListView.MultiSelection)
        box.layout().addWidget(self.relevant_view)
        
        self.distance_view = OWGUI.comboBox(self.controlArea, self, "distance_measure",
                                            box="Distance Measure",
                                            items=[d[0] for d in self.DISTANCE_FUNCTIONS])
        
        OWGUI.rubber(self.controlArea)
        
        box = OWGUI.widgetBox(self.controlArea, "Commit")
        cb = OWGUI.checkBox(box, self, "auto_commit", "Commit on any change",
                            tooltip="Compute and send the distances on any change.",
                            callback=self.commit_if)
        
        b = OWGUI.button(box, self, "Commit",
                         tooltip="Compute the distances and send the output signals.",
                         callback=self.commit,
                         default=True)
        
        OWGUI.setStopper(self, b, cb, "changed_flag", callback=self.commit)
        
        self.groups_box = OWGUI.widgetBox(self.mainArea, "Groups")
        self.groups_scroll_area = QScrollArea()
        self.groups_box.layout().addWidget(self.groups_scroll_area)
        
        self.data = None
        self.partitions = []
        self.matrix = None
        self.split_groups = []
        self._disable_updates = False
        
        self.resize(800, 600)
        
    def clear(self):
        self.data = None
        self.partitions = []
        self.split_groups = []
        self.matrix = None
        self.send("Distances", None)
        
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
        
    def set_data(self, data=None):
        """ Set the input example table.
        """
        self.closeContext()
        self.clear()
        self.data = data
        self.error(0)
        self.warning(0)
        if data and not self.get_suitable_keys(data):
            self.error(0, "Data has no suitable attribute labels.")
            data = None
            
        if data:
            self.info_box.setText("{0} genes\n{1} experiments".format(len(data), len(data.domain)))
            self.update_control()
            self.split_data()
        else:
            self.separate_view.setModel(PyListModel([]))
            self.relevant_view.setModel(PyListModel([]))
            self.groups_scroll_area.setWidget(QWidget())
            self.info_box.setText("No data on input.\n")
            self.commit()
            
    def update_control(self):
        """ Update the control area of the widget. Populate the list
        views with keys from attribute labels.
        
        """
        keys = self.get_suitable_keys(self.data)
         
        model = PyListModel(keys)
        self.separate_view.setModel(model)
        self.connect(self.separate_view.selectionModel(),
                     SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),
                     self.on_separate_key_changed)
        
        model = PyListModel(keys)
        self.relevant_view.setModel(model)
        self.connect(self.relevant_view.selectionModel(),
                     SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),
                     self.on_relevant_key_changed)
        
        self.openContext("", keys)
        
        # Get the selected keys from the open context
        context = self.currentContexts[""]
        separate_keys = getattr(context, "separate_keys", set())
        relevant_keys = getattr(context, "relevant_keys", set())
        
        def select(model, selection_model, selected_items):
            all_items = list(model)
            try:
                indices = [all_items.index(item) for item in selected_items]
            except:
                indices = []
            for ind in indices:
                selection_model.select(model.index(ind), QItemSelectionModel.Select)
                
        self._disable_updates = True
        try:
            select(self.relevant_view.model(),
                   self.relevant_view.selectionModel(),
                   relevant_keys)
            
            select(self.separate_view.model(),
                   self.separate_view.selectionModel(),
                   separate_keys)
        finally:
            self._disable_updates = False
        
    def on_separate_key_changed(self, *args):
        if not self._disable_updates:
            context = self.currentContexts[""]
            context.separate_keys = self.selected_separeate_by_keys()
            self.split_data()
    
    def on_relevant_key_changed(self, *args):
        if not self._disable_updates:
            context = self.currentContexts[""]
            context.relevant_keys = self.selected_relevant_keys()
            self.split_data()
        
    def selected_separeate_by_keys(self):
        """ Return the currently selected separate by keys
        """
        rows = self.separate_view.selectionModel().selectedRows()
        rows = sorted([idx.row() for idx in rows])
        keys = [self.separate_view.model()[row] for row in rows]
        return keys
        
    def selected_relevant_keys(self):
        """ Return the currently selected relevant keys
        """
        rows = self.relevant_view.selectionModel().selectedRows()
        rows = sorted([idx.row() for idx in rows])
        keys = [self.relevant_view.model()[row] for row in rows]
        return keys
    
    def split_data(self):
        """ Split the data and update the Groups widget
        """
        separate_keys = self.selected_separeate_by_keys()
        relevant_keys = self.selected_relevant_keys()
        
        self.warning(0)
        if not separate_keys:
            self.warning(0, "No separate by attribute selected.")

        partitions,uniquepos = separate_by(self.data, separate_keys, consider=relevant_keys)
        partitions = partitions.items()

        all_values = defaultdict(set)
        for a in [ at.attributes for at in self.data.domain.attributes ]:
            for k,v in a.iteritems():
                all_values[k].add(v)

        #sort groups
        pkeys = [ key for key,_ in partitions ]
        types = [ data_type([a[i] for a in pkeys]) for i in range(len(pkeys[0])) ]

        partitions = sorted(partitions, key=lambda x:
                    tuple(types[i](v) for i,v in enumerate(x[0])))

        split_groups = []
        
        # Collect relevant key value pairs for all columns
        relevant_items = None

        for keys, indices in partitions:
            if relevant_items == None:
                relevant_items = [ defaultdict(set) for _ in range(len(indices)) ]
            for i, ind in enumerate(indices):
                if ind is not None:
                    attr = self.data.domain[ind]
                    for key in relevant_keys:
                        relevant_items[i][key].add(attr.attributes[key])

        #those with different values between rows are not relevant
        for d in relevant_items:
            for k,s in d.items():
                if len(s) > 1:
                    del d[k]
                else:
                    d[k] = s.pop()

        def get_attr(attr_index, i):
            if attr_index is None:
                attr = Orange.feature.Continuous(missing_name_gen.next())
                attr.attributes.update(relevant_items[i])
                return attr
            else:
                return self.data.domain[attr_index]
        
        for keys, indices in partitions:
            attrs = [get_attr(attr_index, i) for i, attr_index in enumerate(indices)]
            for attr in attrs:
                attr.attributes.update(zip(separate_keys, keys))
            domain = Orange.data.Domain(attrs, None)
            domain.add_metas(self.data.domain.get_metas().items())
#            newdata = Orange.data.Table(domain)
            split_groups.append((keys, domain))
         
        self.set_groups(separate_keys, split_groups, relevant_keys, relevant_items, all_values, uniquepos)
        
        self.partitions = partitions
        self.split_groups = split_groups
        
        self.commit_if()
#        self.update_distances(separate_keys, partitions, self.data)
        
    def set_groups(self, keys, groups, relevant_keys, relevant_items, all_values, uniquepos):
        """ Set the current data groups and update the Group widget
        """
        layout = QVBoxLayout()
        header_widths = []
        header_views = []
        palette = self.palette()
        all_values = all_values.keys()

        def for_print(rd):
            attrs = []
            for d in rd:
                attr = Orange.feature.Continuous(inactive_name_gen.next())
                attr.attributes.update(d)
                attrs.append(attr)
            return Orange.data.Domain(attrs, None)

        for separatev, domain in [ (None, for_print(relevant_items)) ] + groups:
            label = None
            if separatev != None:
                ann_vals = " <b>|</b> ".join(["<b>{0}</ b> = {1}".format(key,val) \
                     for key, val in zip(keys, separatev)])
                label = QLabel(ann_vals)
            
            model = QStandardItemModel()
            for i, attr in enumerate(domain.attributes):
                item = QStandardItem()
                if separatev != None:
                    up = uniquepos[separatev][i]
                else:
                    up = False if False in [ a[i] for a in uniquepos.values() ] else True
                if str(attr.name).startswith("!!missing "): ## TODO: Change this to not depend on name
                    header_text = ["{0}={1}".format(key, attr.attributes.get(key, "?")) \
                                   for key in all_values if key not in relevant_items[i]]
                    header_text = "\n".join(header_text) if header_text else "Empty"
                    item.setData(QVariant(header_text), Qt.DisplayRole)
                    item.setFlags(Qt.NoItemFlags)
                    item.setData(QVariant(QColor(Qt.red)), Qt.ForegroundRole)
                    item.setData(QVariant(palette.color(QPalette.Disabled, QPalette.Window)), Qt.BackgroundRole)
                    item.setData(QVariant("Missing feature."), Qt.ToolTipRole)
                elif str(attr.name).startswith("!!inactive "):
                    header_text = ["{0}={1}".format(key, attr.attributes.get(key, "?")) \
                                   for key in all_values if key in relevant_items[i]]
                    header_text = "\n".join(header_text) if header_text else "No descriptor"
                    item.setData(QVariant(header_text), Qt.DisplayRole)
                    item.setData(QVariant(palette.color(QPalette.Disabled, QPalette.Window)), Qt.BackgroundRole)
                else:
                    header_text = ["{0}={1}".format(key, attr.attributes.get(key, "?")) \
                                   for key in all_values if key not in relevant_items[i]]
                    header_text = "\n".join(header_text) if header_text else "Empty"
                    item.setData(QVariant(header_text), Qt.DisplayRole)
                    item.setData(QVariant(attr.name), Qt.ToolTipRole)

                if up == False:
                    item.setData(QVariant(QColor(Qt.red)), Qt.ForegroundRole)
                    
                model.setHorizontalHeaderItem(i, item)
            attr_count = len(domain.attributes)
            view = MyHeaderView(Qt.Horizontal)
            view.setResizeMode(QHeaderView.Fixed)
            view.setModel(model)
            hint = view.sizeHint()
            view.setMaximumHeight(hint.height())
            
            widths = [view.sectionSizeHint(i) for i in range(attr_count)]
            header_widths.append(widths)
            header_views.append(view)
            
            if label:
                layout.addWidget(label)
            layout.addWidget(view)
            layout.addSpacing(8)
            
        # Make all header sections the same width
        width_sum = 0
        max_header_count = max([h.count() for h in header_views])
        for i in range(max_header_count):
            max_width = max([w[i] for w in header_widths if i < len(w)] or [0])
            for view in header_views:
                if i < view.count():
                    view.resizeSection(i, max_width)
            width_sum += max_width + 2
                
        for h in header_views:
            h.setMinimumWidth(h.length() + 4)
            
        widget = QWidget()
        widget.setLayout(layout)
        widget.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Maximum)
        layout.activate()
        
        max_width = max(h.length() for h in header_views) + 20
        
        left, top, right, bottom  = self.getContentsMargins()
        widget.setMinimumWidth(width_sum)
        widget.setMinimumWidth(max_width + left + right)
        self.groups_scroll_area.setWidget(widget)
        
    def compute_distances(self, separate_keys, partitions, data):
        """ Compute the distances between genotypes.
        """
        if separate_keys and partitions:
            self.progressBarInit()
            matrix = Orange.core.SymMatrix(len(partitions))
            profiles = [linearize(data, indices) for _, indices in partitions]
            dist_func = self.DISTANCE_FUNCTIONS[self.distance_measure][1]
            from Orange.utils import progress_bar_milestones
            count = (len(profiles) * len(profiles) - 1) / 2
            milestones = progress_bar_milestones(count)
            iter_count = 0
            for i in range(len(profiles)):
                for j in range(i + 1, len(profiles)):
                    matrix[i, j] = dist_func(profiles[i], profiles[j])
                    iter_count += 1
                    if iter_count in milestones:
                        self.progressBarSet(100.0 * iter_count / count)
            self.progressBarFinished()
            
            items = [["{0}={1}".format(key, value) for key, value in zip(separate_keys, values)] \
                      for values, _ in partitions]
            items = [" | ".join(item) for item in items]
            matrix.setattr("items", items)
        else:
            matrix = None
            
        self.matrix = matrix
        
    def commit_if(self):
        if self.auto_commit and self.changed_flag:
            self.commit()
        else:
            self.changed_flag = True
            
    def commit(self):
        separate_keys = self.selected_separeate_by_keys()
        self.compute_distances(separate_keys,
                               self.partitions,
                               self.data)
        
        if self.split_groups:
            all_attrs = []
            for group, domain in self.split_groups: 
                attrs = []
                group_name = " | ".join("{0}={1}".format(*item) for item in \
                                        zip(separate_keys, group))
                for attr in domain.attributes:
                    newattr = clone_attr(attr)
                    newattr.attributes["<GENOTYPE GROUP>"] = group_name # Need a better way to pass the groups to downstream widgets.
                    attrs.append(newattr)
                    
                all_attrs.extend(attrs)
                
            #all_attrs = reduce(add, [list(domain.attributes) for _, domain in self.split_groups], [])
            domain = Orange.data.Domain(all_attrs, self.data.domain.class_var)
            domain.add_metas(self.data.domain.get_metas().items())
            
            data = Orange.data.Table(domain, self.data)
        else:
            data = None
        self.send("Sorted Example Table", data)
        self.send("Distances", self.matrix)
        self.changed_flag = False

if __name__ == "__main__":
    import os, sys
    app = QApplication(sys.argv )
    w = OWGenotypeDistances()
#    data = Orange.data.Table(os.path.expanduser("~/Documents/dicty-express-sample.tab"))
#    data = Orange.data.Table(os.path.expanduser("~/Downloads/tmp.tab"))
    data = Orange.data.Table(os.path.expanduser("~/tgr.tab"))
    w.set_data(data)
    w.show()
    app.exec_()
    w.saveSettings()
    
#    data = Orange.data.Table("tmp.tab")
#    partitions = separate_by(data, [ "genotype" ], consider=["tp", "replicate"]).items()
#    print partitions
#    l1 = linearize(data, partitions[0][1])
#    l2 = linearize(data, partitions[1][1])
#    print  dist_eucl(l1, l2)
#    print  dist_pcorr(l1, l2)

    

