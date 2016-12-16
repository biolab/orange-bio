import sys
import copy
import operator
from functools import reduce

from collections import defaultdict

import numpy

from AnyQt.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QHeaderView, QListView, QScrollArea,
    QSizePolicy, QSpacerItem, QWidgetItem
)
from AnyQt.QtGui import QColor, QPalette, QStandardItemModel, QStandardItem
from AnyQt.QtCore import Qt, QSize, QItemSelectionModel, QItemSelection

import Orange.data
import Orange.misc
from Orange.preprocess import transformation

from Orange.widgets import widget, gui, settings
from Orange.widgets.utils import itemmodels

from ..utils.group import \
    separate_by, data_type, linearize, dist_pcorr, dist_eucl, dist_spearman

from .utils.settings import SetContextHandler


def clone_attr(attr):
    newvar = attr.copy(compute_value=transformation.Identity(attr))
    newvar.attributes = copy.copy(attr.attributes)
    return newvar


def name_gen(a):
    i = 1
    while True:
        yield a + " {0}".format(i)
        i += 1

missing_name_gen = name_gen("!!missing")
inactive_name_gen = name_gen("!!inactive")


class MyHeaderView(QHeaderView):
    def mouseMoveEvent(self, event):
        event.ignore()

    def wheelEvent(self, event):
        event.ignore()


class OWGenotypeDistances(widget.OWWidget):
    name = "Expression Profile Distances"
    description = ("Compute distances between expression profiles of "
                   "different experimental factors.")
    icon = "../widgets/icons/GenotypeDistances.svg"
    priority = 1050

    inputs = [("Data", Orange.data.Table, "set_data")]
    outputs = [("Distances", Orange.misc.DistMatrix),
               ("Sorted Data", Orange.data.Table)]

    settingsHandler = SetContextHandler()

    separate_keys = settings.ContextSetting({})
    relevant_keys = settings.ContextSetting({})

    distance_measure = settings.Setting(0)
    auto_commit = settings.Setting(False)

    DISTANCE_FUNCTIONS = [
        ("Distance from Pearson correlation", dist_pcorr),
        ("Euclidean distance", dist_eucl),
        ("Distance from Spearman correlation", dist_spearman)
    ]

    def __init__(self, parent=None):
        super().__init__(self, parent)

        self.data = None
        self.partitions = []
        self.matrix = None
        self.split_groups = []
        self._disable_updates = False

        ########
        # GUI
        ########

        box = gui.widgetBox(self.controlArea, "Input")

        self.info_box = gui.widgetLabel(box, "No data on input\n")

        box = gui.widgetBox(self.controlArea, "Separate By",
                            addSpace=True)

        self.separate_view = QListView(
            selectionMode=QListView.MultiSelection
        )
        box.layout().addWidget(self.separate_view)

        box = gui.widgetBox(self.controlArea, "Sort By",
                            addSpace=True)
        self.relevant_view = QListView(
            selectionMode=QListView.MultiSelection)

        box.layout().addWidget(self.relevant_view)

        self.distance_view = gui.comboBox(
            self.controlArea, self, "distance_measure",
            box="Distance Measure",
            items=[name for name, _ in self.DISTANCE_FUNCTIONS])

        gui.rubber(self.controlArea)

        gui.auto_commit(self.controlArea, self, "auto_commit", "Commit")
        self.groups_box = gui.widgetBox(self.mainArea, "Groups")
        self.groups_scroll_area = QScrollArea()
        self.groups_box.layout().addWidget(self.groups_scroll_area)

    def sizeHint(self):
        return QSize(800, 600)

    def clear(self):
        self.data = None
        self.partitions = []
        self.split_groups = []
        self.matrix = None

    def get_suitable_keys(self, data):
        """Return suitable attr label keys from the data where the key has at least
        two unique values in the data.

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

    def set_data(self, data=None):
        """Set the input data table.
        """
        self.closeContext()
        self.clear()
        self.error(0)
        self.warning(0)
        if data and not self.get_suitable_keys(data):
            self.error(0, "Data has no suitable column labels.")
            data = None

        self.data = data

        if data:
            self.info_box.setText("{0} genes\n{1} experiments"
                                  .format(len(data), len(data.domain)))
            self.update_control()
            self.split_data()
        else:
            self.separate_view.setModel(itemmodels.PyListModel([]))
            self.relevant_view.setModel(itemmodels.PyListModel([]))
            self.groups_scroll_area.setWidget(QWidget())
            self.info_box.setText("No data on input.\n")
        self.commit()

    def update_control(self):
        """Update the control area of the widget. Populate the list
        views with keys from attribute labels.
        """
        keys = self.get_suitable_keys(self.data)

        model = itemmodels.PyListModel(keys)
        self.separate_view.setModel(model)
        self.separate_view.selectionModel().selectionChanged.connect(
            self.on_separate_key_changed)

        model = itemmodels.PyListModel(keys)
        self.relevant_view.setModel(model)
        self.relevant_view.selectionModel().selectionChanged.connect(
            self.on_relevant_key_changed)

        self.openContext(keys)

        # Get the selected keys from the open context
        separate_keys = self.separate_keys
        relevant_keys = self.relevant_keys

        def select(model, selection_model, selected_items):
            all_items = list(model)
            try:
                indices = [all_items.index(item) for item in selected_items]
            except:
                indices = []
            selection = QItemSelection()
            for ind in indices:
                index = model.index(ind)
                selection.select(index, index)

            selection_model.select(selection, QItemSelectionModel.Select)

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
            self.separate_keys = self.selected_separeate_by_keys()
            self.split_data()

    def on_relevant_key_changed(self, *args):
        if not self._disable_updates:
            self.relevant_keys = self.selected_relevant_keys()
            self.split_data()

    def selected_separeate_by_keys(self):
        """Return the currently selected separate by keys
        """
        rows = self.separate_view.selectionModel().selectedRows()
        rows = sorted([idx.row() for idx in rows])
        keys = [self.separate_view.model()[row] for row in rows]
        return keys

    def selected_relevant_keys(self):
        """Return the currently selected relevant keys
        """
        rows = self.relevant_view.selectionModel().selectedRows()
        rows = sorted([idx.row() for idx in rows])
        keys = [self.relevant_view.model()[row] for row in rows]
        return keys

    def split_data(self):
        """Split the data and update the Groups widget
        """
        separate_keys = self.selected_separeate_by_keys()
        relevant_keys = self.selected_relevant_keys()

        self.warning(0)
        if not separate_keys:
            self.warning(0, "No separate by column selected.")

        partitions, uniquepos = separate_by(
            self.data, separate_keys, consider=relevant_keys)
        partitions = partitions.items()

        all_values = defaultdict(set)
        for a in [at.attributes for at in self.data.domain.attributes]:
            for k, v in a.items():
                all_values[k].add(v)

        # sort groups
        pkeys = [key for key, _ in partitions]
        types = [data_type([a[i] for a in pkeys])
                 for i in range(len(pkeys[0]))]

        partitions = sorted(partitions, key=lambda x:
                    tuple(types[i](v) for i,v in enumerate(x[0])))

        split_groups = []

        # Collect relevant key value pairs for all columns
        relevant_items = None

        for keys, indices in partitions:
            if relevant_items == None:
                relevant_items = [defaultdict(set) for _ in indices]
            for i, ind in enumerate(indices):
                if ind is not None:
                    attr = self.data.domain[ind]
                    for key in relevant_keys:
                        relevant_items[i][key].add(attr.attributes[key])

        #those with different values between rows are not relevant
        for d in relevant_items:
            for k, s in list(d.items()):
                if len(s) > 1:
                    del d[k]
                else:
                    d[k] = s.pop()

        def get_attr(attr_index, i):
            if attr_index is None:
                attr = Orange.data.ContinuousVariable(next(missing_name_gen), 
                    compute_value=lambda x: None)
                attr.attributes.update(relevant_items[i])
                return attr
            else:
                return self.data.domain[attr_index]

        for keys, indices in partitions:
            attrs = [get_attr(attr_index, i)
                     for i, attr_index in enumerate(indices)]
            for attr in attrs:
                attr.attributes.update(zip(separate_keys, keys))
            domain = Orange.data.Domain(attrs, [], self.data.domain.metas)
            split_groups.append((keys, domain))

        self.set_groups(separate_keys, split_groups, relevant_keys,
                        relevant_items, all_values, uniquepos)

        self.partitions = partitions
        self.split_groups = split_groups

        self.commit()

    def set_groups(self, keys, groups, relevant_keys, relevant_items,
                   all_values, uniquepos):
        """Set the current data groups and update the Group widget
        """
        layout = QVBoxLayout()
        header_widths = []
        header_views = []
        palette = self.palette()
        all_values = all_values.keys()

        def for_print(rd):
            attrs = []
            for d in rd:
                attr = Orange.data.ContinuousVariable(next(inactive_name_gen))
                attr.attributes.update(d)
                attrs.append(attr)
            return Orange.data.Domain(attrs, None)

        for separatev, domain in [(None, for_print(relevant_items))] + groups:
            label = None
            if separatev is not None:
                ann_vals = " <b>|</b> ".join(["<b>{0}</b> = {1}".format(key,val) \
                     for key, val in zip(keys, separatev)])
                label = QLabel(ann_vals)

            model = QStandardItemModel()
            for i, attr in enumerate(domain.attributes):
                item = QStandardItem()
                if separatev is not None:
                    isunique = uniquepos[separatev][i]
                else:
                    isunique = all(a[i] for a in uniquepos.values())

                if str(attr.name).startswith("!!missing "):  # TODO: Change this to not depend on name
                    header_text = ["{0}={1}".format(key, attr.attributes.get(key, "?")) \
                                   for key in all_values if key not in relevant_items[i]]
                    header_text = "\n".join(header_text) if header_text else "Empty"
                    item.setData(header_text, Qt.DisplayRole)
                    item.setFlags(Qt.NoItemFlags)
                    item.setData(QColor(Qt.red), Qt.ForegroundRole)
                    item.setData(palette.color(QPalette.Disabled, QPalette.Window), Qt.BackgroundRole)
                    item.setData("Missing feature.", Qt.ToolTipRole)
                elif str(attr.name).startswith("!!inactive "):
                    header_text = ["{0}={1}".format(key, attr.attributes.get(key, "?")) \
                                   for key in all_values if key in relevant_items[i]]
                    header_text = "\n".join(header_text) if header_text else "No descriptor"
                    item.setData(header_text, Qt.DisplayRole)
                    item.setData(palette.color(QPalette.Disabled, QPalette.Window), Qt.BackgroundRole)
                else:
                    header_text = ["{0}={1}".format(key, attr.attributes.get(key, "?")) \
                                   for key in all_values if key not in relevant_items[i]]
                    header_text = "\n".join(header_text) if header_text else "Empty"
                    item.setData(header_text, Qt.DisplayRole)
                    item.setData(attr.name, Qt.ToolTipRole)

                if not isunique:
                    item.setData(QColor(Qt.red), Qt.ForegroundRole)

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

        left, _, right, _ = self.getContentsMargins()
        widget.setMinimumWidth(width_sum)
        widget.setMinimumWidth(max_width + left + right)
        self.groups_scroll_area.setWidget(widget)

    def compute_distances(self, separate_keys, partitions, data):
        """Compute the distances between genotypes.
        """
        if separate_keys and partitions:
            self.progressBarInit()
#             matrix = Orange.misc.DistMatrix(len(partitions))
            matrix = numpy.zeros((len(partitions), len(partitions)))

            profiles = [linearize(data, indices) for _, indices in partitions]
            dist_func = self.DISTANCE_FUNCTIONS[self.distance_measure][1]
#             from Orange.utils import progress_bar_milestones
            count = (len(profiles) * len(profiles) - 1) / 2
#             milestones = progress_bar_milestones(count)
            iter_count = 0
            for i in range(len(profiles)):
                for j in range(i + 1, len(profiles)):
                    matrix[i, j] = dist_func(profiles[i], profiles[j])
                    matrix[j, i] = matrix[i, j]
                    iter_count += 1
#                     if iter_count in milestones:
                    self.progressBarSet(100.0 * iter_count / count)
            self.progressBarFinished()

            items = [["{0}={1}".format(key, value)
                      for key, value in zip(separate_keys, values)]
                     for values, _ in partitions]
            items = [" | ".join(item) for item in items]
#             matrix.setattr("items", items)
            matrix = Orange.misc.DistMatrix(matrix)
        else:
            matrix = None

        self.matrix = matrix

    def commit(self):
        separate_keys = self.selected_separeate_by_keys()
        self.compute_distances(separate_keys,
                               self.partitions,
                               self.data)
        if self.split_groups:
            all_attrs = []
            for group, domain in self.split_groups:
                attrs = []
                group_name = " | ".join("{0}={1}".format(*item) for item in
                                        zip(separate_keys, group))
                for attr in domain.attributes:
                    newattr = clone_attr(attr)
                    newattr.attributes["<GENOTYPE GROUP>"] = group_name # Need a better way to pass the groups to downstream widgets.
                    attrs.append(newattr)

                all_attrs.extend(attrs)

            domain = Orange.data.Domain(all_attrs, self.data.domain.class_vars,
                                        self.data.domain.metas)
            data = Orange.data.Table(domain, self.data)
        else:
            data = None
        self.send("Sorted Data", data)
        self.send("Distances", self.matrix)

    def send_report(self):
        self.report_items((
            ("Separate By", ", ".join(self.selected_separeate_by_keys())),
            ("Sort By", ", ".join(self.selected_relevant_keys())),
            ("Distance Measure", self.DISTANCE_FUNCTIONS[self.distance_measure][0])
        ))
        layout = self.groups_scroll_area.widget().layout()
        html = "<table>"
        for i in range(layout.count()):
            item = layout.itemAt(i)
            if isinstance(item, QSpacerItem):
                html += "<tr><td></td></tr>"
            elif isinstance(item, QWidgetItem):
                hor = item.widget()
                if isinstance(hor, QLabel):
                    label = hor.text()
                    html += "<tr><td><b>%s</b></td></tr>" % label
                elif isinstance(hor, QHeaderView):
                    model = hor.model()
                    content = (model.horizontalHeaderItem(col) for col in range(model.columnCount()))
                    content = (item.text().replace('\n', "<br/>") for item in content)
                    html += "<tr>" + ''.join("<td>{}</td>".format(item) for item in content) + "</tr>"
        html += "</table>"
        self.report_raw("Groups", html)


def test_main(argv=sys.argv):
    import os
    from AnyQt.QtWidgets import QApplication

    if len(argv) > 1:
        filename = argv[1]
    else:
        filename = os.path.expanduser("~/GDS624.tab")
    data = Orange.data.Table(filename)

    app = QApplication(sys.argv)
    w = OWGenotypeDistances()
    w.set_data(data)
    w.show()
    w.raise_()
    r = app.exec_()
    w.saveSettings()
    return r

if __name__ == "__main__":
    sys.exit(test_main())
