import sys
import warnings

from xml.sax.saxutils import escape
from contextlib import contextmanager
from itertools import starmap, chain
from collections import namedtuple, defaultdict

if sys.version_info < (3, 4):
    from singledispatch import singledispatch
else:
    from functools import singledispatch

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import Qt, pyqtSignal as Signal

import Orange.data
from Orange.widgets.utils import itemmodels
from Orange.widgets import gui as _gui


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


@singledispatch
def standarditem_from(obj):
    item = QtGui.QStandardItem(str(obj))
    item.setData(obj, Qt.UserRole)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item


@standarditem_from.register(Orange.data.Variable)
def standarditem_from_var(var):
    item = QtGui.QStandardItem(var.name)
    _, icon = _gui.attributeItem(var)
    item.setIcon(icon)
    item.setToolTip(itemmodels.Variable)


def variable_tooltip(var):
    if isinstance(var, Orange.data.DiscreteVariable):
        return discrete_variable_tooltip(var)
    elif isinstance(Orange.data.ContinuousVariable):
        return continuous_variable_toltip(var)
    elif isinstance(var, Orange.data.StringVariable):
        return string_variable_tooltip(var)


def variable_labels_tooltip(var):
    text = ""
    if var.attributes:
        items = [(escape(key), escape(value))
                 for key, value in var.attributes.items()]
        labels = list(starmap("{!s} = {!s}".format, items))
        text += "<br/>Variable Labels:<br/>"
        text += "<br/>".join(labels)
    return text


def discrete_variable_tooltip(var):
    text = "<b>%s</b><br/>Discrete with %i values: " %\
           (escape(var.name), len(var.values))
    text += ", ".join("%r" % escape(v) for v in var.values)
    text += variable_labels_tooltip(var)
    return text


def continuous_variable_toltip(var):
    text = "<b>%s</b><br/>Continuous" % escape(var.name)
    text += variable_labels_tooltip(var)
    return text


def string_variable_tooltip(var):
    text = "<b>%s</b><br/>String" % escape(var.name)
    text += variable_labels_tooltip(var)
    return text


def python_variable_tooltip(var):
    text = "<b>%s</b><br/>Python" % escape(var.name)
    text += variable_labels_tooltip(var)
    return text


def group_model(rowgroups, columngroups):
    model = QtGui.QStandardItemModel()
    for key, values in rowgroups:
        item = standarditem_from(key)
        values = [standarditem_from(value) for value in values]
        item.appendRows(values)
        model.addRow(item)


class LabelSelectionWidget(QtGui.QWidget):
    """
    A widget for selection of label values.

    The widget displays the contents of the model with two widgets:

    * The top level model items are displayed in a combo box.
    * The children of the current selected top level item are
      displayed in a subordinate list view.

    .. note:: This is not a QAbstractItemView subclass.
    """
    #: Current group/root index has changed.
    groupChanged = Signal(int)
    #: Selection for the current group/root has changed.
    groupSelectionChanged = Signal()

    def __init__(self, parent=None, **kwargs):
        super().__init__(parent, **kwargs)
        self.__model = None
        self.__selectionMode = QtGui.QListView.ExtendedSelection

        self.__currentIndex = -1
        self.__selections = {}

        layout = QtGui.QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)

        def group_box(title):
            box = QtGui.QGroupBox(title)
            box.setFlat(True)
            lay = QtGui.QVBoxLayout()
            lay.setContentsMargins(0, 0, 0, 0)
            box.setLayout(lay)
            return box

        self.labels_combo = QtGui.QComboBox()
        self.values_view = QtGui.QListView(
            selectionMode=self.__selectionMode
        )

        self.labels_combo.currentIndexChanged.connect(
            self.__onCurrentIndexChanged)

        l_box = group_box(self.tr("Label"))
        v_box = group_box(self.tr("Values"))

        l_box.layout().addWidget(self.labels_combo)
        v_box.layout().addWidget(self.values_view)

        layout.addWidget(l_box)
        layout.addWidget(v_box)

        self.setLayout(layout)

        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                           QtGui.QSizePolicy.Expanding)

    def clear(self):
        """
        Clear the widget/model (same as ``setModel(None)``).
        """
        if self.__model is not None:
            self.values_view.selectionModel().clearSelection()
            self.values_view.selectionModel().selectionChanged.disconnect(
                self.__onSelectionChanged)

            self.values_view.setModel(None)
            self.labels_combo.setModel(
                QtGui.QStandardItemModel(self.labels_combo))
            self.__currentIndex = -1
            self.__selections = {}
            self.__model = None

    def setModel(self, model):
        """
        Set the source model for display.

        The model should be a tree model with depth 2.

        Parameters
        ----------
        model : QtCore.QAbstractItemModel
        """
        if model is self.__model:
            return

        self.clear()

        if model is None:
            return

        self.__model = model

        self.values_view.setModel(model)
        self.values_view.setRootIndex(model.index(0, 0))
        self.values_view.selectionModel().selectionChanged.connect(
            self.__onSelectionChanged)

        # will emit the currentIndexChanged (if the model is not empty)
        self.labels_combo.setModel(model)

    def model(self):
        """
        Return the current model.

        Returns
        -------
        model : QtCore.QAbstractItemModel
        """
        return self.__model

    def setCurrentGroupIndex(self, index):
        """
        Set the current selected group/root index row.

        Parameters
        ----------
        index : int
            Group index.
        """
        self.labels_combo.setCurrentIndex(index)

    def currentGroupIndex(self):
        """
        Return the current selected group/root index row.

        Returns
        -------
        row : index
            Current group row index (-1 if there is no current index)
        """
        return self.labels_combo.currentIndex()

    def setSelection(self, selection):
        """
        Set the model item selection.

        Parameters
        ----------
        selection : QtGui.QItemSelection
            Item selection.
        """
        if self.values_view.selectionModel() is not None:
            indices = selection.indexes()
            pind = defaultdict(list)

            for index in indices:
                parent = index.parent()
                if parent.isValid():
                    if parent == self.__model.index(parent.row(),
                                                    parent.column()):
                        pind[parent.row()].append(
                            QtCore.QPersistentModelIndex(index))
                    else:
                        warnings.warn("Die Die Die")
                else:
                    # top level index
                    pass

            self.__selections = pind
            self.__restoreSelection()

    def selection(self):
        """
        Return the item selection.

        Returns
        -------
        selection : QtGui.QItemSelection
        """
        selection = QtGui.QItemSelection()
        if self.__model is None:
            return selection

        for pind in chain(*self.__selections.values()):
            ind = self.__model.index(pind.row(), pind.column(), pind.parent())
            if ind.isValid():
                selection.select(ind, ind)
        return selection

    def currentGroupSelection(self):
        """
        Return the item selection for the current group only.
        """
        if self.values_view.selectionModel() is not None:
            return self.values_view.selectionModel().selection()
        else:
            return QtGui.QItemSelection()

    def __onCurrentIndexChanged(self, index):
        self.__storeSelection(self.__currentIndex,
                              self.values_view.selectedIndexes())

        self.__currentIndex = index
        if self.__model is not None:
            root = self.__model.index(index, 0)
            self.values_view.setRootIndex(root)

            self.__restoreSelection()
        self.groupChanged.emit(index)

    def __onSelectionChanged(self, old, new):
        self.__storeSelection(self.__currentIndex,
                              self.values_view.selectedIndexes())

        self.groupSelectionChanged.emit()

    def __storeSelection(self, groupind, indices):
        # Store current values selection for the current group
        groupind = self.__currentIndex
        indices = [QtCore.QPersistentModelIndex(ind)
                   for ind in self.values_view.selectedIndexes()]
        self.__selections[groupind] = indices

    def __restoreSelection(self):
        # Restore previous selection for root (if available)
        assert self.__model is not None
        groupind = self.__currentIndex
        root = self.__model.index(groupind, 0)
        sel = self.__selections.get(groupind, [])
        indices = [self.__model.index(pind.row(), pind.column(), root)
                   for pind in sel if pind.isValid() and pind.parent() == root]

        selection = QtGui.QItemSelection()
        for ind in indices:
            selection.select(ind, ind)
        self.values_view.selectionModel().select(
            selection, QtGui.QItemSelectionModel.ClearAndSelect)

    def sizeHint(self):
        """Reimplemented from QWidget.sizeHint"""
        return QtCore.QSize(100, 200)


def itemselection(modelindexlist):
    """
    Return an QtCore.QItemSelection from QModelIndex list

    Parameters
    ----------
    modelindexlist : list of QtCore.QModelIndex
        Selected model indices.

    Returns
    -------
    selection : QtCore.QItemSelection
    """
    selection = QtGui.QItemSelection()
    for index in modelindexlist:
        selection.select(index, index)
    return selection
