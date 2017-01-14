import itertools
from collections import defaultdict, namedtuple
from xml.sax.saxutils import escape

import numpy

from AnyQt import QtGui
from AnyQt.QtCore import Qt

import Orange.data
from Orange.widgets import gui

from . import gui as guiutils

ColumnGroup = namedtuple("ColumnGroup", ["name", "key", "values"])
RowGroup = namedtuple("RowGroup", ["name", "var", "values"])


def group_candidates(data):
    items = [attr.attributes.items() for attr in data.domain.attributes]
    items = list(itertools.chain(*items))

    targets = defaultdict(set)
    for label, value in items:
        targets[label].add(value)

    # Need at least 2 distinct values or key
    targets = [(key, sorted(vals)) for key, vals in targets.items() \
               if len(vals) >= 2]

    column_groups = [ColumnGroup(key, key, values)
                     for key, values in sorted(targets)]
    disc_vars = [var for var in data.domain.class_vars + data.domain.metas
                 if isinstance(var, Orange.data.DiscreteVariable)
                 and len(var.values) >= 2]

    row_groups = [RowGroup(var.name, var, var.values)
                  for var in disc_vars]
    return column_groups, row_groups


@guiutils.standarditem_from.register(ColumnGroup)
def standarditem_from_columngroup(colgroup):
    item = QtGui.QStandardItem(colgroup.name)
#     item.setIcon(pkg_path('columnset.svg'))
    item.setToolTip("Split by column label: '{!s}'"
                    .format(escape(colgroup.name)))
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    item.setData(colgroup, Qt.UserRole)
    children = [guiutils.standarditem_from(val)
                for val in colgroup.values]
    item.appendRows(children)
    return item


@guiutils.standarditem_from.register(RowGroup)
def standarditem_from_rowgroup(rowgroup):
    item = QtGui.QStandardItem(rowgroup.name)
    icon, _ = gui.attributeItem(rowgroup.var)
    item.setIcon(icon)
    item.setToolTip(guiutils.variable_tooltip(rowgroup.var))
    item.setData(rowgroup, Qt.UserRole)
    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    children = [guiutils.standarditem_from(val)
                for val in rowgroup.values]
    item.appendRows(children)
    return item


def group_selection_mask(data, group, indices):
    """
    Return the selection masks for the group.
    """
    if isinstance(group, ColumnGroup):
        selected = [group.values[i] for i in indices]
        target = set([(group.key, value) for value in selected])
        I = [bool(set(var.attributes.items()).intersection(target))
             for var in data.domain.attributes]
        return numpy.array(I, dtype=bool)
    elif isinstance(group, RowGroup):
        target = set(indices)
        X, _ = data.get_column_view(group.var)
        I = numpy.zeros_like(X, dtype=bool)
        for i in target:
            I |= X == i
        return I
    else:
        raise TypeError("ColumnGroup or RowGroup expected, got {}"
                        .format(type(group).__name__))
