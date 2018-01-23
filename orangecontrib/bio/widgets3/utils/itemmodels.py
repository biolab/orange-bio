import numbers
from collections import namedtuple
from typing import Sequence

from AnyQt.QtCore import Qt, QSortFilterProxyModel
from AnyQt.QtWidgets import QStyledItemDelegate


class FilterProxyModel(QSortFilterProxyModel):
    """
    A simple filter proxy model with settable filter predicates

    Example
    -------
    >>> proxy = FilterProxyModel()
    >>> proxy.setFilters([
    ...     FilterProxyModel.Filter(0, Qt.DisplayRole, lambda value: value < 1)
    ... ])
    """
    Filter = namedtuple("Filter", ["column", "role", "predicate"])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__filters = []

    def setFilters(self, filters):
        # type: (Sequence[FilterProxyModel.Filter]) -> None
        filters = [FilterProxyModel.Filter(f.column, f.role, f.predicate)
                   for f in filters]
        self.__filters = filters
        self.invalidateFilter()

    def filterAcceptsRow(self, row, parent):
        source = self.sourceModel()

        def apply(f):
            index = source.index(row, f.column, parent)
            data = source.data(index, f.role)
            try:
                return f.predicate(data)
            except (TypeError, ValueError):
                return False

        return all(apply(f) for f in self.__filters)


class RealDelegate(QStyledItemDelegate):
    """
    An Item delegate for displaying numerical columns
    """
    def __init__(self, parent=None, precision=4, **kwargs):
        super().__init__(parent, **kwargs)
        self.precision = precision

    def displayText(self, value, locale):
        if isinstance(value, numbers.Integral):
            return locale.toString(int(value))
        elif isinstance(value, numbers.Real):
            return locale.toString(float(value), "g", self.precision)
        else:
            return super().displayText(value, locale)

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        align = index.data(Qt.TextAlignmentRole)
        data = index.data(Qt.DisplayRole)
        if align is None and isinstance(data, numbers.Real):
            # Right align if the model does not specify otherwise
            option.displayAlignment = Qt.AlignRight | Qt.AlignVCenter
