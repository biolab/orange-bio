import contextlib
import itertools
from itertools import chain

from typing import Iterable, Generator, Tuple

from AnyQt.QtCore import Qt

__all__ = ["group_ranges", "disconnected"]


@contextlib.contextmanager
def disconnected(boundsig, slot, conntype=Qt.UniqueConnection):
    """
    Disconnect a slot from a bound signal

    Parameters
    ----------
    boundsig : pyqtBoundSignal
    slot : pyqtBoundSlot
    conntype : Qt.ConnectionType

    """
    boundsig.disconnect(slot)
    try:
        yield
    finally:
        boundsig.connect(slot, conntype)


def group_ranges(indices):
    # type: (Iterable[int]) -> Iterable[Tuple[int, int]]
    """
    Group consecutive indices into `(start, stop)` tuple ranges.

    Parameters
    ----------
    indices : Iterable[int]
        An iterable yielding indices

    Example
    -------
    >>> list(group_ranges([1, 2, 3, 5, 3, 4]))
    [(1, 4), (5, 6), (3, 5)]
    >>> seq = [1, 2, 3, 3, 3, 5, 6, 9, 8]
    >>> assert seq == list(chain.from_iterable(
    ...        (range(s, e) for s, e in group_ranges(seq))))
    ...
    """
    g = itertools.groupby(enumerate(indices),
                          key=lambda t: t[1] - t[0])
    for _, range_ind in g:
        range_ind = list(range_ind)
        _, start = range_ind[0]
        _, end = range_ind[-1]
        yield start, end + 1
