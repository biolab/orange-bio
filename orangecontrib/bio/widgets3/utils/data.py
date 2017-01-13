"""
Table and other data manipulation utilities
"""
import numbers
from operator import itemgetter
from typing import Sequence, Tuple

import numpy

import Orange.data

ColSpec = Sequence[Tuple[Orange.data.Variable, Sequence[numbers.Real]]]


def append_columns(data, attributes=(), class_vars=(), metas=()):
    # type: (Orange.data.Table, ColSpec, ColSpec, ColSpec) -> Orange.data.Table
    """
    Append a set of columns to a data table.


    Parameters
    ----------
    data : Orange.data.Table
        Primary table.
    attributes : Sequence[Tuple[Orange.data.Variable], Sequence[float]]
        A Sequence of variable and column data tuples to append to the
        `data`.
    class_vars : Sequence[Tuple[Orange.data.Variable], Sequence[float]]
        A Sequence of variable and column data tuples to append to the
        `data` in the
    metas : Sequence[Tuple[Orange.data.Variable], Sequence[float]]
        A Sequence of variable and column data tuples to append to the
        `data`

    Returns
    -------
    data : Orange.data.Table
        A copy of the original `data` input extended with all columns from
        `attributes`, `class_vars`, `metas` parameters

    Note
    ----
    All variables in the original and new columns should be distinct.
    """
    domain = data.domain
    new_attributes = tuple(map(itemgetter(0), attributes))
    new_class_vars = tuple(map(itemgetter(0), class_vars))
    new_metas = tuple(map(itemgetter(0), metas))

    new_domain = Orange.data.Domain(
        domain.attributes + new_attributes,
        domain.class_vars + new_class_vars,
        domain.metas + new_metas
    )

    def ascolumn(array, n):
        # type: (Sequence[float], int) -> numpy.ndarray
        array = numpy.asarray(array)
        if array.ndim < 2:
            array = array.reshape((n, 1))
        return array
    N = len(data)

    attr_cols = [ascolumn(col, N) for _, col in attributes]
    class_cols = [ascolumn(col, N) for _, col in class_vars]
    meta_cols = [ascolumn(col, N) for _, col in metas]

    new_data = data.from_table(new_domain, data)

    for i, (var, col) in enumerate(zip(new_attributes, attr_cols),
                                   start=len(domain.attributes)):
        assert new_data.domain.attributes[i] is var
        new_data.X[:, i] = col.ravel()

    for i, (var, col) in enumerate(zip(new_class_vars, class_cols),
                                   start=len(domain.class_vars)):
        assert new_data.domain.class_vars[i] is var
        new_data._Y[:, i] = col.ravel()

    for i, (var, col) in enumerate(zip(new_metas, meta_cols),
                                   start=len(domain.metas)):
        assert new_data.domain.metas[i] is var
        new_data.metas[:, i] = col.ravel()

    return new_data


import unittest


class TestDataUtils(unittest.TestCase):
    def test_append_columns(self):
        data = Orange.data.Table.from_numpy(
            None,
            numpy.array([[1, 2], [2, 3]]),
            numpy.array([[1], [2]]),
            # numpy.array([[], []])
        )
        dv = Orange.data.DiscreteVariable("D", values=("A", "B"))
        cv = Orange.data.ContinuousVariable("C")
        r = append_columns(data, attributes=[(dv, [0, 1])])
        self.assertEqual(r.domain.attributes, data.domain.attributes + (dv,))
        self.assertEqual(r.domain.class_vars, data.domain.class_vars)
        self.assertEqual(r.domain.metas, data.domain.metas)
        numpy.testing.assert_equal(r.X[:, -1], [0, 1])
        numpy.testing.assert_equal(r.X[:, :-1], data.X)
        numpy.testing.assert_equal(r.Y, data.Y)
        numpy.testing.assert_equal(r.metas, data.metas)

        r = append_columns(data, class_vars=[(cv, [.0, .5])])
        self.assertEqual(r.domain.attributes, data.domain.attributes)
        self.assertEqual(r.domain.class_vars, data.domain.class_vars + (cv,))
        self.assertEqual(r.domain.metas, data.domain.metas)
        numpy.testing.assert_equal(r.X, data.X)
        numpy.testing.assert_equal(r._Y[:, :-1], data._Y)
        numpy.testing.assert_equal(r._Y[:, -1], [.0, .5])
        numpy.testing.assert_equal(r.metas, data.metas)

        r = append_columns(data, metas=[(dv, [0, 1])])
        self.assertEqual(r.domain.attributes, data.domain.attributes)
        self.assertEqual(r.domain.class_vars, data.domain.class_vars)
        self.assertEqual(r.domain.metas, data.domain.metas + (dv,))
        numpy.testing.assert_equal(r.X, data.X)
        numpy.testing.assert_equal(r._Y, r._Y)
        numpy.testing.assert_equal(r.metas[:, :-1], data.metas)
        numpy.testing.assert_equal(r.metas[:, -1], [0, 1])

        data_0 = data[:0]
        r = append_columns(data_0, attributes=[(cv, [])], metas=[(dv, [])])
        self.assertEqual(r.domain.attributes, data.domain.attributes + (cv,))
        self.assertEqual(r.domain.metas, data.domain.metas + (dv,))
