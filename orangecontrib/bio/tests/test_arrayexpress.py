import unittest
import doctest

from .. import arrayexpress


def load_tests(loader, tests, ignore):
    conn = arrayexpress.ArrayExpressConnection()
    foo = type("foo", (object,), {})()
    bar = type("bar", (object,), {})()

    tests.addTests(
        doctest.DocTestSuite(
            arrayexpress,
            extraglobs={"conn": conn, "foo": foo, "bar": bar},
            optionflags=doctest.ELLIPSIS)
    )
    return tests
