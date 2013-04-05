import unittest
import doctest

from .. import api
from ... import obiKEGG


class TestApi(unittest.TestCase):
    pass


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(obiKEGG,
                                        optionflags=doctest.ELLIPSIS))
    return tests
