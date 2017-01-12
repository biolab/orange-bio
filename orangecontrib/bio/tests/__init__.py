import os
import unittest


def suite(loader=None, pattern='test_*.py'):
    """
    Return the default test suite.
    """
    test_dir = os.path.dirname(__file__)
    kegg_dir = os.path.normpath(os.path.join(test_dir, "..", "kegg"))
    if loader is None:
        loader = unittest.TestLoader()
    if pattern is None:
        pattern = 'test_*.py'

    thisdir = os.path.dirname(__file__)

    top_level_dir = os.path.join(thisdir, "..", "..", "..")
    top_level_dir = os.path.realpath(top_level_dir)
    all_tests = [
        loader.discover(test_dir, pattern, top_level_dir),
        loader.discover(kegg_dir, pattern, top_level_dir),
    ]

    return unittest.TestSuite(all_tests)


def load_tests(loader, tests, pattern):
    return suite(loader, pattern)


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
