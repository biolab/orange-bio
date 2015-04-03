from io import StringIO
import doctest

import unittest

from orangecontrib.bio.kegg.entry import parser, fields, DBEntry, entry_decorate


TEST_ENTRY = """\
ENTRY       test_id    something else
NAME        test
DESCRIPTION This is a test's description.
            it spans
            multiple lines
  SUB       This is a description's sub
            section
///
"""


@entry_decorate
class Entry(DBEntry):
    pass


class TestEntry(unittest.TestCase):
    def test_entry(self):
        """
        Test basic DBEntry class.
        """
        entry = Entry(TEST_ENTRY)
        self.assertEqual(entry.entry_key, "test_id")
        self.assertEqual(entry.ENTRY.TITLE, "ENTRY")

        self.assertEqual(str(entry), TEST_ENTRY[:-4])


class TestParser(unittest.TestCase):
    def test_parser(self):
        parse = parser.DBGETEntryParser()
        stream = StringIO(TEST_ENTRY)

        for event, title, text in parse.parse(stream):
            pass


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(parser))
    return tests
