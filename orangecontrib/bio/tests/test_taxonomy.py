
import unittest
import errno

from orangecontrib.bio import taxonomy


def islocal():
    try:
        taxonomy.serverfiles.info(
            taxonomy.Taxonomy.DOMAIN, taxonomy.Taxonomy.FILENAME)
    except OSError as e:
        if e.errno == errno.EEXIST:
            return False
        else:
            raise
    else:
        return True


@unittest.skipIf(not islocal(), "Taxonomy not available")
class TestTaxonomy(unittest.TestCase):

    def test_name(self):
        tax = taxonomy.Taxonomy()
        self.assertEqual(tax["9606"], "Homo sapiens")
        self.assertEqual(tax["4932"], "Saccharomyces cerevisiae")

    def test_other_names(self):
        tax = taxonomy.Taxonomy()
        names = tax.other_names("4932")

        self.assertIn(("brewer's yeast", "common name"), names)
        self.assertIn(("lager beer yeast", "common name"), names)

    def test_search(self):
        tax = taxonomy.Taxonomy()
        res = tax.search("yeast",)
        self.assertIn("4932", res)

        res = tax.search("human", exact=True)
        self.assertEqual(res, ["9606"])

    def test_lineage(self):
        tax = taxonomy.Taxonomy()
        lineage = tax._tax.lineage("9606")
        self.assertEqual(lineage[0], "1")
        self.assertEqual(lineage[-1], "9605")
