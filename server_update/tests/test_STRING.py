from unittest import TestCase
from orangecontrib.bio import ppi


class StringTest(TestCase):

    def setUp(self):
        exclude = ['272634', '5476']
        self.taxids = [idtax for idtax in ppi.STRING.common_taxids() if idtax not in exclude]

    def test_size(self):
        for tid in self.taxids:
            self.string = ppi.STRING(taxid=tid)
            self.assertGreater(len(self.string.ids()), 0)
            self.assertEqual(len(self.string.organisms()), 1)

            self.detailed = ppi.STRINGDetailed(taxid=tid)
            self.assertGreater(len(self.detailed.ids()), 0)
            self.assertEqual(len(self.detailed.organisms()), 1)
