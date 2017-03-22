from unittest import TestCase
from orangecontrib.bio import omim


class OMIMTest(TestCase):

    def test_size(self):
        self.assertGreater(len(omim.diseases()), 5600)
        self.assertGreater(len(omim.gene_diseases()), 8000)
