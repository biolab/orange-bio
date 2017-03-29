from unittest import TestCase
from orangecontrib.bio.ppi import BioGRID


class BioGRIDTest(TestCase):

    def setUp(self):
        self.biogrid = BioGRID()

    def test_size(self):
        """ test if initialization of database was successful
        """
        self.assertGreater(len(self.biogrid.ids()), 60000)  # all protein ids (biogrid_id_interactors)
        self.assertGreater(len(self.biogrid.all_edges()), 1300000)  # all edges
        self.assertGreater(len(self.biogrid.all_edges_annotated()), 1300000)  # all edges annotated
