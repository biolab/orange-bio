from unittest import TestCase
from orangecontrib.bio.gene import homology


class HomologyTest(TestCase):

    def test_homolo(self):
        homology.HomoloGene()

    def test_InParanoid(self):
        homology.InParanoid()
