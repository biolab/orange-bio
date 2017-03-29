from unittest import TestCase
from orangecontrib.bio.ppi import MIPS


class MIPSTest(TestCase):

    def setUp(self):
        self.mips = MIPS()

    def test_size(self):
        self.assertGreater(len([interaction for interaction in self.mips]), 1800)
