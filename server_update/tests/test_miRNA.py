from unittest import TestCase
from orangecontrib.bio.obimiRNA import ids, get_info


class MiRNATest(TestCase):
    # is this correct?
    # "currently only Mus musculus is supported" -> obimiRNA.load_miRNA_microCosm

    def test_miRNAs_library(self):
        self.assertGreater(len(ids()), 700)
        # test if mat_miRNA class is constructed
        for identifier in ids():
            get_info(identifier)
