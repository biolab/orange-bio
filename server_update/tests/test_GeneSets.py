from unittest import TestCase
from server_update import sf_local
from orangecontrib.bio.geneset import filename_parse, load_serverfiles


class GeneSetsTest(TestCase):

    def setUp(self):
        self.files = [fn for domain, fn in sf_local.listfiles('gene_sets') if fn != 'index.pck']

    def test_gene_sets(self):
        for file in self.files:
            hierarchy, organism = filename_parse(file)
            load_serverfiles(hierarchy, organism)
