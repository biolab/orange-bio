from unittest import TestCase
from server_update import sf_local
from orangecontrib.bio import taxonomy
from orangecontrib.bio.widgets3.OWGeneInfo import ncbi_info
from orangecontrib.bio.gene import NCBIGeneInfo


class NCBIGeneInfoTest(TestCase):

    def setUp(self):
        self.listfiles = [fname for domain, fname in sf_local.listfiles("NCBI_geneinfo")]
        self.taxonomy_IDs = sorted(set([name.split(".")[-2] for name in self.listfiles] + NCBIGeneInfo.common_taxids()))
        self.organisms = [(taxonomy.name(tax_id), tax_id) for tax_id in self.taxonomy_IDs]

    def test_GDS5245(self):
        organism, tax_id = ('Homo sapiens', '9606')  # (Organism, Taxon identifier)
        test_genes = ['ADAP1', 'ADAP2']

        n_info = ncbi_info(tax_id, test_genes)
        info_rows = n_info[1]
        self.assertIsInstance(info_rows, list)
        for row in info_rows:
            self.assertIsNotNone(row)

    def test_genes_count(self):
        for organism, tax_id in self.organisms:
            self.assertGreater(len(NCBIGeneInfo(organism)), 0)
