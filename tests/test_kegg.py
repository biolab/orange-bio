import unittest
import obiKEGG2 as kegg

#from obiKEGGservices import 
#from 

#class TestServices(unittest.TestCase):
#    def test_services(self):
#        service = kegg.web_service()
#        service.binfo("gb")
#        service.bfind("gb E-cadherin human")
#        service.bget("eco:b0002 hin:tRNA-Cys-1")
        
        
class TestGenome(unittest.TestCase):
    def test_genome(self):
        genome = kegg.KEGGGenome()
        org_codes = genome.keys()
        
        for code in org_codes[:10]:
            self.assertTrue(code in genome)
            self.assertTrue(genome.has_key(code))
            entry = genome[code]
            self.assertEqual(entry.entry_key, code)
            self.assertIsInstance(entry, genome.ENTRY_TYPE)
            self.assertIsInstance(entry.name, str)
            self.assertIsInstance(entry.taxid, str)
            
            
        self.assertTrue(genome.get_entry("No such name --,,;;';p[[]&&&") is None)
        self.assertRaises(KeyError, genome.__getitem__, ("No such name --,,;;';p[[]&&&"))
        
        self.assertTrue(genome.search("homo sapiens")[0] == "hsa")
        entry = genome['hsa']
        self.assertEqual(entry.taxid, "9606")
        
        
class TestGenes(unittest.TestCase):
    def _tester(self, org):
        genes = kegg.KEGGGenes(org)
        keys = genes.keys()
        for gene in keys[:10]:
            self.assertTrue(gene in genes)
            self.assertTrue(genes.has_key(gene))
            entry = genes[gene]
            self.assertTrue(gene.endswith(entry.entry_key))
            self.assertIsInstance(entry, genes.ENTRY_TYPE)
            self.assertIsInstance(entry.aliases(), list)
            self.assertTrue(all(isinstance(a, basestring) for a in entry.aliases()))
            
    def test_hsa(self):
        self._tester("hsa")
        
    def test_sce(self):
        self._tester("sce") 
        
    
class TestPathways(unittest.TestCase):
    def _tester(self, id):
        pathways = kegg.KEGGPathways()
        
        self.assertRaises(KeyError, pathways.__getitem__, ("--invalid--"))
        pathway = pathways[id]
        self.assertTrue(id.endswith(pathway.entry_key))
        self.assertIsInstance(pathway, pathways.ENTRY_TYPE)
        
        path = kegg.KEGGPathway(id)
        
    def test_1(self):
        self._tester("ec00190")
            
    def test_2(self):
        self._tester("hsa00190")
        
        
class TestOrganism(unittest.TestCase):
    def _tester(self, org):
        self.organism = org = kegg.KEGGOrganism(org)
        genes = org.genes
        self.assertTrue(all(g.startswith(org.org_code) for g in genes))
        pathways = org.pathways()
        pathways_for_genes = org.pathways(with_ids=list(genes)[:5])
        self.assertTrue(all(p.startswith("path:" + org.org_code) \
                            for p in pathways + pathways_for_genes))
        
    def test_hsa(self):
        self._tester("hsa")
        # Test search
        self.assertEqual(self.organism.org_code,
                         kegg.KEGGOrganism("Homo sapiens").org_code)
        
    def test_sce(self):
        self._tester("sce")
        
class TestOld(unittest.TestCase):
    def test(self):
        p = kegg.KEGGPathway("sce00010")
        print p.genes
        print p.reactions
        print p.compounds
        print p.image
        g = kegg.KEGGGenome()
        org = kegg.KEGGOrganism("Homo sapiens")
        print list(org.genes)[:10]
#        org.gene_aliases
        print org.pathways(with_ids=org.genes.keys()[:5])
#        print org.enzymes()
        print org.get_enriched_pathways(org.genes.keys()[:10])
        print org.genematcher