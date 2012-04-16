"""
<name>Databases</name>
<description>Update of systems biology data and knowledge bases.</description>
<contact>Ales Erjavec</contact>
<priority>10</priority>
<icon>icons/UpdateDatabases.png</icon>
"""

from Orange.OrangeWidgets.OWDatabasesUpdate import *

class OWUpdateGenomicsDatabases(OWDatabasesUpdate): 
    def __init__(self, parent=None, signalManager=None, name="Databases", **kwds):
        OWDatabasesUpdate.__init__(self, parent, signalManager, name, domains = \
                ["GO", "KEGG", "MeSH", "Taxonomy", "NCBI_geneinfo", "GEO", 
                 "dictybase", "OMIM", "HomoloGene", "Affy", "miRNA", "gene_sets",
                 "PPI"], **kwds)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWUpdateGenomicsDatabases()
##    app.setMainWidget(w)
    w.show()
    app.exec_()
    w.saveSettings()
