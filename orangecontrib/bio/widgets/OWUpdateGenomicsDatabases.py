"""
<name>Databases</name>
<description>Update of systems biology data and knowledge bases.</description>
<contact>Ales Erjavec</contact>
<priority>10</priority>
<icon>icons/Databases.svg</icon>
"""

from Orange.OrangeWidgets.OWDatabasesUpdate import *

NAME = "Databases"
DESCRIPTION = "Updates local system biology databases."
ICON = "icons/Databases.svg"
PRIORITY = 10

INPUTS = []
OUTPUTS = []

REPLACES = ["_bioinformatics.widgets.OWUpdateGenomicsDatabases.OWUpdateGenomicsDatabases"]


class OWUpdateGenomicsDatabases(OWDatabasesUpdate): 
    def __init__(self, parent=None, signalManager=None, name="Databases", **kwds):
        OWDatabasesUpdate.__init__(self, parent, signalManager, name, domains = \
                ["GO", "MeSH", "Taxonomy", "NCBI_geneinfo", "GEO", 
                 "dictybase", "OMIM", "HomoloGene", "Affy", "miRNA", "gene_sets",
                 "PPI"], **kwds)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = OWUpdateGenomicsDatabases()
##    app.setMainWidget(w)
    w.show()
    app.exec_()
    w.saveSettings()
