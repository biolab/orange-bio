import distutils
import distutils.sysconfig
import sys
import os

if len(sys.argv) <= 1 or len(sys.argv) > 1 and "remove" not in sys.argv[1]:
    try:
        import orngEnviron
        from OWWidget import *
        import OWUpdateGenomicsDatabases
        app = QApplication(sys.argv)
        w = OWUpdateGenomicsDatabases.OWUpdateGenomicsDatabases(wantCloseButton=True, showAll=True, searchString="essential")
        w.show()
        w.exec_()
    except Exception, detail:
        print "Error: ", Exception, detail
        print "Widget OWUpdateGenomicsDatabases failed to run. Databases needed by Orange Genomics add-on can be downloaded and updated only in this widget."
        # need to add support for update from text terminal
