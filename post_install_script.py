import distutils
import distutils.sysconfig
import sys
import os

if len(sys.argv) <= 1 or len(sys.argv) > 1 and "remove" not in sys.argv[1]:
    print "To download and/or update the databases needed by Orange Bioinformatics"
    print "please use the \"Update Genomics Databases\" widget in Orange Canvas"
    print "or use the orngServerFiles module, e.g.:"
    print
    print "    import orngServerFiles"
    print "    orgnServerFiles.consoleupdate()"
    print
