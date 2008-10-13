import distutils
import distutils.sysconfig
import sys
import os
import orngEnviron

try:
    import orngRegistry
    orngRegistry.addWidgetCategory("Orange Genomics", \
        os.path.join(distutils.sysconfig.get_python_lib(),"Genomics","widgets"), \
        "remove" not in sys.argv[1] if len(sys.argv)>1 else 1)
except Exception, detail:
    print "Error: ", Exception, detail
    print "Add-on could not be registered with Orange Canvas. Please, make sure to install Orange and Orange Canvas before installing this add-on."
    print

try:
    from OWWidget import *
    import OWDatabasesUpdate
    app = QApplication(sys.argv)
    w = OWDatabasesUpdate.OWDatabasesUpdate(wantCloseButton=True, showAll=True, searchString="essential")
    w.show()
    w.exec_()
except Exception, detail:
    print "Error: ", Exception, detail
    print "Widget OWDatabasesUpdate failed to run. Databases needed by Orange Genomics add-on can be downloaded and updated only in this widget."
    # need to add support for update from text terminal
