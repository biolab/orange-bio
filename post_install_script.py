import distutils
import distutils.sysconfig
import sys
import os

if len(sys.argv) <= 1 or len(sys.argv) > 1 and "remove" not in sys.argv[1]:
    try:
        import orngEnviron, bla
        from OWWidget import *
        import OWUpdateGenomicsDatabases
        app = QApplication(sys.argv)
        w = OWUpdateGenomicsDatabases.OWUpdateGenomicsDatabases(wantCloseButton=True, showAll=True, searchString="essential")
        w.show()
        w.exec_()
    except Exception, detail:
        ## on windows this is automaticaly run from the distutils installer
        if sys.platform == "win32":
            os.system("start " + sys.prefix+"""\\python.exe -c "import orngServerFiles; orngServerFiles.consoleupdate()" """)
        else:
            import orngServerFiles
            orngServerFiles.consoleupdate()
            
        
