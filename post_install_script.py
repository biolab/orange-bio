import distutils
import distutils.sysconfig
import sys
import os

if len(sys.argv) <= 1 or len(sys.argv) > 1 and "remove" not in sys.argv[1]:
    try:
        ## on windows this is automaticaly run from the distutils installer
        if sys.platform == "win32":
            os.system("start " + sys.prefix+"""\\python.exe -c "import orngServerFiles; orngServerFiles.consoleupdate()" """)
        else:
            import orngServerFiles
            orngServerFiles.consoleupdate()

        print "Database donwload was successful."
        print
        print "For further updates of the databases needed by Orange"
        print "Bioinformatics please use the \"Update Genomics Databases\""
        print "widget in Orange Canvas or use the orngServerFiles module, e.g.:"
        print
        print "    import orngServerFiles"
        print "    orgnServerFiles.consoleupdate()"
        print        
    except Exception, detail:
        print "Database donwload did not succeed! Please, retry later."
        print

