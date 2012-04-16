import distutils
import distutils.sysconfig
import sys
import os

if len(sys.argv) <= 1 or len(sys.argv) > 1 and "remove" not in sys.argv[1]:
    try:
        print "Downloading essential genomics databases"
        if sys.platform == "win32":
            os.system("start " + sys.prefix+"""\\python.exe -c "from Orange.orng import orngServerFiles; orngServerFiles.update_by_tags(tags=['essential'])" """)
        else:
            from Orange.orng import orngServerFiles
            orngServerFiles.update_by_tags(tags=["essential"])
    except Exception, ex:
        print ex
        pass

    print "To download and/or update the databases needed by Orange"
    print "Bioinformatics please use the \"Update Genomics Databases\""
    print "widget in Orange Canvas or use the orngServerFiles module:"
    print
    print "    from Orange.orng import orngServerFiles"
    print "    orgnServerFiles.consoleupdate()"
    print
