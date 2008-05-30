import distutils
import distutils.sysconfig
import sys
import orngRegistry
import os

orngRegistry.addWidgetCategory("Orange Genomics", \
    os.path.join(distutils.sysconfig.get_python_lib(),"Genomics","widgets"), \
    "remove" not in sys.argv[1] if len(sys.argv)>1 else 1)
