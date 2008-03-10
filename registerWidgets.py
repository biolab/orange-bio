import distutils
import distutils.sysconfig
import sys
import orngRegistry

orngRegistry.addWidgetCategory("Orange Genomics", distutils.sysconfig.get_python_lib() + r"\Genomics\widgets", "remove" not in sys.argv[1])
