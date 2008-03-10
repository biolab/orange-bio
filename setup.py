from distutils.core import setup, Extension

setup(name = "Genomics",
      version = "VERSION",
      description = "Genomics extensions for Orange",
      author="University of Ljubljana, AI lab",
      author_email="tomaz.curk@fri.uni-lj.si",
      packages = [ 'widgets', 'data', 'doc' ],
      extra_path="Genomics",
      py_modules = [ 'orngKEGG' ],
      scripts=["registerWidgets.py"]
      )
