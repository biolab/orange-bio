PYTHONINCLUDE = $(shell python -c "import sys; print sys.prefix+'/include/python'+sys.version[:3]")
PYTHONLIB = $(shell python -c "import sys; print sys.prefix+'/lib/python'+sys.version[:3]")
PYTHONSITEPKGS = $(PYTHONLIB)/site-packages
OS = $(shell uname)

COMPILER = g++
CCOMPILER = gcc
LINKER = g++

ifeq ($(OS), Darwin)
  COMPILEOPTIONS = -fPIC -fpermissive -fno-common -w -O3 -I$(PYTHONINCLUDE)
  LINKOPTIONS = -dynamiclib -headerpad_max_install_names -undefined dynamic_lookup -lstdc++ $(LDFLAGS)
else
  COMPILEOPTIONS = -fPIC -fpermissive -w -O3 -I$(PYTHONINCLUDE)
  LINKOPTIONS = -shared -lstdc++ $(LDFLAGS)
endif

all: modules

Makefile.depends: ../setup.py
	python -c "import sys; sys.path=['..']+sys.path; import setup; setup.writeMakeFileDepends();"

include Makefile.depends

%.o: %.cxx
	$(COMPILER) $(COMPILEOPTIONS) $(CXXFLAGS) $(COMPILEOPTIONSMODULES) -c $< -o $@

%.cxx: %.i
	swig -python $<

%.o: %.cpp
	$(COMPILER) $(COMPILEOPTIONS) $(CXXFLAGS) $(COMPILEOPTIONSMODULES) -c $< -o $@

%.o: %.c
	$(CCOMPILER) $(COMPILEOPTIONS) $(CFLAGS) $(COMPILEOPTIONSMODULES) -c $< -o $@

clean:
	rm -f *.o ../*.so
	rm Makefile.depends

.PHONY: all clean
