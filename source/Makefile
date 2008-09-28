PYTHONINCLUDE = $(shell python -c "import sys; print sys.prefix+'/include/python'+sys.version[:3]")
OS = $(shell uname)

COMPILER = g++
CCOMPILER = gcc
LINKER = g++

ifeq ($(OS), Darwin)
  COMPILEOPTIONS = -arch i386 -arch ppc -fPIC -fpermissive -fno-common -w -D_POSIX_C_SOURCE -O3 -I$(PYTHONINCLUDE)
  LINKOPTIONS = -arch i386 -arch ppc -dynamiclib -headerpad_max_install_names -undefined dynamic_lookup -lstdc++
else
  COMPILEOPTIONS = -fPIC -fpermissive -w -O3 -I$(PYTHONINCLUDE)
  LINKOPTIONS = -shared -lstdc++
endif

all: modules

Makefile.depends: ../setup.py
	python -c "import sys; sys.path=['..']+sys.path; import setup; setup.writeMakeFileDepends();"

include Makefile.depends

%.o: %.cxx
	$(COMPILER) $(COMPILEOPTIONS) $(COMPILEOPTIONSMODULES) -c $< -o $@

%.cxx: %.i
	swig -python $<

%.o: %.cpp
	$(COMPILER) $(COMPILEOPTIONS) $(COMPILEOPTIONSMODULES) -c $< -o $@

%.o: %.c
	$(CCOMPILER) $(COMPILEOPTIONS) $(COMPILEOPTIONSMODULES) -c $< -o $@

clean:
	rm -f *.o ../*.so
	rm Makefile.depends

.PHONY: all clean
