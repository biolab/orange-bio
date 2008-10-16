from distutils.core import setup, Extension
import os, glob

# list all documentation files that need to be included
docFiles = []
for (dirp, dirns, n) in os.walk('doc'):
	nr = [n1.replace('\\', '/') for n1 in n]
	dirn = dirp.replace('\\', '/')[4:]
	if len(dirn):
		dirn = dirn + '/'
	docFiles.extend( [dirn + n1r for n1r in nr if '.svn' not in dirp + '/' + n1r] )

modules=[Extension("_GOLib", sources=["source/go.c"])]

def writeMakeFileDepends():
	f = open("Makefile.depends", "wt")

	includePaths = []
	for ext in modules:
		if ext.include_dirs <> []:
			for p in ext.include_dirs:
				includePaths.append( "-I%s" % (p))
	includePaths = " ".join(includePaths)

	if includePaths <> "":
		f.write("COMPILEOPTIONSMODULES = %s\n" % (includePaths))

	f.write("modules:")
	for ext in modules:
		f.write(" %s.so" % (ext.name))
	f.write("\n")

	for ext in modules:
		objs = []
		for s in ext.sources:
			if s[-2:] == '.c' or s[-4:] == '.cpp':
				objfname = os.path.splitext(os.path.join("..", s))[0] + ".o"
				objs.append( objfname)
		objs = " ".join(objs)
		f.write("%s.so: %s\n" % (ext.name, objs))
		f.write("\t$(LINKER) $(LINKOPTIONS) %s -o %s.so" % (objs, os.path.join("..", ext.name)))
	f.close()

if __name__ == "__main__":
	setup(name = "Orange Bioinformatics",
	      version = "1.0",
	      description = "Bioinformatics Add-On for Orange",
	      author="University of Ljubljana, AI lab",
	      maintainer_email="tomaz.curk@fri.uni-lj.si",
	      ext_modules=modules,
	      packages = [ 'widgets', 'doc' ],
	      package_data = {'widgets': ['icons/*.png'], 'doc': docFiles},
	      extra_path="Genomics",
	      py_modules = [ 'go', 'obiKEGG', 'obiGsea', 'obiGeneMatch', 'obiData', 'obiGenomicsUpdate', 'stats', 'pstat', 'obiExpression', 'obiGO', 'obiProb', 'obiAssess', 'obiGeneSets', 'obiMeSH', 'obiDicty' ],
	      scripts=["post_install_script.py"]
	      )
