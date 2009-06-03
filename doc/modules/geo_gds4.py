import orngServerFiles
import glob
print "Cached data files:", len(glob.glob(orngServerFiles.localpath("GEO") + "/GDS*"))
