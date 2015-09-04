import orangecontrib.bio.utils.serverfiles
import glob
import re

filenames = glob.glob(orangecontrib.bio.utils.serverfiles.localpath("GEO") + "/GDS*.soft.gz")
m = re.compile("(GDS[0-9]*).soft")
print("%d data files cached:" % len(filenames))
print(" ".join([m.search(fn).group(1) for fn in filenames]))

