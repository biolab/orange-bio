##interval:7
from common import *

import sys, os
from gzip import GzipFile
import tempfile
from Orange.bio.obiGeneSets  import upload_genesets

upload_genesets(sf_server)


