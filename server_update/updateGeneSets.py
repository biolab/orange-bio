##interval:7
from common import *

import sys, os
from gzip import GzipFile
import tempfile
from orangecontrib.bio.geneset  import upload_genesets

upload_genesets(sf_server)
