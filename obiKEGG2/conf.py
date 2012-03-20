"""
obiKEGG2 configuration 

mostly just caching settings

"""

import os
import ConfigParser
from StringIO import StringIO
from Orange.utils import serverfiles
kegg_dir = serverfiles.localpath("KEGG2")

default = """
[cache]
# path = %(home)s/.obiKEGG/
path = %(kegg_dir)s/
store = sqlite3
kegg_invalidate = always

[service]
transport = urllib2
# transport = requests

"""

# Orange kegg files dir
env = dict(os.environ)
env["kegg_dir"] = kegg_dir

parser = ConfigParser.ConfigParser(env)


parser.readfp(StringIO(default), "default")

# TODO: global settings rc file
parser.read([os.path.expanduser("~/.obiKEGG/rc.cfg")])

params = {}


for section in parser.sections():
    for option in parser.options(section):
        params[section + "." + option] = parser.get(section, option)
