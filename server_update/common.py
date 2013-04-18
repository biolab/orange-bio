import Orange
import Orange.utils.serverfiles as serverfiles
import Orange.utils.environ as environ
import os, sys
import optparse
import gzip, shutil

usage="""usage: %prog [options] [update_script ...]"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-u", "--user", help="User name")
parser.add_option("-p", "--password", help="Password")
parser.add_option("-l", "--log-dir", dest="log_dir", help="Directory to store the logs", default="./")
parser.add_option("-m", "--mailto", help="e-mail the results to EMAIL", metavar="EMAIL", default=None)

option, args = parser.parse_args()

if not option.user or not option.password:
    print "Pass -u username -p password!"
    sys.exit(1)

sf_server = serverfiles.ServerFiles(option.user, option.password)
sf_local = Orange.utils.serverfiles
