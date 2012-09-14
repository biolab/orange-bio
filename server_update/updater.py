import subprocess
import time, glob

from datetime import datetime

from common import *

if not args:
    args = ["updateTaxonomy.py", "updateGO.py", "updateMeSH.py", "updateNCBI_geneinfo.py",
            "updateHomoloGene.py", "updateDictyBase.py", "updatePPI.py"]
    
for script in args:
    log = open(os.path.join(option.log_dir, script + ".log.txt"), "wb")
    p = subprocess.Popen([sys.executable, script, "-u", option.user, "-p", option.password], stdout=log, stderr=log)
    while p.poll() is None:
        time.sleep(3)
    log.write("\n" + script + " exited with exit status %s" % p.poll())
    log.close()
    if option.mailto:
        fromaddr = "orange@fri.uni-lj.si"
        toaddr = option.mailto.split(",")
        msg = open(os.path.join(option.log_dir, script + ".log.txt"), "rb").read()
        msg = "From: %s\r\nTo: %s\r\nSubject: Error running %s update script\r\n\r\n" % (fromaddr, ",".join(toaddr), script) + msg
        try:
            import smtplib
            s = smtplib.SMTP('212.235.188.18', 25)
            s.sendmail(fromaddr, toaddr, msg)
            s.quit()
        except Exception, ex:
            print "Failed to send error report due to:", ex
    

def files_report():
    sf = serverfiles.ServerFiles()
    html = []
    for domain in sf.listdomains():
        if domain not in ["demo", "demo2", "test", "gad"]:
            allinfo = sf.allinfo(domain)
            html += ["<h2>%s</h2>" % domain,
                     "<table><tr><th>Title</th><th>Date</th><th>Filename</th></tr>"] + \
                    ["<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (info["title"], info["datetime"], file) \
                     for file, info in allinfo.items()] + \
                    ["</table>"]
    return "\n".join(html)
  
open(os.path.join(option.log_dir, "serverFiles.html"), "wb").write(files_report())
