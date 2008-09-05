import sys
import go
try:
    if sys.argv[1]=="remove":
        sys.exit(0)
except:
    pass
print "Downloading GO database"
go.downloadGO()
print "Downloading annotation for yeast genome"
go.downloadAnnotation()
