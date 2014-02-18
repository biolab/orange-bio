import Orange

gds = Orange.bio.geo.GDS("GDS1676")
data = gds.getdata(sample_type="infection")
print "Genes: %d, Samples: %d" % (len(data), len(data.domain.attributes))

for a in data.domain.attributes:
    print a.name, a.attributes
