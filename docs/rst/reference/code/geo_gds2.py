from Orange.bio import obiGEO

gds = obiGEO.GDS("GDS1210")

data = gds.getdata(report_genes=True, transpose=False)
print "report_genes=True, transpose=False"
print "Report=Genes, Rows=Genes/Spots"
print "rows=%d cols=%d has_class=%s" % (len(data), len(data.domain.attributes), data.domain.classVar<>None)
print

data = gds.getdata(report_genes=False, transpose=False)
print "report_genes=False, transpose=False"
print "Report=Spots, Rows=Genes/Spots"
print "rows=%d cols=%d has_class=%s" % (len(data), len(data.domain.attributes), data.domain.classVar<>None)
print

data = gds.getdata(report_genes=True, transpose=True)
print "report_genes=True, transpose=True"
print "Report=Genes, Rows=Samples"
print "rows=%d cols=%d has_class=%s" % (len(data), len(data.domain.attributes), data.domain.classVar<>None)
print "Class values:", " ".join([str(cv) for cv in data.domain.classVar.values]) 
print


data = gds.getdata(report_genes=True, transpose=True, sample_type="tissue")
print 'report_genes=True, transpose=True sample_type="tissue"'
print "Report=Genes, Rows=Samples"
print "rows=%d cols=%d has_class=%s" % (len(data), len(data.domain.attributes), data.domain.classVar<>None)
print "Class values:", " ".join([str(cv) for cv in data.domain.classVar.values]) 
print
