import orangecontrib.bio.geo

gds = orangecontrib.bio.geo.GDS("GDS1676")

data = gds.getdata(report_genes=True, transpose=False)
print("report_genes=True, transpose=False")
print("Report=Genes, Rows=Genes/Spots")
print("rows=%d cols=%d has_class=%s" % \
    (len(data), len(data.domain.attributes), data.domain.class_var!=None))
print("")

data = gds.getdata(report_genes=False, transpose=False)
print("report_genes=False, transpose=False")
print("Report=Spots, Rows=Genes/Spots")
print("rows=%d cols=%d has_class=%s" % \
    (len(data), len(data.domain.attributes), data.domain.class_var!=None))
print("")

data = gds.getdata(report_genes=True, transpose=True)
print("report_genes=True, transpose=True")
print("Report=Genes, Rows=Samples")
print("rows=%d cols=%d has_class=%s" % \
    (len(data), len(data.domain.attributes), data.domain.class_var!=None))
print("Class values:" + " ".join([str(cv) for cv in data.domain.class_var.values]) )
print("")

data = gds.getdata(report_genes=True, transpose=True, sample_type="infection")
print('report_genes=True, transpose=True sample_type="infection"')
print("Report=Genes, Rows=Samples")
print("rows=%d cols=%d has_class=%s" % \
    (len(data), len(data.domain.attributes), data.domain.class_var!=None))
print("Class values: " + " ".join([str(cv) for cv in data.domain.class_var.values]) )
print("")
