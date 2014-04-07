import Orange.bio.taxonomy

for taxid in Orange.bio.taxonomy.common_taxids():
    print "%-6s %s" % (taxid, Orange.bio.taxonomy.name(taxid))
