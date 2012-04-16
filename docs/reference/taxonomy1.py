import obiTaxonomy

for taxid in obiTaxonomy.common_taxids():
    print "%-6s %s" % (taxid, obiTaxonomy.name(taxid))
