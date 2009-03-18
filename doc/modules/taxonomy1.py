import obiTaxonomy

for taxid in obiTaxonomy.common_taxids():
    print "%-5s %s" % (taxid, obiTaxonomy.name(taxid))
