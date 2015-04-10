import orangecontrib.bio.taxonomy

for taxid in orangecontrib.bio.taxonomy.common_taxids():
    print("%-6s %s" % (taxid, orangecontrib.bio.taxonomy.name(taxid)))
