from Orange.bio import obiMeSH

names = ['1,3-diallylurea','1,7-octadiene','1,8-nonadiene'
'2,6-Dmhe','2-Dimethylamnoethyl cloride','2-cyanoethyl acrylate','5-fluorouracil']

print "cid  name    smiles  mesh"
print "string   string  string  string"
print ""

d = obiMeSH.obiMeSH()
chem = obiMeSH.pubChemAPI()

for i in names:
    cids = chem.getCIDs(i.strip())
    #print i, cids
    if len(cids) > 0:
        names.remove(i)
        cid = cids[0] # we hope that first CID is the right one
        # MeSH terms
        print cid
        print chem.getMeSHterms([int(cid)])
        terms = chem.getMeSHterms([int(cid)])[cid]
        smiles = chem.getSMILE(int(cid),"cid")
        print cid, "\t", i, "\t", smiles, "\t", terms
        
#print "Manualy annotate following chemicals -> ", names
