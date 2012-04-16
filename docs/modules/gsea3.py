import obiGeneSets
import obiGsea
import orange
import obiGene
import obiGEO

import obiGEO
gds = obiGEO.GDS("GDS10")
data = gds.getdata() 

print "Possible phenotype descriptors:"
print map(lambda x: x[0], obiGsea.allgroups(data).items())

matcher=obiGene.matcher([obiGene.GMKEGG("9606")])

phenVar = "tissue"
geneVar = "gene" #use gene meta variable for gene names

genesets =  obiGeneSets.collections([":kegg:hsa"])
res = obiGsea.runGSEA(data, matcher=matcher, minPart=0.05, geneSets=genesets, 
    permutation="class", n=10, phenVar=phenVar, geneVar=geneVar)

print
print "GSEA results (choosen descriptor: tissue)"
print "%-40s %6s %6s %6s %7s" % ("LABEL", "NES", "FDR", "SIZE", "MATCHED") 
for name,resu in sorted(res.items(), key=lambda x: x[1]["fdr"])[:10]: 
    print "%-40s %6.3f %6.3f %6d %7d" % (name[:30], resu["nes"], resu["fdr"], 
        resu["size"], resu["matched_size"]) 
