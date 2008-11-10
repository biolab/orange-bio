import orange
import obiAssess
import obiGeneSets

gs = obiGeneSets.collections([":kegg:hsa"])
data = orange.ExampleTable("DLBCL.tab")

asl = obiAssess.AssessLearner()
ass = asl(data, "hsa", geneSets=gs)

print "Enrichments for the first example (10 pathways)"
enrichments = ass(data[0])
for patw, enric in sorted(enrichments.items())[:10]:
    print patw, enric

