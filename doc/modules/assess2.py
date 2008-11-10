import orange
import obiAssess
import obiGeneSets

gs = obiGeneSets.collections([":kegg:hsa"])
data = orange.ExampleTable("DLBCL.tab")

asl = obiAssess.AssessLearner()
ass = asl(data, "hsa", geneSets=gs)

def genesetsAsAttributes(data, ass, domain=None):
    """
    Construct new data set with gene sets as attributes from data
    set "data" with assess model "ass".
    """

    ares = {}
    for ex in data:
        cres = ass(ex)
        for name,val in cres.items():
            aresl = ares.get(name, [])
            aresl.append(val)
            ares[name] = aresl

    ares = sorted(ares.items())

    if not domain: #construct new domain instance if needed
        domain = orange.Domain([ orange.FloatVariable(name=name) \
            for name in [ a[0] for a in ares]], data.domain.classVar )

    examples = [ [ b[zap] for a,b in ares ] + \
        [ data[zap][-1] ]   for zap in range(len(data)) ]

    et = orange.ExampleTable(domain, examples)
    return et

tdata = genesetsAsAttributes(data, ass)

print "First 10 attributes of the first example in transformed data set"
for pathw, enric in zip(tdata.domain,tdata[0])[:10]:
    print pathw.name, enric.value
