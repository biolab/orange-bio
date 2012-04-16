import orange
import obiAssess
import obiGeneSets

gs = obiGeneSets.collections([":kegg:hsa"])
data = orange.ExampleTable("DLBCL.tab")

asl = obiAssess.AssessLearner()

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

offer = None

def transformLearningS(data):
    ass = asl(data, "hsa", geneSets=gs)
    et = genesetsAsAttributes(data, ass)

    global offer
    offer = (et.domain, ass) #save assess model

    return et
   
def transformTestingS(data):
    global offer
    if not offer:
        a = fdfsdsdd #exception

    domain, ass = offer
    offer = None

    return genesetsAsAttributes(data, ass, domain)


import orngBayes, orngTest, orngStat
learners = [ orngBayes.BayesLearner() ]

resultsOriginal = orngTest.crossValidation(learners, data, folds=10)
resultsTransformed = orngTest.crossValidation(learners, data, folds=10, 
    pps = [("L", transformLearningS), ("T", transformTestingS)])

print "Original", "CA:", orngStat.CA(resultsOriginal), "AUC:", orngStat.AUC(resultsOriginal)
print "Transformed", "CA:", orngStat.CA(resultsTransformed), "AUC:", orngStat.AUC(resultsTransformed)

