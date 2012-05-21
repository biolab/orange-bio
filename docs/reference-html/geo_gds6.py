import obiGEO
import orange
import orngTest
import orngStat

gds = obiGEO.GDS("GDS2960")
data = gds.getdata(sample_type="disease state", transpose=True)
print "Samples: %d, Genes: %d" % (len(data), len(data.domain.attributes))

learners = [orange.LinearLearner]
results = orngTest.crossValidation(learners, data, folds=10)
print "AUC = %.3f" % orngStat.AUC(results)[0]
