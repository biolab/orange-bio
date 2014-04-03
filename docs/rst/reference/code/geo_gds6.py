import Orange.bio.geo

gds = Orange.bio.geo.GDS("GDS2960")
data = gds.getdata(sample_type="disease state", transpose=True)
print "Samples: %d, Genes: %d" % (len(data), len(data.domain.attributes))

learners = [ Orange.classification.logreg.LibLinearLogRegLearner ]
results = Orange.evaluation.testing.cross_validation(learners, data, folds=10)
print "AUC = %.3f" % Orange.evaluation.scoring.AUC(results)[0]
