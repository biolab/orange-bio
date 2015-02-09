import Orange
import orangecontrib.bio.geo

gds = orangecontrib.bio.geo.GDS("GDS2960")
data = gds.getdata(sample_type="disease state", transpose=True, report_genes=True)
print("Samples: %d, Genes: %d" % (len(data), len(data.domain.attributes)))

if Orange.__version__ > "3":
    learners = [ Orange.classification.LogisticRegressionLearner() ]
    results = Orange.evaluation.testing.CrossValidation(data, learners, k=10)
else:
    learners = [ Orange.classification.logreg.LibLinearLogRegLearner() ]
    results = Orange.evaluation.testing.cross_validation(learners, data, folds=10)

print("AUC = %.3f" % Orange.evaluation.scoring.AUC(results)[0])



