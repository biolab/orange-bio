import Orange
import obiAssess
import Orange.misc
import obiGeneSets
import obiGene
import numpy
from collections import defaultdict
import statc
import stats

def selectGenesetsData(data, matcher, geneSets, minSize=3, maxSize=1000, minPart=0.1, classValues=None):
    """
    Returns gene sets and data which falling under upper criteria.
    """
    gso = obiGsea.GSEA(data, matcher=matcher, classValues=classValues, atLeast=0)
    gso.addGenesets(geneSets)
    okgenesets = gso.selectGenesets(minSize=minSize, maxSize=maxSize, minPart=minPart).keys()
    gsetsnum = gso.to_gsetsnum(okgenesets)
    return gso.data, okgenesets, gsetsnum

class SetSig(object):

    __new__ = Orange.misc._orange__new__(object)

    def __init__(self, matcher, gene_sets, min_size=3, max_size=1000, min_part=0.1, class_values=None):
        self.matcher = matcher
        self.gene_sets = gene_sets
        self.min_size = min_size
        self.max_size = max_size
        self.min_part = min_part
        self.class_values = class_values

    def __call__(self, data, weight_id=None):
        data, oknames, gsetsnum = obiAssess.selectGenesetsData(data, 
            self.matcher, self.gene_sets,
            minSize=self.min_size, maxSize=self.max_size, 
            minPart=self.min_part, classValues=self.class_values)

        def setSig_example_geneset(ex, data):
            """ ex contains only selected genes """

            distances = [ [], [] ]    

            def pearsonr(v1, v2):
                try:
                    return statc.pearsonr(v1, v2)[0]
                except:
                    return numpy.corrcoef([v1, v2])[0,1]

            def pearson(ex1, ex2):
                attrs = range(len(ex1.domain.attributes))
                vals1 = [ ex1[i].value for i in attrs ]
                vals2 = [ ex2[i].value for i in attrs ]
                return pearsonr(vals1, vals2)

            def ttest(ex1, ex2):
                try:
                    return stats.lttest_ind(ex1, ex2)[0]
                except:
                    return 0.0
            
            #maps class value to its index
            classValueMap = dict( [ (val,i) for i,val in enumerate(data.domain.classVar.values) ])
         
            #create distances to all learning data - save or other class
            for c in data:
                distances[classValueMap[c[-1].value]].append(pearson(c, ex))

            return ttest(distances[0], distances[1])

        attributes = []

        for name, gs in gsetsnum.items(): #for each geneset
            #for each gene set: take the attribute subset and work on the attribute subset only
            #only select the subset of genes from the learning data
            at = Orange.feature.Continuous(name=name.id)

            def t(ex, w, gs=gs, ldata=data):
                domain = Orange.data.Domain([ldata.domain.attributes[ai] for ai in gs], ldata.domain.classVar)
                datao = Orange.data.Table(domain, ldata)
                example = Orange.data.Instance(domain, ex) #domains need to be the same
                return setSig_example_geneset(example, datao)
         
            at.get_value_from = t
            attributes.append(at)
       
        newdomain = Orange.data.Domain(attributes, data.domain.class_var)
        return Orange.data.Table(newdomain, data)

if __name__ == "__main__":

    data = Orange.data.Table("iris")
    gsets = obiGeneSets.collections({
        "ALL": ['sepal length', 'sepal width', 'petal length', 'petal width'],
        "f3": ['sepal length', 'sepal width', 'petal length'],
        "l3": ['sepal width', 'petal length', 'petal width'],
        })

    fp = 120
    ldata = Orange.data.Table(data.domain, data[:fp])
    tdata = Orange.data.Table(data.domain, data[fp:])

    matcher = obiGene.matcher([])

    choosen_cv = ["Iris-setosa", "Iris-versicolor"]
    def to_old_dic(d, data):
        ar = defaultdict(list)
        for ex1 in data:
            ex = d(ex1)
            for a,x in zip(d.attributes, ex):
                ar[a.name].append(x.value)
        return ar

    def pp2(ar):
        ol =  sorted(ar.items())
        print '\n'.join([ a + ": " +str(b) for a,b in ol])

    ass = SetSig(ldata, matcher=matcher, gene_sets=gsets, class_values=choosen_cv, min_part=0.0)
    print ass.domain
    ar = to_old_dic(ass.domain, data[:5])
    pp2(ar)
