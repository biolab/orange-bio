import Orange
import obiAssess
import Orange.misc
import obiGeneSets
import obiGene
import numpy
from collections import defaultdict
import stats
import obiGsea
import scipy.stats



def setSig_example_geneset(ex, data):
    """ Gets learning data and example with the same domain, both
    containing only genes from the gene set. """

    distances = [ [], [] ]    

    def pearsonr(v1, v2):
        return numpy.corrcoef([v1, v2])[0,1]

    def pearson(ex1, ex2):
        #leaves undefined elements out

        attrs = range(len(ex1.domain.attributes))
        vals1 = [ ex1[i].value for i in attrs ]
        vals2 = [ ex2[i].value for i in attrs ]

        common = [ True if v1 != "?" and v2 != "?" else False \
            for v1,v2 in zip(vals1,vals2) ]
        vals1 = [ v for v,c in zip(vals1, common) if c ]
        vals2 = [ v for v,c in zip(vals2, common) if c ]

        return numpy.corrcoef([vals1, vals2])[0,1]

    def ttest(ex1, ex2):
        try:
            return stats.lttest_ind(ex1, ex2)[0]
        except:
            return 0.0
    
    #maps class value to its index
    classValueMap = dict( [ (val,i) for i,val in enumerate(data.domain.class_var.values) ])
 
    #create distances to all learning data - save or other class
    for c in data:
        distances[classValueMap[c[-1].value]].append(pearson(c, ex))

    return ttest(distances[0], distances[1])

def mat_ni(data, matcher):
    nm = matcher([at.name for at in data.domain.attributes])
    name_ind = dict((n.name,i) for i,n in enumerate(data.domain.attributes))
    return nm, name_ind

def select_genesets(nm, gene_sets, min_size=3, max_size=1000, min_part=0.1):
    """ Returns a list of gene sets that have sizes in limits """

    def ok_sizes(gs):
        """compares sizes of genesets to limitations"""
        transl = filter(lambda x: x != None, [ nm.umatch(gene) for gene in gs.genes ])
        if len(transl) >= min_size \
            and len(transl) <= max_size \
            and float(len(transl))/len(gs.genes) >= min_part:
            return True
        return False

    return filter(ok_sizes, gene_sets) 

class SetSig(object):

    __new__ = Orange.misc._orange__new__(object)

    def __init__(self, matcher, gene_sets, min_size=3, max_size=1000, min_part=0.1, class_values=None):
        self.matcher = matcher
        self.gene_sets = gene_sets
        self.min_size = min_size
        self.max_size = max_size
        self.min_part = min_part
        self.class_values = class_values
        self._cache = {}

    def _mat_ni(self, data):
        """ With cached gene matchers. """
        if data.domain not in self._cache:
            self._cache[data.domain] = mat_ni(data, self.matcher)
        return self._cache[data.domain]

    def __call__(self, data, weight_id=None):

        data = obiGsea.takeClasses(data, classValues=self.class_values)
        nm,_ =  self._mat_ni(data)
        gene_sets = select_genesets(nm, self.gene_sets, self.min_size, self.max_size, self.min_part)

        attributes = []

        for gs in gene_sets:
            at = Orange.feature.Continuous(name=str(gs))

            def t(ex, w, gs=gs, data=data): #copy od the data
                geneset = list(gs.genes)

                nm, name_ind = self._mat_ni(data)
                nm2, name_ind2 = self._mat_ni(ex)

                genes = [ nm.umatch(gene) for gene in geneset ]
                genes2 = [ nm2.umatch(gene) for gene in geneset ]

                genes, genes2 = zip(*[ (g,g2) for g,g2 in zip(genes, genes2) if g != None])

                domain = Orange.data.Domain([data.domain.attributes[name_ind[gene]] for gene in genes], data.domain.class_var)
                datao = Orange.data.Table(domain, data)

                def vou(ex, gn, indices):
                    """ returns the value or unknown for the given gene name"""
                    if gn not in indices:
                        return "?"
                    else:
                        return ex[indices[gn]].value
                
                #convert the example to the same domain
                exvalues = [ vou(ex, gn, name_ind2) for gn in genes2 ] + [ "?" ]
                example = Orange.data.Instance(domain, exvalues)

                return setSig_example_geneset(example, datao) #only this one is setsig specific
         
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
    matcher = obiGene.matcher([])
    choosen_cv = ["Iris-setosa", "Iris-versicolor"]

    """
    data = Orange.data.Table("DLBCL_200a")
    gsets = obiGeneSets.collections((("KEGG",),"9606"))
    matcher = obiGene.matcher([obiGene.GMKEGG("hsa")])
    choosen_cv = None
    """


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

    ass = SetSig(matcher=matcher, gene_sets=gsets, class_values=choosen_cv, min_part=0.0)
    ass = ass(data)
    ar = to_old_dic(ass.domain, data[:5])
    pp2(ar)
