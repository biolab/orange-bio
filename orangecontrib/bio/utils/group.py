
from collections import defaultdict
from operator import add
import numpy
import math

def data_type(vals):
    try:
        _ = [ int(a) for a in vals ]
        return int
    except ValueError:
        try:
            _ = [ float(a) for a in vals ]
            return float
        except ValueError:
            return lambda x: x

def separate_by(data, separate, ignore=[], consider=None, add_empty=True):
    """
    data - the data - annotations are saved in the at.attributes
    annotatitions: keys of at.attributes  by which to separate
    ignore: ignore values of these annotations
    consider: consider only these annotations
    """
    ignore = set(ignore)

    annotations = [ at.attributes for at in data.domain.attributes ]

    all_values = defaultdict(set)
    for a in annotations:
        for k,v in a.iteritems():
            all_values[k].add(v)

    types = {}
    for k,vals in all_values.iteritems():
        types[k] = data_type(vals)
    
    groups = defaultdict(list)
    for i,a in enumerate(annotations):
        groups[tuple(a[k] for k in separate)].append(i)

    different_in_all = set(k \
        for k,vals in all_values.iteritems() \
        if len(vals) == len(annotations) or len(vals) == 1)

    other_relevant = set(all_values.keys()) - different_in_all - ignore - set(separate)
    if consider != None:
        other_relevant &= set(consider)
    other_relevant = sorted(other_relevant) #TODO how to order them?

    def relevant_vals(annotation):
        if isinstance(annotation, tuple):
            return annotation
        return tuple(types[v](annotation[v]) for v in other_relevant)

    other_relevant_d2 = defaultdict(int) #"multiset" - number
    #of maximum occurances of a relevant value in a group
    for _,g in groups.items():
        d = defaultdict(int)
        for i in g:
            d[relevant_vals(annotations[i])] += 1
        for rv,n in d.items():
            if n > other_relevant_d2[rv]:
                other_relevant_d2[rv] = n
    
    if add_empty: #fill in with "empty" relevant vals
        ngroups = {}
        for g in groups:
            need_to_fill = other_relevant_d2.copy()
            for i in groups[g]:
                need_to_fill[relevant_vals(annotations[i])] -= 1
            add = []
            for rv,num in need_to_fill.items():
                for a in range(num):
                    add.append(rv)
            ngroups[g] = groups[g] + add
        groups = ngroups

    ngroups = {}
    uniquepos = {} #which positions are unique
    for g in groups:
        elements = list(groups[g])

        rv2 = lambda x: relevant_vals(annotations[x] if isinstance(x,int) else x)

        ngroups[g] = map(lambda x: x if isinstance(x,int) else None,
            sorted(elements, key=rv2))

        d = defaultdict(int) #get groups of different relevant values
        for i in elements:
            d[rv2(i)] += 1
        
        uniquepos[g] = map(lambda x: not d[rv2(x)] > 1,
             sorted(elements, key=rv2))
    
    return ngroups, uniquepos

def float_or_none(value):
    return value.value if value.value != "?" else None

def linearize(data, ids):
    """ Returns a list of floats in the data subspace (or None's
    if the values are unknown or not present. """
    l = [ [ None ] * len(data) if id1 == None \
        else [ float_or_none(ex[id1]) for ex in data ] for id1 in ids ]
    l = reduce(add, l)
    return l

def pearson_lists(l1, l2):
    """ Returns pearson correlation between two lists. Ignores elements
    which are None."""
    okvals = [ (a,b) for a,b in zip(l1,l2) if a != None and b != None ]
    return numpy.corrcoef([ [ v[0] for v in okvals], [ v[1] for v in okvals] ])[0,1]

def euclidean_lists(l1, l2):
    """ Returns pearson correlation between two lists. Ignores elements
    which are None."""
    okvals = [ (a,b) for a,b in zip(l1,l2) if a != None and b != None ]
    return math.sqrt( sum((a-b)*(a-b) for a,b in okvals ))

def spearman_lists(l1, l2):
    """ Returns pearson correlation between two lists. Ignores elements
    which are None."""
    import scipy.stats
    okvals = [ (a,b) for a,b in zip(l1,l2) if a != None and b != None ]
    #print okvals, len(okvals)
    return scipy.stats.spearmanr([ v[0] for v in okvals], [ v[1] for v in okvals] )[0]

def dist_spearman(l1, l2):
    return (1.-spearman_lists(l1, l2))/2

def dist_pcorr(l1, l2):
    #normalized to 0..1
    return (1.-pearson_lists(l1, l2))/2

def dist_eucl(l1, l2):
    return euclidean_lists(l1, l2)
